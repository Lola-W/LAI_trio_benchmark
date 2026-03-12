#!/usr/bin/env python3
"""Standardize FLARE local ancestry output for downstream benchmark tasks.

This script converts FLARE `.anc.vcf.gz` output into a long-form interval table
with one row per (segment, sample). Segments are built by run-length encoding
consecutive markers with identical (AN1, AN2) ancestry states.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import gzip
import os
import re
import sys
from dataclasses import dataclass
from datetime import datetime
from tempfile import TemporaryDirectory
from typing import Dict, Iterable, List, Optional, Tuple


@dataclass(frozen=True)
class Trio:
    trio_id: str
    child_id: str
    father_id: str
    mother_id: str


@dataclass
class RunState:
    start: int
    last_pos: int
    h1_code: int
    h2_code: int
    snp_count: int


SEGMENT_HEADER = [
    "Tool",
    "Chromosome",
    "Start",
    "End",
    "Sample_ID",
    "Haplotype_1_Ancestry",
    "Haplotype_2_Ancestry",
    "Haplotype_1_Code",
    "Haplotype_2_Code",
    "Window_SNP_Count",
    "Source_File",
    "Trio_ID",
    "Child_ID",
    "Father_ID",
    "Mother_ID",
    "Sample_Role",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert FLARE anc.vcf.gz output to standardized LAI intervals."
    )

    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument(
        "--flare-dir",
        default="",
        help=(
            "Directory to recursively scan for '*.anc.vcf.gz'. "
            "Each file must have sibling '.model' file with matching prefix."
        ),
    )
    source_group.add_argument(
        "--anc-vcf",
        default="",
        help="Single FLARE anc.vcf.gz file.",
    )

    parser.add_argument(
        "--model",
        default="",
        help="Model file for --anc-vcf mode (default: infer from anc filename).",
    )
    parser.add_argument(
        "--trio-tsv",
        default="",
        help="Optional trio file with 4 columns: trio_id child father mother.",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        help=(
            "Output prefix path. Script writes: "
            "<prefix>.segments.tsv, <prefix>.ancestry_codes.tsv, "
            "and optionally <prefix>.trio_availability.tsv."
        ),
    )
    parser.add_argument(
        "--terminal-end",
        choices=["last_pos_plus_one", "contig_plus_one"],
        default="last_pos_plus_one",
        help=(
            "How to set segment end for the last marker on a chromosome. "
            "default: last_pos_plus_one"
        ),
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="Optional cap on number of anc files to process (0 means no cap).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help=(
            "Number of worker processes for file-level parallelism "
            "(default: 1)."
        ),
    )
    return parser.parse_args()


def load_trios(path: str) -> List[Trio]:
    trios: List[Trio] = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(
                    f"Invalid trio row (expected >=4 fields): {raw.rstrip()}"
                )
            trios.append(
                Trio(
                    trio_id=fields[0],
                    child_id=fields[1],
                    father_id=fields[2],
                    mother_id=fields[3],
                )
            )
    return trios


def build_member_lookup(
    trios: List[Trio],
) -> Dict[str, Tuple[str, str, str, str, str]]:
    """Return sample_id -> (trio_id, child_id, father_id, mother_id, role)."""
    lookup: Dict[str, Tuple[str, str, str, str, str]] = {}
    for t in trios:
        values = {
            t.child_id: (t.trio_id, t.child_id, t.father_id, t.mother_id, "child"),
            t.father_id: (t.trio_id, t.child_id, t.father_id, t.mother_id, "father"),
            t.mother_id: (t.trio_id, t.child_id, t.father_id, t.mother_id, "mother"),
        }
        for sid, payload in values.items():
            if sid not in lookup:
                lookup[sid] = payload
    return lookup


def discover_pairs(
    flare_dir: str,
    anc_vcf: str,
    model: str,
    max_files: int,
) -> List[Tuple[str, str]]:
    pairs: List[Tuple[str, str]] = []

    if anc_vcf:
        anc_path = os.path.abspath(anc_vcf)
        if not os.path.isfile(anc_path):
            raise FileNotFoundError(f"FLARE anc file not found: {anc_path}")
        model_path = os.path.abspath(model) if model else infer_model_path(anc_path)
        if not os.path.isfile(model_path):
            raise FileNotFoundError(
                f"Model file not found for anc file: {anc_path}; expected {model_path}"
            )
        pairs.append((anc_path, model_path))
        return pairs

    root = os.path.abspath(flare_dir)
    if not os.path.isdir(root):
        raise FileNotFoundError(f"FLARE directory not found: {root}")

    for dirpath, _, filenames in os.walk(root):
        for fname in filenames:
            if not fname.endswith(".anc.vcf.gz"):
                continue
            anc_path = os.path.join(dirpath, fname)
            model_path = infer_model_path(anc_path)
            if not os.path.isfile(model_path):
                raise FileNotFoundError(
                    f"Missing model file for {anc_path}; expected {model_path}"
                )
            pairs.append((anc_path, model_path))

    pairs.sort(key=lambda x: x[0])
    if max_files > 0:
        pairs = pairs[:max_files]

    if not pairs:
        raise FileNotFoundError(f"No '*.anc.vcf.gz' files found under: {root}")
    return pairs


def infer_model_path(anc_path: str) -> str:
    if not anc_path.endswith(".anc.vcf.gz"):
        raise ValueError(f"Unexpected anc filename: {anc_path}")
    return anc_path[: -len(".anc.vcf.gz")] + ".model"


def infer_kid_label(anc_path: str) -> str:
    """Infer child/sample label from path for progress logging."""
    parent = os.path.basename(os.path.dirname(anc_path))
    if parent and parent not in {"", ".", ".."}:
        return parent
    base = os.path.basename(anc_path)
    if base.endswith(".anc.vcf.gz"):
        return base[: -len(".anc.vcf.gz")]
    return base


def log_progress(message: str) -> None:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {message}", file=sys.stderr, flush=True)


def parse_model_codes(model_path: str) -> Dict[int, str]:
    with open(model_path, "r", encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle]

    saw_header = False
    for line in lines:
        text = line.strip()
        if not text:
            continue
        if text.startswith("# list of ancestries"):
            saw_header = True
            continue
        if not saw_header:
            continue
        if text.startswith("#"):
            continue
        ancestries = text.split()
        if not ancestries:
            continue
        return {idx: anc for idx, anc in enumerate(ancestries)}

    raise ValueError(f"Failed to parse ancestry list from model file: {model_path}")


def write_code_map(path: str, code_map: Dict[int, str]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Code", "Ancestry"])
        for code in sorted(code_map):
            writer.writerow([code, code_map[code]])


def write_trio_availability(
    path: str, trios: List[Trio], samples_in_output: set[str]
) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "Trio_ID",
                "Child_ID",
                "Father_ID",
                "Mother_ID",
                "Child_In_Output",
                "Father_In_Output",
                "Mother_In_Output",
                "Any_Missing",
            ]
        )
        for t in trios:
            child_hit = int(t.child_id in samples_in_output)
            father_hit = int(t.father_id in samples_in_output)
            mother_hit = int(t.mother_id in samples_in_output)
            any_missing = int(not (child_hit and father_hit and mother_hit))
            writer.writerow(
                [
                    t.trio_id,
                    t.child_id,
                    t.father_id,
                    t.mother_id,
                    child_hit,
                    father_hit,
                    mother_hit,
                    any_missing,
                ]
            )


def parse_int_token(token: str) -> Optional[int]:
    token = token.strip()
    if token in {"", "."}:
        return None
    if "," in token:
        token = token.split(",", 1)[0]
    try:
        return int(token)
    except ValueError:
        return None


def terminal_end(
    mode: str,
    chrom: str,
    last_pos: int,
    contig_lens: Dict[str, int],
) -> int:
    if mode == "contig_plus_one" and chrom in contig_lens:
        return contig_lens[chrom] + 1
    return last_pos + 1


def flush_runs(
    writer: csv.writer,
    runs: Dict[str, RunState],
    chrom: str,
    end_pos: int,
    anc_path: str,
    code_map: Dict[int, str],
    member_lookup: Dict[str, Tuple[str, str, str, str, str]],
) -> int:
    written = 0
    for sample_id, run in runs.items():
        if end_pos <= run.start:
            out_end = run.start + 1
        else:
            out_end = end_pos

        h1_anc = code_map.get(run.h1_code, f"UNK_{run.h1_code}")
        h2_anc = code_map.get(run.h2_code, f"UNK_{run.h2_code}")

        trio_id = ""
        child_id = ""
        father_id = ""
        mother_id = ""
        role = ""
        if sample_id in member_lookup:
            trio_id, child_id, father_id, mother_id, role = member_lookup[sample_id]

        writer.writerow(
            [
                "FLARE",
                chrom,
                run.start,
                out_end,
                sample_id,
                h1_anc,
                h2_anc,
                run.h1_code,
                run.h2_code,
                run.snp_count,
                os.path.abspath(anc_path),
                trio_id,
                child_id,
                father_id,
                mother_id,
                role,
            ]
        )
        written += 1
    runs.clear()
    return written


def process_anc_file(
    anc_path: str,
    code_map: Dict[int, str],
    member_lookup: Dict[str, Tuple[str, str, str, str, str]],
    writer: csv.writer,
    terminal_mode: str,
) -> Tuple[int, int, set[str]]:
    contig_lens: Dict[str, int] = {}
    sample_ids: List[str] = []
    samples_seen: set[str] = set()

    current_chrom = ""
    last_pos_in_chrom = -1
    runs: Dict[str, RunState] = {}

    rows_written = 0
    marker_count = 0

    contig_re = re.compile(r"##contig=<ID=([^,>]+),length=(\d+)")

    with gzip.open(anc_path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line:
                continue

            if line.startswith("##contig="):
                m = contig_re.search(line)
                if m:
                    contig_lens[m.group(1)] = int(m.group(2))
                continue

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                fields = line.split("\t")
                if len(fields) < 10:
                    raise ValueError(
                        f"No sample columns found in FLARE anc file: {anc_path}"
                    )
                sample_ids = fields[9:]
                samples_seen.update(sample_ids)
                continue

            fields = line.split("\t")
            if len(fields) < 10:
                continue

            chrom = fields[0]
            pos = int(fields[1])

            if not current_chrom:
                current_chrom = chrom
            elif chrom != current_chrom:
                end_pos = terminal_end(
                    mode=terminal_mode,
                    chrom=current_chrom,
                    last_pos=last_pos_in_chrom,
                    contig_lens=contig_lens,
                )
                rows_written += flush_runs(
                    writer=writer,
                    runs=runs,
                    chrom=current_chrom,
                    end_pos=end_pos,
                    anc_path=anc_path,
                    code_map=code_map,
                    member_lookup=member_lookup,
                )
                current_chrom = chrom

            fmt = fields[8].split(":")
            try:
                an1_idx = fmt.index("AN1")
                an2_idx = fmt.index("AN2")
            except ValueError as exc:
                raise ValueError(
                    f"AN1/AN2 not found in FORMAT column for {anc_path}"
                ) from exc

            sample_fields = fields[9:]
            if len(sample_fields) != len(sample_ids):
                raise ValueError(
                    f"Sample field count mismatch at {anc_path}:{chrom}:{pos}"
                )

            marker_count += 1
            for i, sample_id in enumerate(sample_ids):
                vals = sample_fields[i].split(":")
                if an1_idx >= len(vals) or an2_idx >= len(vals):
                    h1_code = None
                    h2_code = None
                else:
                    h1_code = parse_int_token(vals[an1_idx])
                    h2_code = parse_int_token(vals[an2_idx])

                if h1_code is None or h2_code is None:
                    if sample_id in runs:
                        rows_written += flush_runs(
                            writer=writer,
                            runs={sample_id: runs[sample_id]},
                            chrom=chrom,
                            end_pos=pos,
                            anc_path=anc_path,
                            code_map=code_map,
                            member_lookup=member_lookup,
                        )
                        runs.pop(sample_id, None)
                    continue

                if sample_id not in runs:
                    runs[sample_id] = RunState(
                        start=pos,
                        last_pos=pos,
                        h1_code=h1_code,
                        h2_code=h2_code,
                        snp_count=1,
                    )
                    continue

                run = runs[sample_id]
                if run.h1_code == h1_code and run.h2_code == h2_code:
                    run.last_pos = pos
                    run.snp_count += 1
                else:
                    rows_written += flush_runs(
                        writer=writer,
                        runs={sample_id: run},
                        chrom=chrom,
                        end_pos=pos,
                        anc_path=anc_path,
                        code_map=code_map,
                        member_lookup=member_lookup,
                    )
                    runs[sample_id] = RunState(
                        start=pos,
                        last_pos=pos,
                        h1_code=h1_code,
                        h2_code=h2_code,
                        snp_count=1,
                    )

            last_pos_in_chrom = pos

    if current_chrom and last_pos_in_chrom >= 0 and runs:
        end_pos = terminal_end(
            mode=terminal_mode,
            chrom=current_chrom,
            last_pos=last_pos_in_chrom,
            contig_lens=contig_lens,
        )
        rows_written += flush_runs(
            writer=writer,
            runs=runs,
            chrom=current_chrom,
            end_pos=end_pos,
            anc_path=anc_path,
            code_map=code_map,
            member_lookup=member_lookup,
        )

    return marker_count, rows_written, samples_seen


def ensure_all_maps_match(
    pairs: Iterable[Tuple[str, str]],
) -> Dict[int, str]:
    canonical: Optional[Dict[int, str]] = None
    canonical_source = ""

    for anc_path, model_path in pairs:
        current = parse_model_codes(model_path)
        if canonical is None:
            canonical = current
            canonical_source = model_path
            continue

        if current != canonical:
            raise ValueError(
                "FLARE model ancestry ordering differs across files. "
                f"Canonical: {canonical_source}; mismatch: {model_path}"
            )

    if canonical is None:
        raise ValueError("No FLARE model map found.")
    return canonical


def process_anc_to_part_file(
    anc_path: str,
    code_map: Dict[int, str],
    member_lookup: Dict[str, Tuple[str, str, str, str, str]],
    terminal_mode: str,
    part_path: str,
) -> Tuple[int, int, List[str]]:
    """Worker helper: process one anc file and write one segment part file."""
    with open(part_path, "w", encoding="utf-8", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(SEGMENT_HEADER)
        markers, rows, samples = process_anc_file(
            anc_path=anc_path,
            code_map=code_map,
            member_lookup=member_lookup,
            writer=writer,
            terminal_mode=terminal_mode,
        )
    return markers, rows, sorted(samples)


def merge_part_files(segments_path: str, part_paths_in_order: List[str]) -> None:
    with open(segments_path, "w", encoding="utf-8", newline="") as out_handle:
        out_writer = csv.writer(out_handle, delimiter="\t")
        out_writer.writerow(SEGMENT_HEADER)
        for part_path in part_paths_in_order:
            with open(part_path, "r", encoding="utf-8") as in_handle:
                # Skip per-part header
                next(in_handle, None)
                for line in in_handle:
                    out_handle.write(line)


def main() -> int:
    args = parse_args()

    trio_path = os.path.abspath(args.trio_tsv) if args.trio_tsv else ""
    out_prefix = os.path.abspath(args.out_prefix)
    if args.threads < 1:
        print("ERROR: --threads must be >= 1", file=sys.stderr)
        return 1

    if trio_path and not os.path.isfile(trio_path):
        print(f"ERROR: Trio file not found: {trio_path}", file=sys.stderr)
        return 1

    out_dir = os.path.dirname(out_prefix) or "."
    os.makedirs(out_dir, exist_ok=True)

    try:
        pairs = discover_pairs(
            flare_dir=args.flare_dir,
            anc_vcf=args.anc_vcf,
            model=args.model,
            max_files=args.max_files,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    try:
        code_map = ensure_all_maps_match(pairs)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    trios: List[Trio] = []
    member_lookup: Dict[str, Tuple[str, str, str, str, str]] = {}
    if trio_path:
        trios = load_trios(trio_path)
        member_lookup = build_member_lookup(trios)

    segments_path = f"{out_prefix}.segments.tsv"
    codes_path = f"{out_prefix}.ancestry_codes.tsv"
    trio_availability_path = f"{out_prefix}.trio_availability.tsv"

    total_markers = 0
    total_rows = 0
    all_samples: set[str] = set()
    effective_threads = min(args.threads, len(pairs))
    try:
        if effective_threads == 1:
            with open(segments_path, "w", encoding="utf-8", newline="") as out_handle:
                writer = csv.writer(out_handle, delimiter="\t")
                writer.writerow(SEGMENT_HEADER)
                for anc_path, _ in pairs:
                    markers, rows, samples = process_anc_file(
                        anc_path=anc_path,
                        code_map=code_map,
                        member_lookup=member_lookup,
                        writer=writer,
                        terminal_mode=args.terminal_end,
                    )
                    total_markers += markers
                    total_rows += rows
                    all_samples.update(samples)
        else:
            with TemporaryDirectory(prefix="flare_integrate_") as tmpdir:
                futures = {}
                with ProcessPoolExecutor(max_workers=effective_threads) as pool:
                    for idx, (anc_path, _) in enumerate(pairs):
                        part_path = os.path.join(tmpdir, f"segments_part_{idx:06d}.tsv")
                        fut = pool.submit(
                            process_anc_to_part_file,
                            anc_path,
                            code_map,
                            member_lookup,
                            args.terminal_end,
                            part_path,
                        )
                        futures[fut] = (idx, anc_path, part_path)

                    results_by_idx = {}
                    for fut in as_completed(futures):
                        idx, anc_path, part_path = futures[fut]
                        try:
                            markers, rows, samples = fut.result()
                        except Exception as exc:
                            raise RuntimeError(
                                f"Failed processing anc file: {anc_path}"
                            ) from exc
                        results_by_idx[idx] = (markers, rows, set(samples), part_path)

                part_paths: List[str] = []
                for idx in range(len(pairs)):
                    markers, rows, samples, part_path = results_by_idx[idx]
                    total_markers += markers
                    total_rows += rows
                    all_samples.update(samples)
                    part_paths.append(part_path)

                merge_part_files(segments_path=segments_path, part_paths_in_order=part_paths)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    write_code_map(codes_path, code_map)
    if trios:
        write_trio_availability(trio_availability_path, trios, all_samples)

    print("FLARE integration complete")
    print(f"  Anc files processed: {len(pairs)}")
    print(f"  Threads used: {effective_threads}")
    print(f"  Markers parsed: {total_markers}")
    print(f"  Samples observed: {len(all_samples)}")
    print(f"  Output rows: {total_rows}")
    print(f"  Segments table: {segments_path}")
    print(f"  Ancestry code map: {codes_path}")
    if trios:
        print(f"  Trio availability: {trio_availability_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
