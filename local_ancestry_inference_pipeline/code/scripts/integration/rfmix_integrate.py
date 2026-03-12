#!/usr/bin/env python3
"""Standardize RFMix MSP output for downstream benchmark tasks.

This script converts an RFMix `.msp.tsv` file into a long-form interval table
with one row per (window, sample). It also writes ancestry code mappings and,
if a trio file is provided, a trio availability report for Mendelian checks.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Tuple


@dataclass(frozen=True)
class Trio:
    trio_id: str
    child_id: str
    father_id: str
    mother_id: str


@dataclass(frozen=True)
class SampleCols:
    sample_id: str
    hap0_col: int
    hap1_col: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert RFMix MSP output to standardized LAI intervals."
    )
    parser.add_argument(
        "--rfmix-msp",
        required=True,
        help="Path to RFMix MSP file (e.g., test_chr22.msp.tsv).",
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
    return parser.parse_args()


def parse_subpop_map(first_line: str) -> Dict[int, str]:
    prefix = "#Subpopulation order/codes:"
    if not first_line.startswith(prefix):
        raise ValueError(
            "MSP first line does not start with '#Subpopulation order/codes:'."
        )

    tail = first_line.split(":", 1)[1].strip()
    entries = [x.strip() for x in tail.split("\t") if x.strip()]
    code_to_ancestry: Dict[int, str] = {}
    for entry in entries:
        if "=" not in entry:
            continue
        ancestry, code_str = entry.split("=", 1)
        ancestry = ancestry.strip()
        code = int(code_str.strip())
        code_to_ancestry[code] = ancestry

    if not code_to_ancestry:
        raise ValueError("Failed to parse ancestry code map from MSP first line.")
    return code_to_ancestry


def parse_sample_columns(header_line: str) -> List[SampleCols]:
    cols = header_line.rstrip("\n").split("\t")
    if cols and cols[0].startswith("#"):
        cols[0] = cols[0].lstrip("#")

    if len(cols) < 8:
        raise ValueError("MSP header has too few columns to contain haplotype calls.")

    sample_state: Dict[str, Dict[str, int]] = {}
    for idx in range(6, len(cols)):
        name = cols[idx]
        if "." not in name:
            continue
        sample_id, hap = name.rsplit(".", 1)
        if hap not in {"0", "1"}:
            continue
        if sample_id not in sample_state:
            sample_state[sample_id] = {"first_idx": idx}
        sample_state[sample_id][hap] = idx
        if idx < sample_state[sample_id]["first_idx"]:
            sample_state[sample_id]["first_idx"] = idx

    result: List[SampleCols] = []
    for sample_id, state in sample_state.items():
        if "0" not in state or "1" not in state:
            continue
        result.append(
            SampleCols(sample_id=sample_id, hap0_col=state["0"], hap1_col=state["1"])
        )

    result.sort(key=lambda x: sample_state[x.sample_id]["first_idx"])
    if not result:
        raise ValueError("No paired sample haplotype columns found in MSP header.")
    return result


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


def convert_msp(
    msp_path: str,
    segments_path: str,
    code_map: Dict[int, str],
    sample_cols: List[SampleCols],
    member_lookup: Dict[str, Tuple[str, str, str, str, str]],
) -> Tuple[int, int, set[str]]:
    out_rows = 0
    window_count = 0
    sample_ids = {s.sample_id for s in sample_cols}

    with open(msp_path, "r", encoding="utf-8") as in_handle, open(
        segments_path, "w", encoding="utf-8", newline=""
    ) as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(
            [
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
        )

        # Skip first two header lines.
        next(in_handle)
        next(in_handle)

        reader = csv.reader(in_handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if row[0].startswith("#"):
                continue

            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            n_snps = int(row[5])
            window_count += 1

            for sc in sample_cols:
                h1_code = int(row[sc.hap0_col])
                h2_code = int(row[sc.hap1_col])
                h1_anc = code_map.get(h1_code, f"UNK_{h1_code}")
                h2_anc = code_map.get(h2_code, f"UNK_{h2_code}")

                trio_id = ""
                child_id = ""
                father_id = ""
                mother_id = ""
                role = ""
                if sc.sample_id in member_lookup:
                    trio_id, child_id, father_id, mother_id, role = member_lookup[
                        sc.sample_id
                    ]

                writer.writerow(
                    [
                        "RFMix",
                        chrom,
                        start,
                        end,
                        sc.sample_id,
                        h1_anc,
                        h2_anc,
                        h1_code,
                        h2_code,
                        n_snps,
                        os.path.abspath(msp_path),
                        trio_id,
                        child_id,
                        father_id,
                        mother_id,
                        role,
                    ]
                )
                out_rows += 1

    return window_count, out_rows, sample_ids


def main() -> int:
    args = parse_args()
    msp_path = os.path.abspath(args.rfmix_msp)
    trio_path = os.path.abspath(args.trio_tsv) if args.trio_tsv else ""
    out_prefix = os.path.abspath(args.out_prefix)

    if not os.path.isfile(msp_path):
        print(f"ERROR: MSP file not found: {msp_path}", file=sys.stderr)
        return 1
    if trio_path and not os.path.isfile(trio_path):
        print(f"ERROR: Trio file not found: {trio_path}", file=sys.stderr)
        return 1

    out_dir = os.path.dirname(out_prefix) or "."
    os.makedirs(out_dir, exist_ok=True)

    with open(msp_path, "r", encoding="utf-8") as handle:
        first_line = handle.readline().rstrip("\n")
        second_line = handle.readline().rstrip("\n")

    code_map = parse_subpop_map(first_line)
    sample_cols = parse_sample_columns(second_line)

    trios: List[Trio] = []
    member_lookup: Dict[str, Tuple[str, str, str, str, str]] = {}
    if trio_path:
        trios = load_trios(trio_path)
        member_lookup = build_member_lookup(trios)

    segments_path = f"{out_prefix}.segments.tsv"
    codes_path = f"{out_prefix}.ancestry_codes.tsv"
    trio_availability_path = f"{out_prefix}.trio_availability.tsv"

    window_count, out_rows, samples_in_output = convert_msp(
        msp_path=msp_path,
        segments_path=segments_path,
        code_map=code_map,
        sample_cols=sample_cols,
        member_lookup=member_lookup,
    )
    write_code_map(codes_path, code_map)

    if trios:
        write_trio_availability(trio_availability_path, trios, samples_in_output)

    print("RFMix integration complete")
    print(f"  MSP input: {msp_path}")
    print(f"  Windows parsed: {window_count}")
    print(f"  Samples in MSP: {len(sample_cols)}")
    print(f"  Output rows: {out_rows}")
    print(f"  Segments table: {segments_path}")
    print(f"  Ancestry code map: {codes_path}")
    if trios:
        print(f"  Trio availability: {trio_availability_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
