#!/usr/bin/env python3
"""Merge standardized RFMix/FLARE LAI outputs and prebuild task 2/3 inputs.

Inputs are `*.segments.tsv` produced by:
- rfmix_integrate.py
- flare_integrate.py

Outputs:
- unified segments table
- task-2 Mendelian atomic intervals + summary
- task-3 pairwise cross-tool atomic intervals + summary
"""

from __future__ import annotations

import argparse
import csv
import itertools
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple


SEGMENT_COLUMNS = [
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


@dataclass(frozen=True)
class SegRec:
    tool: str
    chrom: str
    start: int
    end: int
    sample_id: str
    h1_anc: str
    h2_anc: str
    h1_code: Optional[int]
    h2_code: Optional[int]
    window_snp_count: Optional[int]
    source_file: str
    trio_id: str
    child_id: str
    father_id: str
    mother_id: str
    sample_role: str
    input_file: str


@dataclass
class Cursor:
    idx: int = 0


@dataclass(frozen=True)
class TrioMember:
    trio_id: str
    child_id: str
    father_id: str
    mother_id: str
    role: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge standardized RFMix/FLARE outputs and build task 2/3 inputs."
    )
    parser.add_argument(
        "--rfmix-segments",
        nargs="+",
        default=[],
        help="One or more RFMix segment TSV files from rfmix_integrate.py.",
    )
    parser.add_argument(
        "--rfmix-dir",
        default="",
        help=(
            "Directory of pre-standardized RFMix segment TSVs "
            "(loads *.segments.tsv, e.g. data/RFMix/out)."
        ),
    )
    parser.add_argument(
        "--rfmix-trio-tsv",
        default="",
        help=(
            "Optional trio TSV used to backfill trio metadata onto RFMix rows "
            "that already have segment calls but blank trio columns."
        ),
    )
    parser.add_argument(
        "--flare-segments",
        nargs="+",
        default=[],
        help="One or more FLARE segment TSV files from flare_integrate.py.",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix path for merged and task2/task3 files.",
    )
    return parser.parse_args()


def parse_int_or_none(value: str) -> Optional[int]:
    text = value.strip()
    if text in {"", ".", "NA", "nan", "None"}:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def natural_sort_key(text: str) -> Tuple[object, ...]:
    parts = re.split(r"(\d+)", text)
    key: List[object] = []
    for part in parts:
        if not part:
            continue
        key.append(int(part) if part.isdigit() else part)
    return tuple(key)


def discover_segment_files(directory: str) -> List[str]:
    path = os.path.abspath(directory)
    if not os.path.isdir(path):
        raise FileNotFoundError(f"RFMix directory not found: {path}")

    segment_paths = [
        os.path.join(path, name)
        for name in os.listdir(path)
        if name.endswith(".segments.tsv")
    ]
    segment_paths.sort(key=lambda p: natural_sort_key(os.path.basename(p)))
    if not segment_paths:
        raise FileNotFoundError(f"No *.segments.tsv files found in: {path}")
    return segment_paths


def load_trio_lookup(path: str) -> Dict[str, TrioMember]:
    trio_path = os.path.abspath(path)
    if not os.path.isfile(trio_path):
        raise FileNotFoundError(f"Trio file not found: {trio_path}")

    lookup: Dict[str, TrioMember] = {}
    with open(trio_path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(
                    f"Invalid trio row in {trio_path} (expected >=4 fields): {raw.rstrip()}"
                )

            trio_id, child_id, father_id, mother_id = fields[:4]
            members = {
                child_id: TrioMember(
                    trio_id=trio_id,
                    child_id=child_id,
                    father_id=father_id,
                    mother_id=mother_id,
                    role="child",
                ),
                father_id: TrioMember(
                    trio_id=trio_id,
                    child_id=child_id,
                    father_id=father_id,
                    mother_id=mother_id,
                    role="father",
                ),
                mother_id: TrioMember(
                    trio_id=trio_id,
                    child_id=child_id,
                    father_id=father_id,
                    mother_id=mother_id,
                    role="mother",
                ),
            }
            for sample_id, member in members.items():
                lookup.setdefault(sample_id, member)

    return lookup


def load_segments(
    paths: Iterable[str],
    rfmix_trio_lookup: Optional[Dict[str, TrioMember]] = None,
) -> Tuple[List[SegRec], int]:
    records: List[SegRec] = []
    backfilled_rows = 0
    for raw_path in paths:
        path = os.path.abspath(raw_path)
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Segment file not found: {path}")

        with open(path, "r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = [c for c in SEGMENT_COLUMNS if c not in reader.fieldnames]
            if missing:
                raise ValueError(
                    f"Missing columns in {path}: {', '.join(missing)}"
                )

            for row in reader:
                start = int(row["Start"])
                end = int(row["End"])
                if end <= start:
                    continue

                tool = row["Tool"].strip()
                sample_id = row["Sample_ID"].strip()
                trio_id = row["Trio_ID"].strip()
                child_id = row["Child_ID"].strip()
                father_id = row["Father_ID"].strip()
                mother_id = row["Mother_ID"].strip()
                sample_role = row["Sample_Role"].strip()

                if (
                    tool == "RFMix"
                    and rfmix_trio_lookup
                    and sample_id in rfmix_trio_lookup
                ):
                    member = rfmix_trio_lookup[sample_id]
                    updated = False
                    if not trio_id:
                        trio_id = member.trio_id
                        updated = True
                    if not child_id:
                        child_id = member.child_id
                        updated = True
                    if not father_id:
                        father_id = member.father_id
                        updated = True
                    if not mother_id:
                        mother_id = member.mother_id
                        updated = True
                    if not sample_role:
                        sample_role = member.role
                        updated = True
                    if updated:
                        backfilled_rows += 1

                records.append(
                    SegRec(
                        tool=tool,
                        chrom=row["Chromosome"].strip(),
                        start=start,
                        end=end,
                        sample_id=sample_id,
                        h1_anc=row["Haplotype_1_Ancestry"].strip(),
                        h2_anc=row["Haplotype_2_Ancestry"].strip(),
                        h1_code=parse_int_or_none(row["Haplotype_1_Code"]),
                        h2_code=parse_int_or_none(row["Haplotype_2_Code"]),
                        window_snp_count=parse_int_or_none(row["Window_SNP_Count"]),
                        source_file=row["Source_File"].strip(),
                        trio_id=trio_id,
                        child_id=child_id,
                        father_id=father_id,
                        mother_id=mother_id,
                        sample_role=sample_role,
                        input_file=path,
                    )
                )
    return records, backfilled_rows


def sort_segments(segments: List[SegRec]) -> None:
    segments.sort(
        key=lambda r: (
            r.tool,
            r.sample_id,
            r.chrom,
            r.start,
            r.end,
            r.h1_anc,
            r.h2_anc,
        )
    )


def write_unified(path: str, segments: List[SegRec]) -> None:
    header = SEGMENT_COLUMNS + ["Input_File"]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        for r in segments:
            writer.writerow(
                [
                    r.tool,
                    r.chrom,
                    r.start,
                    r.end,
                    r.sample_id,
                    r.h1_anc,
                    r.h2_anc,
                    "" if r.h1_code is None else r.h1_code,
                    "" if r.h2_code is None else r.h2_code,
                    "" if r.window_snp_count is None else r.window_snp_count,
                    r.source_file,
                    r.trio_id,
                    r.child_id,
                    r.father_id,
                    r.mother_id,
                    r.sample_role,
                    r.input_file,
                ]
            )


def state_at(segments: List[SegRec], cursor: Cursor, pos: int) -> Optional[SegRec]:
    n = len(segments)
    i = cursor.idx
    while i < n and segments[i].end <= pos:
        i += 1
    cursor.idx = i
    if i < n and segments[i].start <= pos < segments[i].end:
        return segments[i]
    return None


def build_task2(
    segments: List[SegRec],
    intervals_path: str,
    summary_path: str,
) -> Tuple[int, int]:
    # Group by tool + trio + chrom and split by role.
    grouped: Dict[Tuple[str, str, str, str, str, str], Dict[str, List[SegRec]]] = {}
    for r in segments:
        if not r.trio_id:
            continue
        if r.sample_role not in {"child", "father", "mother"}:
            continue
        key = (r.tool, r.trio_id, r.chrom, r.child_id, r.father_id, r.mother_id)
        if key not in grouped:
            grouped[key] = {"child": [], "father": [], "mother": []}
        grouped[key][r.sample_role].append(r)

    for key in grouped:
        for role in ("child", "father", "mother"):
            grouped[key][role].sort(key=lambda x: (x.start, x.end))

    summary_stats: Dict[Tuple[str, str, str], Dict[str, int]] = defaultdict(
        lambda: {
            "intervals": 0,
            "bp_total": 0,
            "bp_evaluable": 0,
            "bp_violation_ordered": 0,
            "bp_violation_unordered": 0,
        }
    )

    written_rows = 0
    with open(intervals_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "Tool",
                "Trio_ID",
                "Chromosome",
                "Start",
                "End",
                "BP_Length",
                "Child_ID",
                "Father_ID",
                "Mother_ID",
                "Child_H1_Ancestry",
                "Child_H2_Ancestry",
                "Father_H1_Ancestry",
                "Father_H2_Ancestry",
                "Mother_H1_Ancestry",
                "Mother_H2_Ancestry",
                "Evaluable",
                "Violation_OrderedParents",
                "Violation_UnorderedParents",
            ]
        )

        for (tool, trio_id, chrom, child_id, father_id, mother_id), by_role in grouped.items():
            role_segments = by_role
            bounds = set()
            for role in ("child", "father", "mother"):
                for seg in role_segments[role]:
                    bounds.add(seg.start)
                    bounds.add(seg.end)

            if len(bounds) < 2:
                continue
            b = sorted(bounds)
            cursors = {"child": Cursor(), "father": Cursor(), "mother": Cursor()}

            for i in range(len(b) - 1):
                start = b[i]
                end = b[i + 1]
                if end <= start:
                    continue
                bp = end - start

                child = state_at(role_segments["child"], cursors["child"], start)
                father = state_at(role_segments["father"], cursors["father"], start)
                mother = state_at(role_segments["mother"], cursors["mother"], start)

                evaluable = int(child is not None and father is not None and mother is not None)
                violation_ordered = ""
                violation_unordered = ""

                ch1 = child.h1_anc if child else ""
                ch2 = child.h2_anc if child else ""
                fa1 = father.h1_anc if father else ""
                fa2 = father.h2_anc if father else ""
                mo1 = mother.h1_anc if mother else ""
                mo2 = mother.h2_anc if mother else ""

                if evaluable:
                    ordered_ok = (ch1 in {fa1, fa2}) and (ch2 in {mo1, mo2})
                    swapped_ok = (ch1 in {mo1, mo2}) and (ch2 in {fa1, fa2})
                    violation_ordered = int(not ordered_ok)
                    violation_unordered = int(not (ordered_ok or swapped_ok))

                writer.writerow(
                    [
                        tool,
                        trio_id,
                        chrom,
                        start,
                        end,
                        bp,
                        child_id,
                        father_id,
                        mother_id,
                        ch1,
                        ch2,
                        fa1,
                        fa2,
                        mo1,
                        mo2,
                        evaluable,
                        violation_ordered,
                        violation_unordered,
                    ]
                )
                written_rows += 1

                stat_key = (tool, trio_id, child_id)
                stats = summary_stats[stat_key]
                stats["intervals"] += 1
                stats["bp_total"] += bp
                if evaluable:
                    stats["bp_evaluable"] += bp
                    stats["bp_violation_ordered"] += bp * int(violation_ordered)
                    stats["bp_violation_unordered"] += bp * int(violation_unordered)

    with open(summary_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "Tool",
                "Trio_ID",
                "Child_ID",
                "Intervals",
                "BP_Total",
                "BP_Evaluable",
                "BP_Violation_OrderedParents",
                "BP_Violation_UnorderedParents",
                "ViolationRate_OrderedParents",
                "ViolationRate_UnorderedParents",
            ]
        )
        for (tool, trio_id, child_id), st in sorted(summary_stats.items()):
            bp_eval = st["bp_evaluable"]
            rate_ord = "" if bp_eval == 0 else st["bp_violation_ordered"] / bp_eval
            rate_unord = "" if bp_eval == 0 else st["bp_violation_unordered"] / bp_eval
            writer.writerow(
                [
                    tool,
                    trio_id,
                    child_id,
                    st["intervals"],
                    st["bp_total"],
                    st["bp_evaluable"],
                    st["bp_violation_ordered"],
                    st["bp_violation_unordered"],
                    rate_ord,
                    rate_unord,
                ]
            )

    return len(summary_stats), written_rows


def build_task3(
    segments: List[SegRec],
    intervals_path: str,
    summary_path: str,
) -> Tuple[int, int]:
    # sample/chrom/tool -> segments
    sc_tool: Dict[Tuple[str, str], Dict[str, List[SegRec]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for r in segments:
        sc_tool[(r.sample_id, r.chrom)][r.tool].append(r)

    for key in sc_tool:
        for tool in sc_tool[key]:
            sc_tool[key][tool].sort(key=lambda x: (x.start, x.end))

    summary_stats: Dict[Tuple[str, str, str, str], Dict[str, int]] = defaultdict(
        lambda: {
            "intervals": 0,
            "bp_total": 0,
            "bp_evaluable": 0,
            "bp_match_ordered": 0,
            "bp_match_unordered": 0,
        }
    )

    written_rows = 0
    with open(intervals_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "Sample_ID",
                "Chromosome",
                "Tool_A",
                "Tool_B",
                "Start",
                "End",
                "BP_Length",
                "A_H1_Ancestry",
                "A_H2_Ancestry",
                "B_H1_Ancestry",
                "B_H2_Ancestry",
                "Evaluable",
                "Match_OrderedHaps",
                "Match_UnorderedHaps",
            ]
        )

        for (sample_id, chrom), tool_map in sorted(sc_tool.items()):
            tools = sorted(tool_map.keys())
            if len(tools) < 2:
                continue

            for tool_a, tool_b in itertools.combinations(tools, 2):
                seg_a = tool_map[tool_a]
                seg_b = tool_map[tool_b]
                bounds = set()
                for s in seg_a:
                    bounds.add(s.start)
                    bounds.add(s.end)
                for s in seg_b:
                    bounds.add(s.start)
                    bounds.add(s.end)

                if len(bounds) < 2:
                    continue
                b = sorted(bounds)
                cur_a = Cursor()
                cur_b = Cursor()

                for i in range(len(b) - 1):
                    start = b[i]
                    end = b[i + 1]
                    if end <= start:
                        continue
                    bp = end - start

                    sa = state_at(seg_a, cur_a, start)
                    sb = state_at(seg_b, cur_b, start)
                    evaluable = int(sa is not None and sb is not None)

                    a1 = sa.h1_anc if sa else ""
                    a2 = sa.h2_anc if sa else ""
                    b1 = sb.h1_anc if sb else ""
                    b2 = sb.h2_anc if sb else ""

                    m_ord = ""
                    m_unord = ""
                    if evaluable:
                        m_ord = int(a1 == b1 and a2 == b2)
                        m_unord = int(sorted((a1, a2)) == sorted((b1, b2)))

                    writer.writerow(
                        [
                            sample_id,
                            chrom,
                            tool_a,
                            tool_b,
                            start,
                            end,
                            bp,
                            a1,
                            a2,
                            b1,
                            b2,
                            evaluable,
                            m_ord,
                            m_unord,
                        ]
                    )
                    written_rows += 1

                    stat_key = (sample_id, chrom, tool_a, tool_b)
                    st = summary_stats[stat_key]
                    st["intervals"] += 1
                    st["bp_total"] += bp
                    if evaluable:
                        st["bp_evaluable"] += bp
                        st["bp_match_ordered"] += bp * int(m_ord)
                        st["bp_match_unordered"] += bp * int(m_unord)

    with open(summary_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "Sample_ID",
                "Chromosome",
                "Tool_A",
                "Tool_B",
                "Intervals",
                "BP_Total",
                "BP_Evaluable",
                "BP_Match_OrderedHaps",
                "BP_Match_UnorderedHaps",
                "Concordance_OrderedHaps",
                "Concordance_UnorderedHaps",
            ]
        )
        for key in sorted(summary_stats):
            st = summary_stats[key]
            bp_eval = st["bp_evaluable"]
            c_ord = "" if bp_eval == 0 else st["bp_match_ordered"] / bp_eval
            c_unord = "" if bp_eval == 0 else st["bp_match_unordered"] / bp_eval
            writer.writerow(
                [
                    key[0],
                    key[1],
                    key[2],
                    key[3],
                    st["intervals"],
                    st["bp_total"],
                    st["bp_evaluable"],
                    st["bp_match_ordered"],
                    st["bp_match_unordered"],
                    c_ord,
                    c_unord,
                ]
            )

    return len(summary_stats), written_rows


def main() -> int:
    args = parse_args()

    rfmix_segment_paths = list(args.rfmix_segments)
    if args.rfmix_dir:
        try:
            rfmix_segment_paths.extend(discover_segment_files(args.rfmix_dir))
        except Exception as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 1

    if not rfmix_segment_paths and not args.flare_segments:
        print(
            "ERROR: provide RFMix input via --rfmix-segments/--rfmix-dir and/or FLARE input via --flare-segments",
            file=sys.stderr,
        )
        return 1

    out_prefix = os.path.abspath(args.out_prefix)
    out_dir = os.path.dirname(out_prefix) or "."
    os.makedirs(out_dir, exist_ok=True)

    rfmix_trio_lookup: Dict[str, TrioMember] = {}
    if args.rfmix_trio_tsv:
        try:
            rfmix_trio_lookup = load_trio_lookup(args.rfmix_trio_tsv)
        except Exception as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 1

    try:
        segments, backfilled_rows = load_segments(
            rfmix_segment_paths + args.flare_segments,
            rfmix_trio_lookup=rfmix_trio_lookup,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if not segments:
        print("ERROR: no valid segment rows loaded", file=sys.stderr)
        return 1

    sort_segments(segments)

    unified_path = f"{out_prefix}.unified_segments.tsv"
    task2_intervals_path = f"{out_prefix}.task2_intervals.tsv"
    task2_summary_path = f"{out_prefix}.task2_summary.tsv"
    task3_intervals_path = f"{out_prefix}.task3_pairwise_intervals.tsv"
    task3_summary_path = f"{out_prefix}.task3_pairwise_summary.tsv"

    write_unified(unified_path, segments)
    task2_groups, task2_rows = build_task2(
        segments=segments,
        intervals_path=task2_intervals_path,
        summary_path=task2_summary_path,
    )
    task3_groups, task3_rows = build_task3(
        segments=segments,
        intervals_path=task3_intervals_path,
        summary_path=task3_summary_path,
    )

    print("Merge + task input build complete")
    print(f"  Loaded segment rows: {len(segments)}")
    if args.rfmix_dir:
        print(f"  RFMix files discovered from directory: {len(rfmix_segment_paths)}")
    if rfmix_trio_lookup:
        print(f"  RFMix rows with trio metadata backfilled: {backfilled_rows}")
    print(f"  Unified table: {unified_path}")
    print(f"  Task2 intervals: {task2_intervals_path}")
    print(f"  Task2 summary groups: {task2_groups} (rows={task2_rows})")
    print(f"  Task3 pairwise intervals: {task3_intervals_path}")
    print(f"  Task3 pairwise groups: {task3_groups} (rows={task3_rows})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
