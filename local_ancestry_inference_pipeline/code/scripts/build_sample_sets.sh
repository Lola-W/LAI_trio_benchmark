#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PED="${ROOT_DIR}/1kGP.3202_samples.pedigree_info.txt"
OUT_DIR="${ROOT_DIR}/code/meta"

if [[ ! -f "${PED}" ]]; then
  echo "Missing pedigree file: ${PED}" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"

# 1) Complete trios: child has both parents in pedigree.
awk 'BEGIN{OFS="\t"}
NR==1{next}
$2!="0" && $3!="0" {
  n++;
  trio_id=sprintf("TRIO%04d", n);
  print trio_id,$1,$2,$3
}' "${PED}" > "${OUT_DIR}/trios.complete.tsv"

# 2) Children from complete trios.
awk 'NR==1{next} $2!="0" && $3!="0" {print $1}' "${PED}" > "${OUT_DIR}/samples.children.tsv"

# 3) Unique parent IDs referenced by complete trios (first-seen order).
awk 'NR==1{next}
$2!="0" && $3!="0" {
  if(!seen[$2]++){print $2}
  if(!seen[$3]++){print $3}
}' "${PED}" > "${OUT_DIR}/samples.parents.tsv"

# 4) Unique trio members = children + referenced parents.
awk 'NR==1{next}
$2!="0" && $3!="0" {
  if(!seen[$1]++){print $1}
  if(!seen[$2]++){print $2}
  if(!seen[$3]++){print $3}
}' "${PED}" > "${OUT_DIR}/samples.trio_members.tsv"

# 5) All pedigree sample IDs.
awk 'NR==1{next}{print $1}' "${PED}" > "${OUT_DIR}/samples.all_3202.tsv"

# 6) Strict masking:
#    query   = all trio members
#    ref     = all samples not in query
cp "${OUT_DIR}/samples.trio_members.tsv" "${OUT_DIR}/query.strict.tsv"
awk 'NR==FNR {in_query[$1]=1; next} !in_query[$1] {print $1}' \
  "${OUT_DIR}/query.strict.tsv" \
  "${OUT_DIR}/samples.all_3202.tsv" > "${OUT_DIR}/reference.strict.tsv"

# 7) Child-only masking:
#    query   = complete-trio children
#    ref     = all samples not in query (parents retained in reference)
cp "${OUT_DIR}/samples.children.tsv" "${OUT_DIR}/query.child_only.tsv"
awk 'NR==FNR {in_query[$1]=1; next} !in_query[$1] {print $1}' \
  "${OUT_DIR}/query.child_only.tsv" \
  "${OUT_DIR}/samples.all_3202.tsv" > "${OUT_DIR}/reference.child_only.tsv"

# 8) QC summaries.
{
  printf "metric\tvalue\n"
  awk 'END{printf "count.complete_trios\t%s\n", NR}' "${OUT_DIR}/trios.complete.tsv"
  awk 'END{printf "count.children\t%s\n", NR}' "${OUT_DIR}/samples.children.tsv"
  awk 'END{printf "count.unique_parents\t%s\n", NR}' "${OUT_DIR}/samples.parents.tsv"
  awk 'END{printf "count.unique_trio_members\t%s\n", NR}' "${OUT_DIR}/samples.trio_members.tsv"
  awk 'END{printf "count.reference.strict\t%s\n", NR}' "${OUT_DIR}/reference.strict.tsv"
  awk 'END{printf "count.reference.child_only\t%s\n", NR}' "${OUT_DIR}/reference.child_only.tsv"
} > "${OUT_DIR}/summary.counts.tsv"

# 9) Sanity checks: query/reference sets must be disjoint.
awk 'NR==FNR{q[$1]=1;next} ($1 in q){bad=1} END{exit bad}' \
  "${OUT_DIR}/query.strict.tsv" "${OUT_DIR}/reference.strict.tsv"
awk 'NR==FNR{q[$1]=1;next} ($1 in q){bad=1} END{exit bad}' \
  "${OUT_DIR}/query.child_only.tsv" "${OUT_DIR}/reference.child_only.tsv"

echo "Wrote preprocessing outputs to: ${OUT_DIR}"
awk 'NR==1 || NR>1 {print}' "${OUT_DIR}/summary.counts.tsv"

