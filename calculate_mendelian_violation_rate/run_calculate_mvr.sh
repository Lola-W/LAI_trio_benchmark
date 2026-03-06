#!/bin/bash

PYTHON_SCRIPT="/home/wzhang/pop_gen/LAI_trio_benchmark/calculate_mendelian_violation_rate/calculate_mvr.py"
INTERVALS_FILE="/home/wzhang/pop_gen/data/merged_chr22.task2_intervals.tsv"
POP_FILE="/home/wzhang/pop_gen/data/sample_to_pop_superpop.tsv"
PED_FILE="/home/wzhang/pop_gen/data/1kGP.3202_samples.pedigree_info.txt"
OUTPUT_FILE="/home/wzhang/pop_gen/data/mendelian_violation_rates_output.tsv"

echo "Starting Mendelian Violation Rate calculation..."
echo "Using intervals: $INTERVALS_FILE"
echo "Using population map: $POP_FILE"
echo "Using pedigree: $PED_FILE"
echo "------------------------------------------------"

python3 "$PYTHON_SCRIPT" "$INTERVALS_FILE" "$POP_FILE" "$PED_FILE" -o "$OUTPUT_FILE"

echo "------------------------------------------------"
echo "Bash script finished! Check $OUTPUT_FILE for your results."