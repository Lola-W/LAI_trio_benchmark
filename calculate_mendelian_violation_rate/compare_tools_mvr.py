import pandas as pd
import argparse
import sys

def compare_tools_mvr(intervals_file, pop_file, ped_file, output_file="tool_comparison.csv"):
    df = pd.read_csv(intervals_file, sep='\t')
    pop_df = pd.read_csv(pop_file, sep='\t')
    pop_map = dict(zip(pop_df['sample_id'], pop_df['population']))
    
    ped_df = pd.read_csv(ped_file, sep=' ')
    ped_map = ped_df.set_index('sampleID')[['fatherID', 'motherID']].to_dict('index')

    def is_violation(row):
        child_h1 = row['Child_H1_Ancestry']
        child_h2 = row['Child_H2_Ancestry']
        child_id = row['Child_ID']
        
        if child_id in ped_map and ped_map[child_id]['fatherID'] != '0':
            father_id = ped_map[child_id]['fatherID']
            mother_id = ped_map[child_id]['motherID']
        else:
            father_id = row.get('Father_ID')
            mother_id = row.get('Mother_ID')
            
        father_pop = pop_map.get(father_id)
        mother_pop = pop_map.get(mother_id)
        
        if pd.isna(father_pop) or pd.isna(mother_pop):
            return pd.NA
            
        match_ordered = (child_h1 == father_pop) and (child_h2 == mother_pop)
        match_unordered = (child_h1 == mother_pop) and (child_h2 == father_pop)
        
        return not (match_ordered or match_unordered)

    df_evaluable = df.dropna(subset=['Child_H1_Ancestry', 'Child_H2_Ancestry']).copy()
    df_evaluable['Is_Violation'] = df_evaluable.apply(is_violation, axis=1)
    df_evaluable = df_evaluable.dropna(subset=['Is_Violation'])
    df_evaluable['Is_Violation'] = df_evaluable['Is_Violation'].astype(bool)
    df_evaluable.to_csv("/home/wzhang/pop_gen/data/evaluable_segments_with_violations.csv", index=False)

    child_tool_totals = df_evaluable.groupby(['Tool', 'Child_ID'])['BP_Length'].sum().rename('Total_Evaluable_Length')
    violations_df = df_evaluable[df_evaluable['Is_Violation']]
    child_tool_violations = violations_df.groupby(['Tool', 'Child_ID'])['BP_Length'].sum().reindex(child_tool_totals.index, fill_value=0).rename('Violation_Length')
    child_tool_violations.to_csv("/home/wzhang/pop_gen/data/child_tool_violations.csv", index=True)

    summary_df = pd.concat([child_tool_totals, child_tool_violations], axis=1)
    summary_df['Violation_Rate'] = summary_df['Violation_Length'] / summary_df['Total_Evaluable_Length']
    summary_df.reset_index(inplace=True)

    tool_summary = summary_df.groupby('Tool').agg(
        Num_Children=('Child_ID', 'nunique'),
        Mean_Violation_Rate=('Violation_Rate', 'mean'),
        Total_Evaluable_BP=('Total_Evaluable_Length', 'sum'),
        Total_Violation_BP=('Violation_Length', 'sum')
    ).reset_index()
    
    tool_summary['Overall_Violation_Rate'] = tool_summary['Total_Violation_BP'] / tool_summary['Total_Evaluable_BP']

    print(tool_summary.to_string(index=False))
    tool_summary.to_csv(output_file, index=False)
    print(f"\nTool comparison saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare Mendelian Violation Rates across tools using Pedigree and Population mapping.")
    
    parser.add_argument("intervals_file", help="Path to the input intervals TSV (e.g., merged_chr22.task2_intervals.tsv)")
    parser.add_argument("pop_file", help="Path to the population mapping TSV (e.g., sample_to_pop_superpop.tsv)")
    parser.add_argument("ped_file", help="Path to the pedigree file (e.g., 1kGP.3202_samples.pedigree_info.txt)")
    parser.add_argument("-o", "--output", default="mendelian_violation_rates_output.tsv", help="Path to save the output")
    
    args = parser.parse_args()
    
    compare_tools_mvr(args.intervals_file, args.pop_file, args.ped_file, args.output)