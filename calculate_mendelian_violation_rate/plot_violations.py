import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

def plot_violation_rates(df, bin_size=1_000_000, output_file=None):
    print("Cleaning data and extracting valid coordinates...")
    
    # Ensure Is_Violation is boolean
    if 'Is_Violation' in df.columns:
        # Handle cases where pandas might read it as a string 'True'/'False' instead of boolean
        df['Is_Violation'] = df['Is_Violation'].replace({'True': True, 'False': False, '1': True, '0': False}).astype(bool)
    else:
        print("Error: 'Is_Violation' column not found in data.")
        sys.exit(1)
    
    # Force Start and End columns to be numeric (turns bad data into NaNs)
    df['Start'] = pd.to_numeric(df['Start'], errors='coerce')
    df['End'] = pd.to_numeric(df['End'], errors='coerce')
    
    # Drop any rows where Start or End became NaN due to formatting issues
    df = df.dropna(subset=['Start', 'End'])
    
    # Ensure they are integers
    df['Start'] = df['Start'].astype(int)
    df['End'] = df['End'].astype(int)
        
    # Determine the start and end of the chromosome based on data
    min_pos = df['Start'].min()
    max_pos = df['End'].max()
    
    # Create bins
    bins = np.arange(min_pos, max_pos + bin_size, bin_size)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    
    tools = df['Tool'].dropna().unique()
    if len(tools) == 0:
        print("No tools found in the 'Tool' column.")
        sys.exit(1)
        
    print(f"Processing data into {bin_size / 1e6} Mb bins for tools: {', '.join(tools)}...")
    
    # Dictionary to store violation rates for each tool
    tool_rates = {tool: np.zeros(len(bins)-1) for tool in tools}
    
    for tool in tools:
        tool_df = df[df['Tool'] == tool]
        
        # Arrays to accumulate total bp and violated bp per bin
        total_bp_per_bin = np.zeros(len(bins)-1)
        violation_bp_per_bin = np.zeros(len(bins)-1)
        
        for _, row in tool_df.iterrows():
            start, end = row['Start'], row['End']
            is_violation = row['Is_Violation']
            
            # Find which bins this segment overlaps (cap indices to valid range)
            start_idx = max(0, np.searchsorted(bins, start, side='right') - 1)
            end_idx = min(len(bins)-1, np.searchsorted(bins, end, side='left'))
            
            for i in range(start_idx, end_idx):
                bin_start = bins[i]
                bin_end = bins[i+1]
                
                # Calculate overlap between segment and bin
                overlap_len = max(0, min(end, bin_end) - max(start, bin_start))
                
                total_bp_per_bin[i] += overlap_len
                if is_violation:
                    violation_bp_per_bin[i] += overlap_len
                    
        # Calculate rate (avoid division by zero)
        with np.errstate(divide='ignore', invalid='ignore'):
            rate = np.where(total_bp_per_bin > 0, violation_bp_per_bin / total_bp_per_bin, 0)
            
        tool_rates[tool] = rate

    # --- Plotting ---
    plt.figure(figsize=(14, 6))
    
    for tool in tools:
        plt.plot(bin_centers / 1e6, tool_rates[tool], marker='o', linestyle='-', alpha=0.8, label=tool)
        
    plt.title(f'Mendelian Violation Rate along Chromosome (Bin size: {bin_size // 1000} kb)', fontsize=14)
    plt.xlabel('Chromosome Position (Mb)', fontsize=12)
    plt.ylabel('Violation Rate (Violated BP / Total BP)', fontsize=12)
    plt.legend(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    
    # Save or show
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot successfully saved to {output_file}")
    else:
        plt.show()

if __name__ == "__main__":
    # Setup Argument Parser
    parser = argparse.ArgumentParser(description="Plot Mendelian violation rates from a CSV file.")
    
    # Define arguments
    parser.add_argument("-i", "--input", required=True, 
                        help="Path to the input results file (CSV format expected).")
    parser.add_argument("-b", "--bin_size", type=int, default=1_000_000, 
                        help="Bin size in base pairs (default: 1000000 for 1Mb).")
    parser.add_argument("-o", "--output", default=None, 
                        help="Path to save the output plot (e.g., 'violation_plot.png').")

    # Parse arguments
    args = parser.parse_args()

    # Load data and run
    print(f"Loading data from {args.input}...")
    try:
        # Note the default separator is now a comma (,) for CSV
        data = pd.read_csv(args.input, sep=',')
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)
        
    plot_violation_rates(data, bin_size=args.bin_size, output_file=args.output)