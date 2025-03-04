import matplotlib.pyplot as plt
import numpy as np
import os
import re
from collections import defaultdict

def generate_convergence_plot(filepath, output_dir):
    """
    Generate a convergence plot from a run file
    
    Parameters:
    filepath (str): Full path to the run file
    output_dir (str): Directory to save the plot
    """
    # Extract run number from filename
    filename = os.path.basename(filepath)
    match = re.search(r'run(\d+)\.dat', filename)
    if not match:
        print(f"Warning: Couldn't extract run number from {filename}")
        run_number = "unknown"
    else:
        run_number = match.group(1)
    
    # Read data
    with open(filepath, 'r') as file:
        lines = file.readlines()
    
    # Skip header and parse data
    data_lines = [line.strip() for line in lines if line.strip() and not line.startswith('Run')]
    
    # Parse the data
    iterations = []
    mean_values = []
    best_values = []
    
    for line in data_lines:
        parts = line.split()
        if len(parts) >= 7:  # Ensure we have enough columns
            iterations.append(int(parts[1]))  # RI column
            mean_values.append(float(parts[2]))  # Mean column
            best_values.append(float(parts[5]))  # Best column
    
    if not iterations:
        print(f"Warning: No valid data found in {filepath}")
        return
    
    # Create figure
    plt.figure(figsize=(12, 6))
    
    # Plot both lines
    plt.plot(iterations, mean_values, 'b-', label='Mean', linewidth=2)
    plt.plot(iterations, best_values, 'r-', label='Best', linewidth=2)
    
    # Add horizontal grid lines
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Set labels and title
    plt.xlabel('Reinforcement Iterations (RI)', fontsize=12)
    plt.ylabel('Value', fontsize=12)
    plt.title(f'Convergence Plot for Run {run_number}', fontsize=14)
    plt.legend(fontsize=12)
    
    # Set y-axis limits with padding
    ymin = min(min(mean_values), min(best_values)) - 1
    ymax = max(max(mean_values), max(best_values)) + 1
    plt.ylim(ymin, ymax)
    
    # Add annotations for final values
    plt.annotate(f'Final Mean: {mean_values[-1]:.2f}', 
                 xy=(iterations[-1], mean_values[-1]),
                 xytext=(iterations[-1]-20, mean_values[-1]+1),
                 arrowprops=dict(arrowstyle='->'),
                 fontsize=10)
    
    plt.annotate(f'Final Best: {best_values[-1]:.0f}', 
                 xy=(iterations[-1], best_values[-1]),
                 xytext=(iterations[-1]-20, best_values[-1]-2),
                 arrowprops=dict(arrowstyle='->'),
                 fontsize=10)
    
    plt.tight_layout()
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Save the figure
    output_path = os.path.join(output_dir, f'convergence_plot_run{run_number}.png')
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    return output_path


def process_directory_with_run_files(directory, run_files):
    """
    Process a directory containing run files and generate plots
    
    Parameters:
    directory (str): Path to the directory
    run_files (list): List of run files in the directory
    """
    print(f"\nProcessing directory: {directory}")
    print(f"Found {len(run_files)} run files")
    
    # Create an output directory for plots within the same directory as the run files
    output_dir = os.path.join(directory, "convergence_plots")
    
    # Process each file
    saved_plots = []
    for filename in run_files:
        filepath = os.path.join(directory, filename)
        print(f"  Processing {filename}...")
        output_path = generate_convergence_plot(filepath, output_dir)
        if output_path:
            saved_plots.append(output_path)
    
    print(f"Generated {len(saved_plots)} plots in {output_dir}")
    return len(saved_plots)


def recursive_search_and_plot():
    """
    Recursively walk through directories starting from current directory,
    looking for groups of run files to process
    """
    total_directories_processed = 0
    total_plots_generated = 0
    
    print("Starting recursive search for run files...")
    
    # Walk through all directories from current directory
    for root, dirs, files in os.walk('.'):
        # Look for run files in this directory
        run_files = [f for f in files if re.match(r'run\d+\.dat', f)]
        
        # If run files exist, process this directory
        if run_files:
            # Sort files numerically
            run_files.sort(key=lambda x: int(re.search(r'run(\d+)\.dat', x).group(1)))
            
            # Process this directory
            plots_generated = process_directory_with_run_files(root, run_files)
            
            total_directories_processed += 1
            total_plots_generated += plots_generated
    
    print("\nSummary:")
    print(f"Total directories processed: {total_directories_processed}")
    print(f"Total convergence plots generated: {total_plots_generated}")
    
    if total_directories_processed == 0:
        print("No directories with run files were found. Make sure files are named as 'run[number].dat'")


if __name__ == "__main__":
    recursive_search_and_plot()