import matplotlib.pyplot as plt
import numpy as np
import os
import re

def generate_convergence_plot(filename, output_dir="plots"):
    """
    Generate a convergence plot from a run file
    
    Parameters:
    filename (str): Path to the run file
    output_dir (str): Directory to save the plot
    """
    # Extract run number from filename
    match = re.search(r'run(\d+)\.dat', filename)
    if not match:
        print(f"Warning: Couldn't extract run number from {filename}")
        run_number = "unknown"
    else:
        run_number = match.group(1)
    
    # Read data
    with open(filename, 'r') as file:
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
    
    print(f"Plot saved: {output_path}")


def process_all_runs():
    """
    Find all run files in the current directory and generate plots for each
    """
    # Find all run files in the current directory
    run_files = [f for f in os.listdir('.') if re.match(r'run\d+\.dat', f)]
    
    if not run_files:
        print("No run files found. Make sure files are named as 'run[number].dat'")
        return
    
    # Sort files numerically
    run_files.sort(key=lambda x: int(re.search(r'run(\d+)\.dat', x).group(1)))
    
    print(f"Found {len(run_files)} run files. Generating plots...")
    
    # Process each file
    for filename in run_files:
        print(f"Processing {filename}...")
        generate_convergence_plot(filename)
    
    print("All plots generated successfully!")


if __name__ == "__main__":
    process_all_runs()