import os
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse

def extract_best_fitness(file_path):
    """
    Extract the best fitness values from a single experiment file.
    Each file contains multiple runs of the algorithm.
    """
    best_fitness_values = []
    
    with open(file_path, 'r') as file:
        content = file.read()
        
        # Find all occurrences of "The best fitness is X"
        matches = re.findall(r'The best fitness is (\d+)', content)
        
        # Convert to integers
        best_fitness_values = [int(value) for value in matches]
        
    return best_fitness_values

def process_directory(directory_path):
    """
    Recursively scan a directory tree for exp.dat files and extract best fitness values.
    """
    results = {}
    
    # Walk through the directory tree
    for root, _, files in os.walk(directory_path):
        if 'exp.dat' in files:
            file_path = os.path.join(root, 'exp.dat')
            
            # Use the parent directory name as the experiment identifier
            experiment_name = os.path.basename(root)
            
            # Extract best fitness values
            fitness_values = extract_best_fitness(file_path)
            
            # Store results if we found any fitness values
            if fitness_values:
                results[experiment_name] = fitness_values
    
    return results

def create_boxplots(results, output_path='fitness_boxplots.png'):
    """
    Create box plots from the extracted fitness values.
    """
    if not results:
        print("No fitness data found.")
        return
    
    # Set up the plot
    plt.figure(figsize=(12, 8))
    
    # Create a boxplot
    sns.set_style("whitegrid")
    
    # Prepare data for boxplot
    data = []
    labels = []
    
    for exp_name, fitness_values in results.items():
        data.append(fitness_values)
        labels.append(exp_name)
    
    # Create the boxplot
    ax = sns.boxplot(data=data, palette="Set3")
    
    # Set labels and title
    plt.xlabel('Experiment')
    plt.ylabel('Best Fitness Value')
    plt.title('Distribution of Best Fitness Values Across Experiments')
    
    # Set x-axis tick labels
    ax.set_xticklabels(labels, rotation=45, ha='right')
    
    # Add a grid
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Tight layout to ensure everything fits
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(output_path)
    print(f"Boxplot saved to {output_path}")
    
    # Show the plot
    plt.show()

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process experiment files and create fitness boxplots.')
    parser.add_argument('directory', nargs='?', default='.', help='Directory containing experiment files (default: current directory)')
    parser.add_argument('--output', '-o', default='fitness_boxplots.png', help='Output file path for the boxplot image')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Process directory and create boxplots
    results = process_directory(args.directory)
    
    if not results:
        print(f"No exp.dat files found in {args.directory}")
        return
    
    print(f"Found {len(results)} experiment directories with exp.dat files.")
    for exp_name, values in results.items():
        print(f"  - {exp_name}: {len(values)} fitness values (max: {max(values)}, min: {min(values)})")
    
    create_boxplots(results, args.output)

if __name__ == "__main__":
    main()