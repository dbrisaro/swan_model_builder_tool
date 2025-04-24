#!/usr/bin/env python3
"""Run SWAN simulations for each grid"""

import os
import sys
import yaml
import platform
from pathlib import Path

def read_config(config_file):
    """Read experiment configuration from YAML file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def get_grid_names(config):
    """Get grid names in order from config file."""
    return [config['grids'][grid]['name'] for grid in config['grids']]

def find_swn_files(sim_dir, grid_names):
    """Find .swn files for each grid."""
    swn_files = {}
    for grid_name in grid_names:
        # Look for files containing the grid name
        matching_files = list(Path(sim_dir).glob(f"*{grid_name}*.swn"))
        if matching_files:
            swn_files[grid_name] = str(matching_files[0])
        else:
            print(f"Warning: No .swn file found for grid {grid_name}")
    return swn_files

def run_simulations(sim_dir, swn_files):
    """Run SWAN simulations for each grid."""
    # Change to simulation directory
    os.chdir(sim_dir)
    
    # Print files that will be executed
    print("\nWill execute the following .swn files:")
    for grid_name, swn_file in swn_files.items():
        print(f"- {swn_file}")
    
    # Ask for confirmation
    confirm = input("\nContinue with execution? (y/n): ")
    if confirm.lower() != 'y':
        print("Execution cancelled")
        return
    
    # Determine the appropriate swanrun command based on OS
    is_windows = platform.system() == 'Windows'
    
    # Run simulations in order
    for grid_name, swn_file in swn_files.items():
        print(f"\nRunning SWAN simulation for {swn_file}...")
        if is_windows:
            # Remove .swn extension for Windows
            file_name = Path(swn_file).stem
            os.system(f"swanrun {file_name}")
        else:
            os.system(f"swanrun -input {swn_file}")

def main():
    # Read configuration
    config_file = Path('/Users/daniela/Documents/swan/swan_experiments/experiments_specs.txt')
    config = read_config(config_file)
    
    # Get output directory
    base_dir = Path('/Users/daniela/Documents/swan/swan_experiments')
    output_dir = base_dir / config['output']['directory']
    sim_dir = output_dir / 'SWAN/03_simulation'
    
    # Get grid names
    grid_names = get_grid_names(config)
    print(f"Found grid names in config: {' '.join(grid_names)}")
    
    # Find .swn files
    swn_files = find_swn_files(sim_dir, grid_names)
    
    # Run simulations
    run_simulations(sim_dir, swn_files)

if __name__ == '__main__':
    main() 