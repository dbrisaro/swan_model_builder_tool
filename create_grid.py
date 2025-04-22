"""Create SWAN grid files"""

import os
import sys
import yaml
from pathlib import Path
import geopandas as gpd
import pandas as pd
from swan_functions import create_rectangular_grid

def parse_experiments_specs(specs_file):
    """Parse experiments specifications file."""
    with open(specs_file, 'r') as f:
        return yaml.safe_load(f)

def main():
    # Parse experiments specifications
    specs_file = Path('/Users/daniela/Documents/swan/swan_experiments/experiments_specs.txt')
    specs = parse_experiments_specs(specs_file)
    
    # Get grid parameters
    regional_grid = specs['grids']['regional']
    transition_grid = specs['grids']['transition']
    
    # Create output directory
    base_dir = Path('/Users/daniela/Documents/swan/swan_experiments')
    output_dir = base_dir / specs['output']['directory'] / 'QGIS'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create regional grid
    regional = create_rectangular_grid(
        regional_grid['bounds']['lon_min'],
        regional_grid['bounds']['lon_max'],
        regional_grid['bounds']['lat_min'],
        regional_grid['bounds']['lat_max'],
        regional_grid['resolution']['dx'],
        regional_grid['resolution']['dy'],
        regional_grid['name']
    )
    regional['grid_type'] = 'regional'  # Add type field
    
    # Create transition grid
    transition = create_rectangular_grid(
        transition_grid['bounds']['lon_min'],
        transition_grid['bounds']['lon_max'],
        transition_grid['bounds']['lat_min'],
        transition_grid['bounds']['lat_max'],
        transition_grid['resolution']['dx'],
        transition_grid['resolution']['dy'],
        transition_grid['name']
    )
    transition['grid_type'] = 'transition'  # Add type field
    
    # Combine grids into a single GeoDataFrame
    combined_grids = gpd.GeoDataFrame(pd.concat([regional, transition], ignore_index=True))
    
    # Save combined grids
    combined_grids.to_file(output_dir / 'swan_grids.shp')
    print(f"Grids saved to: {output_dir / 'swan_grids.shp'}")

if __name__ == '__main__':
    main() 