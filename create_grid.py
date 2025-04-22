import sys
import os
import re
import yaml
from pathlib import Path

# Add the functions directory to the path
sys.path.append('functions/')

# Import necessary libraries
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon

# Import grid function
from swan_functions import create_rectangular_grid

def parse_experiments_specs(file_path):
    """
    Parse the experiments_specs.txt file to extract grid parameters.
    
    Parameters:
    -----------
    file_path : str or Path
        Path to the experiments_specs.txt file
    
    Returns:
    --------
    dict
        Dictionary containing grid parameters for each grid type
    """
    with open(file_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Initialize parameters dictionary
    params = {}
    
    # Extract regional grid parameters
    regional = config['grids']['regional']
    params['regional'] = {
        'lon_min': float(regional['bounds']['lon_min']),
        'lon_max': float(regional['bounds']['lon_max']),
        'lat_min': float(regional['bounds']['lat_min']),
        'lat_max': float(regional['bounds']['lat_max']),
        'dx': float(regional['resolution']['dx']),
        'dy': float(regional['resolution']['dy']),
        'name': regional['name'],
        'rotation': 0
    }
    
    # Extract transition grid parameters
    transition = config['grids']['transition']
    params['transition'] = {
        'lon_min': float(transition['bounds']['lon_min']),
        'lon_max': float(transition['bounds']['lon_max']),
        'lat_min': float(transition['bounds']['lat_min']),
        'lat_max': float(transition['bounds']['lat_max']),
        'dx': float(transition['resolution']['dx']),
        'dy': float(transition['resolution']['dy']),
        'name': transition['name'],
        'rotation': 0
    }
    
    return params

def create_swan_grid(grid_params, output_dir):
    """
    Creates a SWAN grid based on the provided parameters.
    
    Parameters:
    -----------
    grid_params : dict
        Dictionary containing grid parameters:
        - lon_min : float
            Minimum longitude of the grid
        - lon_max : float
            Maximum longitude of the grid
        - lat_min : float
            Minimum latitude of the grid
        - lat_max : float
            Maximum latitude of the grid
        - dx : float
            Grid resolution in x-direction
        - dy : float
            Grid resolution in y-direction
        - name : str
            Name of the grid
        - rotation : float, optional
            Rotation angle in degrees (default: 0)
    output_dir : str or Path
        Directory where the grid will be saved
    
    Returns:
    --------
    gdf : GeoDataFrame
        The created grid as a GeoDataFrame
    """
    # Calculate number of cells based on resolution
    nx = int((grid_params['lon_max'] - grid_params['lon_min']) / grid_params['dx'])
    ny = int((grid_params['lat_max'] - grid_params['lat_min']) / grid_params['dy'])
    
    # Create the grid
    grid = create_rectangular_grid(
        lon_min=grid_params['lon_min'],
        lon_max=grid_params['lon_max'],
        lat_min=grid_params['lat_min'],
        lat_max=grid_params['lat_max'],
        x_len=nx,
        y_len=ny,
        name=grid_params['name'],
        rotation=grid_params.get('rotation', 0)
    )
    
    # Create a polygon for the grid boundary
    polygon = Polygon([
        (grid_params['lon_min'], grid_params['lat_min']),
        (grid_params['lon_max'], grid_params['lat_min']),
        (grid_params['lon_max'], grid_params['lat_max']),
        (grid_params['lon_min'], grid_params['lat_max']),
        (grid_params['lon_min'], grid_params['lat_min'])
    ])
    
    # Convert SwanGrid to GeoDataFrame
    gdf = gpd.GeoDataFrame(
        {
            'name': [grid_params['name']],
            'dx': [grid_params['dx']],
            'dy': [grid_params['dy']],
            'nx': [nx],
            'ny': [ny],
            'rotation': [grid_params.get('rotation', 0)],
            'geometry': [polygon]
        },
        crs='EPSG:4326'  # WGS84
    )
    
    # Save the grid as shapefile
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    output_file = output_dir / f"{grid_params['name']}.shp"
    gdf.to_file(output_file, driver='ESRI Shapefile')
    
    print(f"\nGrid '{grid_params['name']}' created successfully!")
    print(f"Number of cells: {nx} x {ny}")
    print(f"Resolution: {grid_params['dx']}° x {grid_params['dy']}°")
    print(f"Shapefile saved to: {output_file}")
    
    return gdf

if __name__ == "__main__":
    # Define paths
    specs_file = Path('/Users/daniela/Documents/swan/swan_experiments/experiments_specs.txt')
    output_dir = Path('/Users/daniela/Documents/swan/swan_experiments/run_climatology_southern_peru/QGIS')
    
    # Parse parameters from experiments_specs.txt
    grid_params = parse_experiments_specs(specs_file)
    
    # Create all grids
    for grid_type, params in grid_params.items():
        create_swan_grid(params, output_dir) 