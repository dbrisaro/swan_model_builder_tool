"""Request bathymetry data from ETOPO"""

import os
import sys
import yaml
import json
import numpy as np
from pathlib import Path
from erddapy import ERDDAP
import rasterio
from functions.data import generate_date_lists, generate_filename

def read_experiment_config(config_file):
    """Read experiment configuration from YAML file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)
    
def main(config):

    # Get parameters from config
    time_config = config['time']
    start_date = time_config['start']
    end_date = time_config['end']
    frequency = time_config['frequency']

    # Get bounds from regional grid
    regional_grid = config['grids']['regional']['bounds']
    bounds = {
        "lon_min": regional_grid['lon_min'],
        "lon_max": regional_grid['lon_max'],
        "lat_min": regional_grid['lat_min'],
        "lat_max": regional_grid['lat_max']
    }
    bounds_request = {
        "longitude>=": float(np.mod(regional_grid['lon_min'],360)),
        "longitude<=": float(np.mod(regional_grid['lon_max'],360)),
        "latitude>=": regional_grid['lat_min'],
        "latitude<=": regional_grid['lat_max']
    }

    # Base request parameters for bathymetry
    e = ERDDAP(
        server="https://oceanwatch.pifsc.noaa.gov/erddap",
        protocol="griddap",
    )
    dataset = "ETOPO_2022_v1_15s"
    e.dataset_id = dataset

    # Update bounds
    e.constraints.update(bounds_request)  

    base_path = Path(config['base']['path'])
    output_dir = base_path / config['output']['directory'] / config['output']['data']['bathy']
    output_dir.mkdir(parents=True, exist_ok=True)
    filename = output_dir / generate_filename('bathy_data', frequency, start_date, end_date, bounds)

    # Print request details
    print("\nETOPO Request Details:")
    print("===================")
    print(f"Dataset: {dataset}")
    print("\nRequest parameters:")
    print(json.dumps(e.constraints, indent=2))
    print(f"\nOutput file: {filename}")
    
    # Ask for confirmation
    response = input("\nDo you want to proceed with this request? (yes/no): ")
    if response.lower() not in ('Y', 'y', 'yes', 'YES', 'S', 's'):
        print("Request cancelled.")
        return
    
    # Request data
    ds = e.to_xarray()
    ds['z'].squeeze().rename({
        'latitude': 'y', 
        'longitude': 'x'
    }).rio.write_crs("EPSG:4326").rio.to_raster(filename)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run SWAN config generator.')
    parser.add_argument('--config', type=str, required=True, help='Path to the experiment specifications file')
    args = parser.parse_args()

    main(args.config)