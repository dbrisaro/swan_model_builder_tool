"""Request wind data from ERA5"""

import os
import sys
import cdsapi
import yaml
import json
from pathlib import Path
from functions.data import generate_date_lists, generate_filename

def read_experiment_config(config_file):
    """Read experiment configuration from YAML file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def main():
    # Read configuration
    config_file = Path('/Users/daniela/Documents/swan/swan_experiments/experiments_specs.txt')
    config = read_experiment_config(config_file)
    
    # API key and client setup
    c = cdsapi.Client()
    
    # Get parameters from config
    time_config = config['time']
    start_date = time_config['start']
    end_date = time_config['end']
    frequency = time_config['frequency']
    
    # Get bounds from regional grid
    regional_grid = config['grids']['regional']['bounds']
    bounds = {
        'lon_min': regional_grid['lon_min'],
        'lon_max': regional_grid['lon_max'],
        'lat_min': regional_grid['lat_min'],
        'lat_max': regional_grid['lat_max']
    }
    
    # Generate dates
    dates = generate_date_lists(start_date, end_date, frequency)
    
    # Get unique years and months
    years = sorted(list(set(d[:4] for d in dates)))
    months = sorted(list(set(d[5:7] for d in dates)))
    
    # Base request parameters for wind
    request = {
        "variable": [
            "10m_u_component_of_wind",
            "10m_v_component_of_wind"
        ],
        "year": years,
        "month": months,
        "area": [
            bounds['lat_max'],
            bounds['lon_min'],
            bounds['lat_min'],
            bounds['lon_max']
        ],
        "format": "netcdf"
    }
    
    # Adjust request based on frequency
    if frequency == 'monthly':
        request["product_type"] = ["monthly_averaged_reanalysis"]
        request["time"] = [
            "00:00", "01:00", "02:00",
            "03:00", "04:00", "05:00",
            "06:00", "07:00", "08:00",
            "09:00", "10:00", "11:00",
            "12:00", "13:00", "14:00",
            "15:00", "16:00", "17:00",
            "18:00", "19:00", "20:00",
            "21:00", "22:00", "23:00"
        ]
        dataset = 'reanalysis-era5-single-levels-monthly-means'
    else:  # hourly
        request["product_type"] = ["reanalysis"]
        request["day"] = sorted(list(set(d[8:10] for d in dates)))
        request["time"] = [
            "00:00", "01:00", "02:00",
            "03:00", "04:00", "05:00",
            "06:00", "07:00", "08:00",
            "09:00", "10:00", "11:00",
            "12:00", "13:00", "14:00",
            "15:00", "16:00", "17:00",
            "18:00", "19:00", "20:00",
            "21:00", "22:00", "23:00"
        ]
        dataset = 'reanalysis-era5-single-levels'
    
    # Get output directory from config
    base_dir = Path('/Users/daniela/Documents/swan/swan_experiments')
    output_dir = base_dir / config['output']['directory'] / config['output']['data']['wind']
    output_dir.mkdir(parents=True, exist_ok=True)
    filename = output_dir / generate_filename('wind_data', frequency, start_date, end_date, bounds)
    
    # Print request details
    print("\nERA5 Request Details:")
    print("===================")
    print(f"Dataset: {dataset}")
    print("\nRequest parameters:")
    print(json.dumps(request, indent=2))
    print(f"\nOutput file: {filename}")
    
    # Ask for confirmation
    response = input("\nDo you want to proceed with this request? (yes/no): ")
    if response.lower() != 'yes':
        print("Request cancelled.")
        return
    
    # Request data
    c.retrieve(
        dataset,
        request,
        str(filename)
    )

if __name__ == '__main__':
    main() 