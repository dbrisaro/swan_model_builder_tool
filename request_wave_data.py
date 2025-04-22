"""Request wave data from ERA5"""

import os
import sys
import cdsapi
from datetime import datetime, timedelta
import yaml
from pathlib import Path
from functions.data import generate_date_lists, generate_filename

# Add the parent directory to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.swan_functions import SwanModel, SwanConfig

def read_experiment_config(config_file):
    """Read experiment configuration from YAML file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def main():
    # API key and client setup
    c = cdsapi.Client()
    
    # Parameters
    variable = 'wave_height'
    frequency = 'daily'
    start_date = '2020-01-01'
    end_date = '2020-12-31'
    bounds = {
        'lon_min': -85.0,
        'lon_max': -70.0,
        'lat_min': -20.0,
        'lat_max': -10.0
    }
    
    # Generate dates and filename
    dates = generate_date_lists(start_date, end_date, frequency)
    filename = generate_filename(variable, frequency, start_date, end_date, bounds)
    
    # Request data
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'variable': 'significant_height_of_combined_wind_waves_and_swell',
            'year': [d[:4] for d in dates],
            'month': [d[5:7] for d in dates],
            'day': [d[8:10] for d in dates],
            'time': '00:00',
            'area': [
                bounds['lat_max'],
                bounds['lon_min'],
                bounds['lat_min'],
                bounds['lon_max']
            ],
            'format': 'netcdf'
        },
        filename
    )

if __name__ == '__main__':
    main() 