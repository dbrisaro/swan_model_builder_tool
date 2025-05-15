import os
import sys
import yaml
from pathlib import Path
from swan_functions import SwanBuilder
from functions.mdatetime import *
from functions.data import generate_filename

def check_file_exists(file_path, file_type):
    """Check if a file exists and return its path if it does."""
    if not os.path.exists(file_path):
        print(f"Error: {file_type} file not found at {file_path}")
        return False
    return True

def format_date_for_swan(date_obj):
    """Convert date to DD/MM/YYYY HH:MM:SS format."""
    from datetime import datetime
    if isinstance(date_obj, str):
        date_obj = datetime.strptime(date_obj, '%Y-%m-%d')
    return date_obj.strftime('%d/%m/%Y 00:00:00')

def format_date_for_filename(date_obj):
    """Convert date to YYYYMMDD format for filenames."""
    return date_obj.strftime('%Y%m%d')

def main(config):

    start_date = config['analysis']['dates']['start']
    end_date = config['analysis']['dates']['end']
    regional_grid = config['grids']['regional'] 
    lon_min = regional_grid['bounds']['lon_min']
    lon_max = regional_grid['bounds']['lon_max']
    lat_min = regional_grid['bounds']['lat_min']
    lat_max = regional_grid['bounds']['lat_max']

    base_path = Path(config['base']['path'])
    output_dir = Path(config['output']['directory'])
    wind_dir = config['output']['data']['wind']
    wave_dir = config['output']['data']['wave']
    bathy_dir = config['output']['data']['bathy']
    # Get frequency from config
    frequency = config['time']['frequency']
    filename_wave = generate_filename(
        'wave_data', frequency, start_date, end_date,
        {'lon_min': lon_min, 'lon_max': lon_max, 'lat_min': lat_min, 'lat_max': lat_max}
    )
    filename_wind = generate_filename(
        'wind_data', frequency, start_date, end_date,
        {'lon_min': lon_min, 'lon_max': lon_max, 'lat_min': lat_min, 'lat_max': lat_max}
    )

    # Generar nombres de archivo usando generate_filename
    # filename_bathy = f"gebco_{start_date_file}_{lat_min}_{lat_max}_{lon_min}_{lon_max}.tif"
    
    
    output_dir_wave = base_path / output_dir / wave_dir
    output_path_wave = output_dir_wave / filename_wave

    output_dir_wind = base_path / output_dir / wind_dir
    output_path_wind = output_dir_wind / filename_wind

    # output_dir_bathy = base_path / output_dir / bathy_dir
    # output_path_bathy = output_dir_bathy / filename_bathy

    # Check if all required files exist
    # bathy_path = output_dir / 'DATA/BATHY' / bathy_file

    wind_path = output_path_wind
    wave_path = output_path_wave
    grid_path = base_path / output_dir / 'QGIS/swan_grids.shp'

    print(f"Checking required files:")
    # print(f"Bathymetry: {bathy_path}")
    print(f"Wind: {wind_path}")
    print(f"Wave: {wave_path}")
    print(f"Grid: {grid_path}")

    if not all([
        # check_file_exists(bathy_path, "Bathymetry"),
        check_file_exists(wind_path, "Wind"),
        check_file_exists(wave_path, "Wave"),
        check_file_exists(grid_path, "Grid")
    ]):
        print("Error: One or more required files are missing. Please check the paths above.")
        sys.exit(1)
    
    # print(output_dir)
    # # Initialize builder
    # builder = SwanBuilder(
    #     rootFolder=str(output_dir / 'SWAN'),
    #     templateSource=str(Path(__file__).parent / 'template'),
    #     configSource=str(output_dir / 'SWAN/CONFIG.ini'),
    #     gridSource=str(grid_path),
    #     bottomSource=[str(bathy_path)],
    #     windSource=str(wind_path),
    #     waveSource=str(wave_path)
    # )
    
    # # Build and run
    # models = builder.buildRun(
    #     timeStart=datenum(start_date),
    #     timeEnd=datenum(end_date)
    # )
    # return models
    return None

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run SWAN config generator.')
    parser.add_argument('--config', type=str, required=True, help='Path to the experiment specifications file')
    args = parser.parse_args()

    main(args.config)