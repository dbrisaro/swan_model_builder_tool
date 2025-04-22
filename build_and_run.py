import os
import sys
import yaml
from pathlib import Path
from swan_functions import SwanBuilder

def read_experiment_config(config_file):
    """Read experiment configuration from YAML file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def generate_file_name(prefix, start_date, end_date, lat_min, lat_max, lon_min, lon_max, extension):
    """Generate file name using grid boundaries and dates."""
    if prefix == 'gebco':
        return f"{prefix}_{start_date}_n{lat_max:.1f}_s{lat_min:.1f}_w{lon_min:.1f}_e{lon_max:.1f}{extension}"
    else:
        return f"{prefix}_data_monthly_{start_date[:4]}-{start_date[4:6]}-{start_date[6:]}_{end_date[:4]}-{end_date[4:6]}-{end_date[6:]}_{lon_min:.1f}_{lon_max:.1f}_{lat_min:.1f}_{lat_max:.1f}{extension}"

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

def main():
    # Read configuration
    config_file = Path('/Users/daniela/Documents/swan/swan_experiments/experiments_specs.txt')
    config = read_experiment_config(config_file)
    
    # Get grid boundaries from regional grid
    regional_grid = config['grids']['regional']
    lat_min = regional_grid['bounds']['lat_min']
    lat_max = regional_grid['bounds']['lat_max']
    lon_min = regional_grid['bounds']['lon_min']
    lon_max = regional_grid['bounds']['lon_max']
    
    # Get dates and format them for SWAN
    start_date = format_date_for_swan(config['time']['start'])
    end_date = format_date_for_swan(config['time']['end'])
    
    # Generate file names (using original dates for filenames)
    bathy_file = generate_file_name('gebco', '2024', '2024', lat_min, lat_max, lon_min, lon_max, '.tif')
    wind_file = generate_file_name('wind', 
                                 format_date_for_filename(config['time']['start']), 
                                 format_date_for_filename(config['time']['end']), 
                                 lat_min, lat_max, lon_min, lon_max, '.nc')
    wave_file = generate_file_name('wave', 
                                 format_date_for_filename(config['time']['start']), 
                                 format_date_for_filename(config['time']['end']), 
                                 lat_min, lat_max, lon_min, lon_max, '.nc')
    
    # Set up paths
    base_dir = Path('/Users/daniela/Documents/swan/swan_experiments')
    output_dir = base_dir / config['output']['directory']
    
    # Check if all required files exist
    bathy_path = output_dir / 'DATA/BATHY' / bathy_file
    wind_path = output_dir / 'DATA/WIND' / wind_file
    wave_path = output_dir / 'DATA/WAVE' / wave_file
    grid_path = output_dir / 'QGIS/swan_grids.shp'
    
    print(f"Checking required files:")
    print(f"Bathymetry: {bathy_path}")
    print(f"Wind: {wind_path}")
    print(f"Wave: {wave_path}")
    print(f"Grid: {grid_path}")
    
    if not all([
        check_file_exists(bathy_path, "Bathymetry"),
        check_file_exists(wind_path, "Wind"),
        check_file_exists(wave_path, "Wave"),
        check_file_exists(grid_path, "Grid")
    ]):
        print("Error: One or more required files are missing. Please check the paths above.")
        sys.exit(1)
    
    # Initialize builder
    builder = SwanBuilder(
        rootFolder=str(output_dir / 'SWAN'),
        templateSource=str(Path('/Users/daniela/Documents/swan/swan_model_builder_tool/template')),
        configSource=str(output_dir / 'SWAN/CONFIG.ini'),
        gridSource=str(grid_path),
        bottomSource=[str(bathy_path)],
        windSource=str(wind_path),
        waveSource=str(wave_path)
    )
    
    # Build and run
    builder.buildRun(
        timeStart=start_date,
        timeEnd=end_date
    )

if __name__ == '__main__':
    main() 