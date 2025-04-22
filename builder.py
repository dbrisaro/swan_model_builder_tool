import os
import sys
import yaml
from pathlib import Path
from swan_functions import SwanBuilder
from datetime import datetime

# Add the parent directory to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.swan_functions import SwanModel, SwanConfig
from grid.create_grid import create_swan_grid
from data.request_wave_data import main as request_wave_data
from data.request_wind_data import main as request_wind_data

def read_config(config_file):
    """Read experiment configuration from YAML file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def create_directory_structure(base_dir):
    """Create the required directory structure for the simulation."""
    dirs = [
        'DATA/BATHY',
        'DATA/WAVE',
        'DATA/WIND',
        'QGIS',
        'SWAN/01_geometry',
        'SWAN/02_bc_dbase',
        'SWAN/03_simulation',
        'SWAN/04_results'
    ]
    
    print("\nChecking directory structure...")
    all_dirs_exist = True
    for dir_path in dirs:
        full_path = base_dir / dir_path
        if full_path.exists():
            print(f"✓ Directory exists: {full_path}")
        else:
            all_dirs_exist = False
            print(f"✗ Directory missing: {full_path}")
    
    if all_dirs_exist:
        print("\nAll required directories already exist.")
        return True
    
    print("\nCreating missing directories...")
    for dir_path in dirs:
        full_path = base_dir / dir_path
        if not full_path.exists():
            full_path.mkdir(parents=True, exist_ok=True)
            print(f"Created directory: {full_path}")
    
    return True

def check_data_files(base_dir, config):
    """Check if wind and wave data files already exist."""
    # Generate expected filenames
    wave_file = base_dir / f"DATA/WAVE/wave_{config['time']['frequency']}_{config['time']['start']}_{config['time']['end']}_{config['grids']['regional']['bounds']['lat_min']}_{config['grids']['regional']['bounds']['lat_max']}_{config['grids']['regional']['bounds']['lon_min']}_{config['grids']['regional']['bounds']['lon_max']}.nc"
    wind_file = base_dir / f"DATA/WIND/wind_{config['time']['frequency']}_{config['time']['start']}_{config['time']['end']}_{config['grids']['regional']['bounds']['lat_min']}_{config['grids']['regional']['bounds']['lat_max']}_{config['grids']['regional']['bounds']['lon_min']}_{config['grids']['regional']['bounds']['lon_max']}.nc"
    
    files_exist = {
        'wave': wave_file.exists(),
        'wind': wind_file.exists()
    }
    
    if files_exist['wave']:
        print(f"✓ Wave data file already exists: {wave_file}")
    else:
        print(f"✗ Wave data file missing: {wave_file}")
    
    if files_exist['wind']:
        print(f"✓ Wind data file already exists: {wind_file}")
    else:
        print(f"✗ Wind data file missing: {wind_file}")
    
    return files_exist

def main():
    # Read configuration
    config_file = Path('/Users/daniela/Documents/swan/swan_experiments/experiments_specs.txt')
    config = read_config(config_file)
    
    # Set up base directory
    base_dir = Path(f"/Users/daniela/Documents/swan/swan_experiments/{config['output']['directory']}")
    base_dir.mkdir(parents=True, exist_ok=True)
    
    print("\n=== Starting SWAN Simulation Builder ===")
    print(f"Creating simulation in: {base_dir}")
    
    # Step 1: Create directory structure
    print("\nStep 1: Checking/Creating directory structure...")
    if not create_directory_structure(base_dir):
        print("Error with directory structure. Exiting.")
        return
    
    # Step 2: Check if data files already exist
    print("\nStep 2: Checking for existing data files...")
    files_exist = check_data_files(base_dir, config)
    
    # Step 3: Initialize SwanBuilder
    print("\nStep 3: Initializing SWAN Builder...")
    try:
        builder = SwanBuilder(
            rootFolder=str(base_dir / 'SWAN'),
            templateSource='/Users/daniela/Documents/swan/swan_model_builder_tool/template',
            configSource=str(base_dir / 'SWAN/CONFIG.ini'),
            gridSource=str(base_dir / f"QGIS/QGIS_{config['region']}_GRIDS.shp"),
            bottomSource=[str(base_dir / 'DATA/BATHY/GEBCO_01_Apr_2025_c2f0ebf57bc0/gebco_2024_n-8.0_s-10.0_w-80.0_e-78.0.tif')],
            windSource=str(base_dir / f"DATA/WIND/wind_{config['time']['frequency']}_{config['time']['start']}_{config['time']['end']}_{config['grids']['regional']['bounds']['lat_min']}_{config['grids']['regional']['bounds']['lat_max']}_{config['grids']['regional']['bounds']['lon_min']}_{config['grids']['regional']['bounds']['lon_max']}.nc"),
            waveSource=str(base_dir / f"DATA/WAVE/wave_{config['time']['frequency']}_{config['time']['start']}_{config['time']['end']}_{config['grids']['regional']['bounds']['lat_min']}_{config['grids']['regional']['bounds']['lat_max']}_{config['grids']['regional']['bounds']['lon_min']}_{config['grids']['regional']['bounds']['lon_max']}.nc")
        )
        print("SWAN Builder initialized successfully")
    except Exception as e:
        print(f"Error initializing SWAN Builder: {e}")
        return
    
    # Step 4: Generate grids
    print("\nStep 4: Generating grids...")
    try:
        builder.generate_grids()
        print("Grids generated successfully")
    except Exception as e:
        print(f"Error generating grids: {e}")
        return
    
    # Step 5: Generate CONFIG.ini
    print("\nStep 5: Generating CONFIG.ini...")
    try:
        builder.generate_config()
        print("CONFIG.ini generated successfully")
    except Exception as e:
        print(f"Error generating CONFIG.ini: {e}")
        return
    
    # Step 6: Request wave data (if needed)
    if not files_exist['wave']:
        print("\nStep 6: Requesting wave data...")
        try:
            builder.request_wave_data()
            print("Wave data requested successfully")
        except Exception as e:
            print(f"Error requesting wave data: {e}")
            return
    else:
        print("\nStep 6: Skipping wave data request (file already exists)")
    
    # Step 7: Request wind data (if needed)
    if not files_exist['wind']:
        print("\nStep 7: Requesting wind data...")
        try:
            builder.request_wind_data()
            print("Wind data requested successfully")
        except Exception as e:
            print(f"Error requesting wind data: {e}")
            return
    else:
        print("\nStep 7: Skipping wind data request (file already exists)")
    
    # Step 8: Build and run SWAN simulation
    print("\nStep 8: Building and running SWAN simulation...")
    try:
        # Convert dates from YYYY-MM-DD to DD/MM/YYYY format
        start_date = datetime.strptime(config['time']['start'], '%Y-%m-%d')
        end_date = datetime.strptime(config['time']['end'], '%Y-%m-%d')
        
        # Format dates as strings
        start_str = start_date.strftime('%d/%m/%Y %H:%M:%S')
        end_str = end_date.strftime('%d/%m/%Y %H:%M:%S')
        
        print(f"Running simulation from {start_str} to {end_str}")
        builder.buildRun(timeStart=start_str, timeEnd=end_str)
        print("SWAN simulation built and run successfully")
    except Exception as e:
        print(f"Error building and running SWAN simulation: {e}")
        return
    
    print("\n=== SWAN Simulation Builder completed successfully ===")
    print(f"Simulation directory: {base_dir}")
    print("\nSimulation results are available in the SWAN/04_results directory")

if __name__ == "__main__":
    main() 