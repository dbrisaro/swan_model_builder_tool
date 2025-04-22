import os
import sys
import yaml
from pathlib import Path
from datetime import datetime

# Add the parent directory to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.swan_functions import SwanModel, SwanConfig
from simulation.builder import SwanBuilder
from create_grid import create_swan_grid
from request_wave_data import main as request_wave_data
from request_wind_data import main as request_wind_data

def load_config(config_file='experiments_specs.txt'):
    """Load configuration from YAML file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def main():
    # Load configuration
    config = load_config()
    
    # Create output directory structure
    output_dir = f"/Users/daniela/Documents/swan/swan_experiments/{config['output']['directory']}"
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Create grids
    print("Creating grids...")
    regional_grid = create_swan_grid(
        lon_min=config['grids']['regional']['lon_min'],
        lon_max=config['grids']['regional']['lon_max'],
        lat_min=config['grids']['regional']['lat_min'],
        lat_max=config['grids']['regional']['lat_max'],
        n_cells_x=config['grids']['regional']['n_cells_x'],
        n_cells_y=config['grids']['regional']['n_cells_y'],
        grid_name=config['grids']['regional']['name'],
        rotation=config['grids']['regional']['rotation'],
        output_dir=os.path.join(output_dir, '01_geometry')
    )
    
    transition_grid = create_swan_grid(
        lon_min=config['grids']['transition']['lon_min'],
        lon_max=config['grids']['transition']['lon_max'],
        lat_min=config['grids']['transition']['lat_min'],
        lat_max=config['grids']['transition']['lat_max'],
        n_cells_x=config['grids']['transition']['n_cells_x'],
        n_cells_y=config['grids']['transition']['n_cells_y'],
        grid_name=config['grids']['transition']['name'],
        rotation=config['grids']['transition']['rotation'],
        output_dir=os.path.join(output_dir, '01_geometry')
    )
    
    # Step 2: Request wave data
    print("Requesting wave data...")
    request_wave_data()
    
    # Step 3: Request wind data
    print("Requesting wind data...")
    request_wind_data()
    
    # Step 4: Build and run SWAN simulation
    print("Building SWAN simulation...")
    builder = SwanBuilder(
        rootFolder=output_dir,
        templateSource='template',
        configSource=os.path.join(output_dir, 'CONFIG.ini'),
        gridSource=os.path.join(output_dir, '01_geometry'),
        bottomSource=None,  # Add if needed
        windSource=os.path.join(output_dir, '02_bc_dbase', 'WIND'),
        waveSource=os.path.join(output_dir, '02_bc_dbase', 'WAVE')
    )
    
    # Convert date strings to datetime objects
    time_start = datetime.strptime(config['time']['start_date'], '%Y-%m-%d')
    time_end = datetime.strptime(config['time']['end_date'], '%Y-%m-%d')
    
    # Build and run simulation
    builder.buildRun(
        timeStart=time_start,
        timeEnd=time_end,
        prefix=config['output']['prefix'],
        suffix=config['output']['suffix']
    )
    
    print("Simulation setup complete!")

if __name__ == "__main__":
    main() 