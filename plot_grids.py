"""Plot SWAN grids"""

import os
import sys
import yaml
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pathlib import Path
from functions.grid import create_grid_lines

def read_experiment_config(config_file):
    """Read experiment configuration from YAML file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def main():
    # Read configuration
    config_file = Path('/Users/daniela/Documents/swan/swan_experiments/experiments_specs.txt')
    config = read_experiment_config(config_file)
    
    # Get grid parameters
    regional_grid = config['grids']['regional']
    transition_grid = config['grids']['transition']
    
    # Create figure with cartopy projection
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    
    # Add coastlines and land
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    
    # Create and plot regional grid lines
    regional_lines = create_grid_lines({
        'lon_min': regional_grid['bounds']['lon_min'],
        'lon_max': regional_grid['bounds']['lon_max'],
        'lat_min': regional_grid['bounds']['lat_min'],
        'lat_max': regional_grid['bounds']['lat_max'],
        'dx': regional_grid['resolution']['dx'],
        'dy': regional_grid['resolution']['dy'],
        'name': regional_grid['name'],
        'rotation': 0.0
    })
    regional_lines.plot(ax=ax, color='red', linewidth=1.0, alpha=0.7)
    
    # Create and plot transition grid lines
    transition_lines = create_grid_lines({
        'lon_min': transition_grid['bounds']['lon_min'],
        'lon_max': transition_grid['bounds']['lon_max'],
        'lat_min': transition_grid['bounds']['lat_min'],
        'lat_max': transition_grid['bounds']['lat_max'],
        'dx': transition_grid['resolution']['dx'],
        'dy': transition_grid['resolution']['dy'],
        'name': transition_grid['name'],
        'rotation': 0.0
    })
    transition_lines.plot(ax=ax, color='blue', linewidth=1.0, alpha=0.7)
    
    # Set plot limits
    ax.set_xlim(regional_grid['bounds']['lon_min'] - 0.1, regional_grid['bounds']['lon_max'] + 0.1)
    ax.set_ylim(regional_grid['bounds']['lat_min'] - 0.1, regional_grid['bounds']['lat_max'] + 0.1)
    
    # Add grid
    ax.gridlines(draw_labels=True)
    
    # Add title and legend
    plt.title('SWAN Grids')
    ax.plot([], [], color='red', label='Regional Grid')
    ax.plot([], [], color='blue', label='Transition Grid')
    ax.legend()
    
    # Save plot
    output_dir = Path('/Users/daniela/Documents/swan/swan_experiments') / config['output']['directory'] / 'QGIS'
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / 'swan_grids.png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    main() 