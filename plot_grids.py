"""Plot SWAN grids"""

import os
import yaml
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from functions.grid import create_grid_lines

def load_config(config_file='/Users/daniela/Documents/swan/swan_experiments/experiments_specs.txt'):
    """Load configuration from YAML file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def plot_grids():
    """Plot SWAN grids with map background"""
    # Load configuration
    config = load_config()
    
    # Create grid parameters dictionary for regional grid
    regional_params = {
        'lon_min': config['grids']['regional']['bounds']['lon_min'],
        'lon_max': config['grids']['regional']['bounds']['lon_max'],
        'lat_min': config['grids']['regional']['bounds']['lat_min'],
        'lat_max': config['grids']['regional']['bounds']['lat_max'],
        'dx': config['grids']['regional']['resolution']['dx'],
        'dy': config['grids']['regional']['resolution']['dy'],
        'name': config['grids']['regional']['name'],
        'rotation': 0.0
    }
    
    # Create grid parameters dictionary for transition grid
    transition_params = {
        'lon_min': config['grids']['transition']['bounds']['lon_min'],
        'lon_max': config['grids']['transition']['bounds']['lon_max'],
        'lat_min': config['grids']['transition']['bounds']['lat_min'],
        'lat_max': config['grids']['transition']['bounds']['lat_max'],
        'dx': config['grids']['transition']['resolution']['dx'],
        'dy': config['grids']['transition']['resolution']['dy'],
        'name': config['grids']['transition']['name'],
        'rotation': 0.0
    }
    
    # Create grid lines
    regional_grid = create_grid_lines(regional_params)
    transition_grid = create_grid_lines(transition_params)
    
    # Set plot limits with some margin
    lon_min = min(regional_params['lon_min'], transition_params['lon_min'])
    lon_max = max(regional_params['lon_max'], transition_params['lon_max'])
    lat_min = min(regional_params['lat_min'], transition_params['lat_min'])
    lat_max = max(regional_params['lat_max'], transition_params['lat_max'])
    margin = 0.5
    
    # Create figure with Cartopy projection
    proj = ccrs.PlateCarree()
    fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={'projection': proj})
    
    # Add map features
    ax.add_feature(cfeature.LAND, facecolor='lightgray', alpha=0.5)
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue', alpha=0.5)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    
    # Plot grids
    for line in regional_grid.geometry:
        ax.plot(*line.xy, color='blue', linewidth=0.5, alpha=0.7, transform=proj)
    for line in transition_grid.geometry:
        ax.plot(*line.xy, color='red', linewidth=0.5, alpha=0.7, transform=proj)
    
    # Add custom legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='blue', linewidth=1, label='Regional Grid'),
        Line2D([0], [0], color='red', linewidth=1, label='Transition Grid')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    # Set extent
    ax.set_extent([lon_min - margin, lon_max + margin, 
                   lat_min - margin, lat_max + margin], 
                  crs=proj)
    
    # Add title
    plt.title('SWAN Model Grids')
    
    # Save plot in QGIS directory
    output_dir = os.path.join('/Users/daniela/Documents/swan/swan_experiments',
                             config['output']['directory'],
                             'QGIS')
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, 'swan_grids.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    plot_grids() 