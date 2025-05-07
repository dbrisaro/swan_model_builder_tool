"""Plot SWAN grids"""

import os
import sys
import yaml
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pathlib import Path
import geopandas as gpd
import numpy as np
import configparser
import pandas as pd
from math import cos, sin, radians

def parse_specs(config):
    # Get rotation
    rotation = config['rotation']
    
    # Get regional grid specs
    regional = config['REGIONAL GRID']
    regional_specs = {
        'lon_min': float(regional['lon_min']),
        'lon_max': float(regional['lon_max']),
        'lat_min': float(regional['lat_min']),
        'lat_max': float(regional['lat_max']),
        'x_len': float(regional['x_len']),
        'y_len': float(regional['y_len']),
        'rotation': float(rotation)
    }
    
    # Get transition grid specs
    transition = config['TRANSITION GRID']
    transition_specs = {
        'lon_min': float(transition['lon_min']),
        'lon_max': float(transition['lon_max']),
        'lat_min': float(transition['lat_min']),
        'lat_max': float(transition['lat_max']),
        'x_len': float(transition['x_len']),
        'y_len': float(transition['y_len']),
        'rotation': float(rotation)
    }
    
    return regional_specs, transition_specs

def get_rotated_grid_points(lon_min, lon_max, lat_min, lat_max, x_len, y_len, rotation):
    center_lon = (lon_min + lon_max) / 2
    center_lat = (lat_min + lat_max) / 2
    lons = np.arange(lon_min, lon_max + x_len, x_len)
    lats = np.arange(lat_min, lat_max + y_len, y_len)
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    points = np.column_stack([lon_grid.flatten(), lat_grid.flatten()])
    angle = radians(rotation)
    rot_points = []
    for lon, lat in points:
        x = lon - center_lon
        y = lat - center_lat
        x_rot = x * cos(angle) - y * sin(angle)
        y_rot = x * sin(angle) + y * cos(angle)
        rot_points.append((x_rot + center_lon, y_rot + center_lat))
    return np.array(rot_points)

def main(config):
    # Create figure with cartopy projection
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    
    # Add coastlines and land
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    
    # Get rotation
    rotation = config['rotation']

    # Regional grid
    reg = config['grids']['regional']
    reg_points = get_rotated_grid_points(
        reg['bounds']['lon_min'], reg['bounds']['lon_max'],
        reg['bounds']['lat_min'], reg['bounds']['lat_max'],
        reg['resolution']['dx'], reg['resolution']['dy'],
        rotation=rotation
    )
    ax.scatter(reg_points[:,0], reg_points[:,1], color='red', s=2, alpha=0.7, 
               label='Regional Grid', transform=ccrs.PlateCarree())

    # Transition grid
    trans = config['grids']['transition']
    trans_points = get_rotated_grid_points(
        trans['bounds']['lon_min'], trans['bounds']['lon_max'],
        trans['bounds']['lat_min'], trans['bounds']['lat_max'],
        trans['resolution']['dx'], trans['resolution']['dy'],
        rotation=rotation
    )
    ax.scatter(trans_points[:,0], trans_points[:,1], color='blue', s=2, alpha=0.7, 
               label='Transition Grid', transform=ccrs.PlateCarree())

    # Set plot limits with some padding
    ax.set_extent([
        reg['bounds']['lon_min'] - 1, 
        reg['bounds']['lon_max'] + 1,
        reg['bounds']['lat_min'] - 1, 
        reg['bounds']['lat_max'] + 1
    ], crs=ccrs.PlateCarree())

    # Add gridlines
    ax.gridlines(draw_labels=True)
    
    ax.set_title('SWAN Grids')
    ax.legend()

    base_path = Path(config['base']['path'])
    output_dir = base_path / config['output']['directory'] / 'QGIS'
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / 'swan_grids.png', dpi=300, bbox_inches='tight')
    print(f"\nFigure saved to: {output_dir / 'swan_grids.png'}")
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create and plot SWAN grids')
    parser.add_argument('--config', required=True, help='Path to the configuration file')
    args = parser.parse_args()
    
    main(args.config)
