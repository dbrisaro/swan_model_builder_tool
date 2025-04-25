#!/usr/bin/env python3
"""Process SWAN results from .mat files"""

import os
import sys
import yaml
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from netCDF4 import Dataset
from other_functions import convertMat2Nc, extractSwanTs
from functions.mdatetime import *
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def read_config(config_file):
    """Read experiment configuration from YAML file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def process_mat_files(results_dir, output_dir):
    """Process all .mat files in the results directory."""
    # Create output directory if it doesn't exist
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Process each .mat file
    for mat_file in Path(results_dir).glob('*.mat'):
        print(f"Processing {mat_file.name}...")
        
        # Convert to NetCDF
        nc_file = output_dir / f"{mat_file.stem}.nc"
        convertMat2Nc(str(mat_file), str(nc_file), isRotated=False, isNonSpherical=False)
        
        # Extract time series at specific points (if points file exists)
        points_file = Path('points.csv')  # You should create this file with your points of interest
        if points_file.exists():
            ts_file = output_dir / f"{mat_file.stem}_timeseries.csv"
            extractSwanTs(str(nc_file), str(ts_file), str(points_file))

def calculate_statistics(nc_file):
    """Calculate basic statistics from NetCDF file."""
    with Dataset(nc_file, 'r') as nc:
        stats = {}
        for var in nc.variables:
            if nc[var].ndim == 3:  # Time-varying variables
                data = nc[var][:]
                stats[var] = {
                    'mean': np.nanmean(data),
                    'max': np.nanmax(data),
                    'min': np.nanmin(data),
                    'std': np.nanstd(data)
                }
        return stats

def plot_results(nc_file, output_dir):
    """Generate basic plots from NetCDF file."""
    with Dataset(nc_file, 'r') as nc:
        # Create plots directory
        plots_dir = Path(output_dir) / 'plots'
        plots_dir.mkdir(exist_ok=True)
        
        # Plot time series for each variable
        for var in nc.variables:
            if nc[var].ndim == 3:  # Time-varying variables
                plt.figure(figsize=(10, 6))
                data = nc[var][:]
                time = nc['time'][:]
                
                # Plot mean over space
                plt.plot(time, np.nanmean(data, axis=(1, 2)))
                plt.title(f'{var} Time Series')
                plt.xlabel('Time')
                plt.ylabel(nc[var].units)
                
                # Save plot
                plt.savefig(plots_dir / f'{var}_timeseries.png')
                plt.close()

def plot_maps(nc_file, output_dir, time_indices=None):
    """Generate maps for each variable at different times."""
    with Dataset(nc_file, 'r') as nc:
        # Create maps directory
        maps_dir = Path(output_dir) / 'maps'
        maps_dir.mkdir(exist_ok=True)
        
        # Get coordinates
        lon = nc['longitude'][:]
        lat = nc['latitude'][:]
        
        # If no specific time indices provided, use first, middle and last time step
        if time_indices is None:
            n_times = len(nc['time'])
            time_indices = [0, n_times//2, -1]
        
        # Plot maps for each variable
        for var in nc.variables:
            if nc[var].ndim == 3:  # Time-varying variables
                data = nc[var][:]
                
                for t_idx in time_indices:
                    # Create figure with cartopy projection
                    fig = plt.figure(figsize=(12, 8))
                    ax = plt.axes(projection=ccrs.PlateCarree())
                    
                    # Add map features
                    ax.add_feature(cfeature.LAND)
                    ax.add_feature(cfeature.COASTLINE)
                    ax.add_feature(cfeature.BORDERS, linestyle=':')
                    
                    # Plot data
                    im = ax.pcolormesh(lon, lat, data[t_idx], 
                                     transform=ccrs.PlateCarree(),
                                     cmap='viridis')
                    
                    # Add colorbar
                    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05)
                    cbar.set_label(nc[var].units)
                    
                    # Add title and labels
                    time_str = nc['time'][t_idx]
                    plt.title(f'{var} at time {time_str}')
                    ax.set_xlabel('Longitude')
                    ax.set_ylabel('Latitude')
                    
                    # Save plot
                    plt.savefig(maps_dir / f'{var}_map_t{t_idx}.png', 
                              bbox_inches='tight', dpi=300)
                    plt.close()

def main():
    # Read configuration
    config_file = Path('/Users/daniela/Documents/swan/swan_experiments/experiments_specs.txt')
    config = read_config(config_file)
    
    # Set up paths
    base_dir = Path('/Users/daniela/Documents/swan/swan_experiments')
    results_dir = base_dir / config['output']['directory'] / 'SWAN/04_results'
    output_dir = base_dir / config['output']['directory'] / 'processed_results'
    
    # Process results
    process_mat_files(results_dir, output_dir)
    
    # Calculate statistics and generate plots for each NetCDF file
    for nc_file in output_dir.glob('*.nc'):
        stats = calculate_statistics(str(nc_file))
        print(f"\nStatistics for {nc_file.name}:")
        for var, values in stats.items():
            print(f"{var}:")
            for stat, value in values.items():
                print(f"  {stat}: {value:.2f}")
        
        # Generate time series plots
        plot_results(str(nc_file), output_dir)
        
        # Generate maps
        plot_maps(str(nc_file), output_dir)

if __name__ == '__main__':
    main()
