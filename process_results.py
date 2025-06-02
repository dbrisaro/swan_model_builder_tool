#!/usr/bin/env python3
"""Process SWAN results from .mat files"""

import yaml
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
from other_functions import convertMat2Nc
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def check_nc_file(nc_file):
    """Print information about the NetCDF file contents"""
    print(f"\nChecking file: {nc_file}")
    
    # Open dataset
    ds = xr.open_dataset(nc_file)
    
    # Print basic information
    print("\nDimensions:")
    for dim_name, size in ds.sizes.items():
        print(f"  {dim_name}: {size}")
    
    print("\nVariables:")
    for var in ds.variables:
        print(f"  {var}:")
        print(f"    shape: {ds[var].shape}")
        print(f"    dtype: {ds[var].dtype}")
        if 'units' in ds[var].attrs:
            print(f"    units: {ds[var].attrs['units']}")
    
    # Check if Xp and Yp are 2D (which they should be for rotated grids)
    print("\nGrid structure:")
    print(f"  Xp dimensions: {ds.Xp.dims}")
    print(f"  Yp dimensions: {ds.Yp.dims}")
    
    # Plot first timestep of Hsig with coordinates to check rotation
    plt.figure(figsize=(10, 8))
    plt.pcolormesh(ds.Xp, ds.Yp, ds.Hsig[0])
    plt.colorbar(label='Hsig [m]')
    plt.title(f'First timestep Hsig - {Path(nc_file).stem}')
    check_png_path = Path(nc_file).parent / f'check_{Path(nc_file).stem}.png'
    plt.savefig(check_png_path)
    plt.close()
    
    ds.close()

def process_mat_files(results_dir, output_dir):
    """Convert .mat files to .nc format"""
    # Create output directory if it doesn't exist
    output_dir = Path(output_dir) / 'files'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Process each .mat file
    for mat_file in results_dir.glob('*.mat'):
        nc_file = output_dir / f"{mat_file.stem}.nc"
        print(f"Converting {mat_file} to {nc_file} (rotated: True)")
        convertMat2Nc(str(mat_file), str(nc_file), isRotated=True, isNonSpherical=True)
        
        # Check the generated NC file
        check_nc_file(nc_file)

def plot_daily_mean_hsig(nc_file, output_dir, start_date, end_date):
    """Generate maps of daily mean significant wave height."""
    # Open the dataset with xarray
    ds = xr.open_dataset(nc_file)
    
    # Convert time to datetime64
    ds['time'] = pd.to_datetime(ds.time.values, unit='s')
    
    # Select time period
    start = pd.to_datetime(start_date, dayfirst=True)
    end = pd.to_datetime(end_date)
    ds_period = ds.sel(time=slice(start, end))
    
    # Check if Hsig exists
    if 'Hsig' not in ds_period:
        ds.close()
        return
    # Calculate daily mean
    daily_mean = ds_period['Hsig'].resample(time='1D').mean()
    
    # Get global min and max for colorbar limits
    vmin = daily_mean.min().values
    vmax = daily_mean.max().values
    
    # Create maps directory
    maps_dir = Path(output_dir) / 'maps' / f"{start_date}_to_{end_date}"
    maps_dir.mkdir(parents=True, exist_ok=True)
    
    # Calculate number of rows and columns for subplots
    n_days = len(daily_mean.time)
    n_cols = 4  # Increased number of columns
    n_rows = (n_days + n_cols - 1) // n_cols  # Calculate number of rows needed
    
    # Create a figure with subplots in a grid
    fig = plt.figure(figsize=(20, 4*n_rows))
    gs = plt.GridSpec(n_rows, n_cols, figure=fig)
    
    # Create a list to store all the plots
    plots = []
    
    # Plot each day
    for i, day in enumerate(daily_mean.time):
        row = i // n_cols
        col = i % n_cols
        ax = fig.add_subplot(gs[row, col], projection=ccrs.PlateCarree())
        
        # Get data for this day
        data = daily_mean.sel(time=day).values
        
        # For rotated grids, Xp and Yp are 2D arrays
        x = ds.Xp.values
        y = ds.Yp.values
        
        # Plot data with fixed colorbar limits
        im = ax.pcolormesh(x, y, data,
                          transform=ccrs.PlateCarree(),
                          cmap='viridis',
                          vmin=vmin,
                          vmax=vmax,
                          shading='nearest')  # Changed to 'nearest' to handle same dimensions
        plots.append(im)
        
        # Add map features - LAND should be solid color to mask data
        ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=100)
        ax.add_feature(cfeature.COASTLINE, zorder=101)
        ax.add_feature(cfeature.BORDERS, linestyle=':', zorder=102)
        
        # Add lat/lon labels
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        
        # Add title aligned to the left
        date_str = pd.to_datetime(day.values).strftime('%Y-%m-%d')
        ax.set_title(f'Daily Mean Hsig - {date_str}', loc='left', pad=10)
    
    # Add a single colorbar for all plots
    cbar_ax = fig.add_axes([0.25, 0.02, 0.5, 0.02])  # [left, bottom, width, height]
    cbar = fig.colorbar(plots[0], cax=cbar_ax, orientation='horizontal')
    cbar.set_label('Hsig [m]')
    
    # Adjust layout and save
    grid_name = Path(nc_file).stem.split('_')[0]
    output_path = maps_dir / f'{grid_name}_daily_mean_hsig.png'
    plt.savefig(output_path,
               bbox_inches='tight', dpi=300)
    plt.close()
    
    # Close the dataset
    ds.close()

def extract_time_series(nc_file, output_dir, start_date, end_date, points):
    """Extract and plot time series at specific points."""
    # Open the dataset with xarray
    ds = xr.open_dataset(nc_file)
    
    # Convert time to datetime64
    ds['time'] = pd.to_datetime(ds.time.values, unit='s')
    
    # Select time period
    start = pd.to_datetime(start_date)
    end = pd.to_datetime(end_date)
    ds_period = ds.sel(time=slice(start, end))
    
    # Create time series directory
    ts_dir = Path(output_dir) / 'time_series' / f"{start_date}_to_{end_date}"
    ts_dir.mkdir(parents=True, exist_ok=True)
    
    # Get coordinates (Xp y Yp ya son 2D)
    lon_grid = ds['Xp'].values  # (4, 5)
    lat_grid = ds['Yp'].values  # (4, 5)

    # Create figure for time series
    fig, ax = plt.subplots(figsize=(12, 6))

    # Extract and plot time series for each point
    for point in points:
        # Find nearest grid point
        dist = np.sqrt((lon_grid - point['lon'])**2 + (lat_grid - point['lat'])**2)
        idx = np.unravel_index(np.argmin(dist), dist.shape)
        print(f"Point: {point['name']}, idx: {idx}, idx[0] max: {lon_grid.shape[0]-1}, idx[1] max: {lon_grid.shape[1]-1}")
        # Extract time series (ajustado el orden de los índices)
        ts = ds_period['Hsig'].isel(Yp=idx[0], Xp=idx[1])
        # Plot time series
        ax.plot(ts.time, ts.values, label=point['name'])
    
    # Customize plot
    ax.set_xlabel('Time')
    ax.set_ylabel('Hsig [m]')
    ax.set_title('Significant Wave Height Time Series')
    ax.legend()
    ax.grid(True)
    
    # Save figure
    grid_name = Path(nc_file).stem.split('_')[0]
    plt.savefig(ts_dir / f'{grid_name}_time_series.png',
               bbox_inches='tight', dpi=300)
    plt.close()
    
    # Close the dataset
    ds.close()

def main(config):
    # Read dates and points from config['analysis']
    start_date = config['analysis']['dates']['start']
    end_date = config['analysis']['dates']['end']
    points = config['analysis']['points']
    
    # Set up paths using base_path from config
    base_path = Path(config['base']['path'])
    results_dir = base_path / config['output']['directory'] / 'SWAN/04_results'
    output_dir = base_path / config['output']['directory'] / 'SWAN/05_processed_results'
    files_dir = output_dir / 'files'
    
    # Process results
    process_mat_files(results_dir, output_dir)
    
    # Generate maps and time series for each NetCDF file
    for nc_file in files_dir.glob('*.nc'):
        
        # Calculate number of days between start and end dates
        days_diff = (pd.to_datetime(end_date) - pd.to_datetime(start_date)).days
        if days_diff < 20:
            plot_daily_mean_hsig(str(nc_file), output_dir, start_date, end_date)
        extract_time_series(str(nc_file), output_dir, start_date, end_date, points)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run SWAN results processor.')
    parser.add_argument('--config', type=str, required=True, help='YAML config file for one run.')
    args = parser.parse_args()

    main(args.config)
