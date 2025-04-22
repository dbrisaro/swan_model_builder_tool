"""Functions for grid operations"""

import numpy as np
import geopandas as gpd
from shapely.geometry import LineString

def create_rectangular_grid(lon_min, lon_max, lat_min, lat_max, x_len, y_len, name='REGIONAL', rotation=0.0):
    """Create a rectangular grid for SWAN model
    
    Parameters
    ----------
    lon_min, lon_max : float
        Minimum and maximum longitude
    lat_min, lat_max : float
        Minimum and maximum latitude
    x_len, y_len : int
        Number of grid points in x and y directions
    name : str, optional
        Name of the grid, by default 'REGIONAL'
    rotation : float, optional
        Grid rotation in degrees, by default 0.0
    
    Returns
    -------
    gpd.GeoDataFrame
        Grid as a GeoDataFrame
    """
    # Calculate grid spacing
    dx = (lon_max - lon_min) / (x_len - 1)
    dy = (lat_max - lat_min) / (y_len - 1)
    
    # Create grid coordinates
    x = np.linspace(lon_min, lon_max, x_len)
    y = np.linspace(lat_min, lat_max, y_len)
    
    # Create grid lines
    lines = []
    for i in range(x_len):
        lines.append(LineString([(x[i], lat_min), (x[i], lat_max)]))
    for j in range(y_len):
        lines.append(LineString([(lon_min, y[j]), (lon_max, y[j])]))
    
    # Create GeoDataFrame
    grid = gpd.GeoDataFrame(geometry=lines)
    grid['name'] = name
    grid['rotation'] = rotation
    
    return grid

def create_grid_lines(grid_params):
    """Create grid lines from grid parameters
    
    Parameters
    ----------
    grid_params : dict
        Grid parameters including:
        - lon_min, lon_max : float
            Minimum and maximum longitude
        - lat_min, lat_max : float
            Minimum and maximum latitude
        - dx, dy : float
            Grid spacing in x and y directions
        - name : str
            Name of the grid
        - rotation : float
            Grid rotation in degrees
    
    Returns
    -------
    gpd.GeoDataFrame
        Grid lines as a GeoDataFrame
    """
    # Calculate number of cells
    x_len = int((grid_params['lon_max'] - grid_params['lon_min']) / grid_params['dx']) + 1
    y_len = int((grid_params['lat_max'] - grid_params['lat_min']) / grid_params['dy']) + 1
    
    # Create grid coordinates
    x = np.linspace(grid_params['lon_min'], grid_params['lon_max'], x_len)
    y = np.linspace(grid_params['lat_min'], grid_params['lat_max'], y_len)
    
    # Create grid lines
    lines = []
    
    # Create horizontal lines
    for j in range(y_len):
        line = LineString([(x[0], y[j]), (x[-1], y[j])])
        lines.append(line)
    
    # Create vertical lines
    for i in range(x_len):
        line = LineString([(x[i], y[0]), (x[i], y[-1])])
        lines.append(line)
    
    # Create GeoDataFrame
    grid = gpd.GeoDataFrame(geometry=lines, crs='EPSG:4326')
    grid['name'] = grid_params['name']
    grid['rotation'] = grid_params['rotation']
    
    return grid 