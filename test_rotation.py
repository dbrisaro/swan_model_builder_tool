import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, radians

def create_rotated_grid(lon_min, lon_max, lat_min, lat_max, dx, dy, rotation=0.0):
    """
    Creates a grid of points and its rotated version.
    
    Parameters:
        lon_min, lon_max: Longitude bounds
        lat_min, lat_max: Latitude bounds
        dx, dy: Grid spacing in degrees
        rotation: Rotation angle in degrees
    
    Returns:
        original_points: List of (lon, lat) tuples for original grid
        rotated_points: List of (lon, lat) tuples for rotated grid
        center: (lon, lat) tuple for center point
    """
    # Calculate center point
    center_lon = (lon_min + lon_max) / 2
    center_lat = (lat_min + lat_max) / 2
    
    # Create original grid points
    lons = np.arange(lon_min, lon_max + dx, dx)
    lats = np.arange(lat_min, lat_max + dy, dy)
    
    # Create meshgrid of points
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    original_points = list(zip(lon_grid.flatten(), lat_grid.flatten()))
    
    # Convert rotation to radians
    angle_rad = radians(rotation)
    
    # Rotate each point around the center
    rotated_points = []
    for lon, lat in original_points:
        # Translate to origin
        x = lon - center_lon
        y = lat - center_lat
        
        # Rotate
        x_rot = x * cos(angle_rad) - y * sin(angle_rad)
        y_rot = x * sin(angle_rad) + y * cos(angle_rad)
        
        # Translate back
        new_lon = x_rot + center_lon
        new_lat = y_rot + center_lat
        
        rotated_points.append((new_lon, new_lat))
    
    return original_points, rotated_points, (center_lon, center_lat)

def plot_grids(original_points, rotated_points, center):
    """Plot both original and rotated grids."""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Plot original grid
    x_orig, y_orig = zip(*original_points)
    ax.plot(x_orig, y_orig, 'bo', markersize=2, alpha=0.5, label='Original Grid')
    
    # Plot rotated grid
    x_rot, y_rot = zip(*rotated_points)
    ax.plot(x_rot, y_rot, 'ro', markersize=2, alpha=0.5, label='Rotated Grid')
    
    # Plot center point
    ax.plot(center[0], center[1], 'ko', label='Center')
    
    # Print some example points
    print("\nExample original grid points:")
    for i, (lon, lat) in enumerate(original_points[:5]):
        print(f"  Point {i+1}: ({lon:.4f}, {lat:.4f})")
    
    print("\nExample rotated grid points:")
    for i, (lon, lat) in enumerate(rotated_points[:5]):
        print(f"  Point {i+1}: ({lon:.4f}, {lat:.4f})")
    
    print(f"\nCenter point: ({center[0]:.4f}, {center[1]:.4f})")
    print(f"\nTotal number of grid points: {len(original_points)}")
    
    ax.set_aspect('equal')
    ax.grid(True)
    ax.legend()
    plt.savefig('rotation_test.png')
    plt.close()

# Example coordinates from Peru transition grid
lon_min, lon_max = -80.5, -75.5
lat_min, lat_max = -14.5, -11.0
dx, dy = 0.5, 0.5  # Grid spacing in degrees
rotation = 45.0

# Create and plot grids
original_points, rotated_points, center = create_rotated_grid(
    lon_min, lon_max, lat_min, lat_max, dx, dy, rotation
)
plot_grids(original_points, rotated_points, center) 