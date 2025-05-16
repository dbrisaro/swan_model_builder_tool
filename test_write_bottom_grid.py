import os
import numpy as np
import xarray as xr
import rasterio
from rasterio.transform import from_origin
from pathlib import Path
from scipy.interpolate import griddata

def write_bottom_grid(inFiles, outFile, x, y, faces=None):
    """
    Inspecciona archivos raster (GeoTIFF o NetCDF) y escribe un archivo de texto con la batimetría interpolada en los puntos (x, y).
    Si el archivo es .nc, interpola directamente desde el NetCDF (sin convertir a .tif).
    Soporta longitudes en 0-360 o -180 a 180 automáticamente.
    """
    zp0 = np.zeros((x.size,), dtype=np.float64)*np.nan
    for inFile in inFiles:
        if str(inFile).endswith('.nc'):
            ds = xr.open_dataset(inFile)
            # Usar los nombres correctos: 'z', 'longitude', 'latitude'
            data = ds['z'].values
            lon = ds['longitude'].values
            lat = ds['latitude'].values
            # --- Manejo de sistemas de longitud ---
            try:
                # Si las longitudes del archivo están en 0-360 y x está en -180 a 180, convierte x a 0-360
                if lon.min() > 0 and np.any(x < 0):
                    x_interp = x.copy()
                    x_interp[x_interp < 0] += 360
                # Si las longitudes del archivo están en -180 a 180 y x está en 0-360, convierte x a -180 a 180
                elif lon.min() < 0 and np.any(x > 180):
                    x_interp = x.copy()
                    x_interp[x_interp > 180] -= 360
                else:
                    x_interp = x
            except Exception as e:
                x_interp = x
            # Si data es 2D (lat, lon), hacer meshgrid
            if data.ndim == 2:
                lon2d, lat2d = np.meshgrid(lon, lat)
                points = np.column_stack((lon2d.flatten(), lat2d.flatten()))
                values = data.flatten()
            else:
                points = np.column_stack((lon, lat))
                values = data
            # Interpola en los puntos de la grilla
            interp_vals = griddata(points, values, (x_interp.flatten(), y.flatten()), method='linear')
            zp1 = interp_vals
        else:
            with rasterio.open(inFile) as src:
                vals = []
                for xi, yi in zip(x.flatten(), y.flatten()):
                    for val in src.sample([(xi, yi)]):
                        vals.append(val[0])
                zp1 = np.array(vals)
        zp0[np.isnan(zp1) == 0] = zp1[np.isnan(zp1) == 0]

    # Escribir el resultado a un archivo de texto
    zp0 = np.reshape(zp0, x.shape)
    with open(os.path.abspath(outFile), 'w') as f:
        np.savetxt(f, -1*zp0, '%8.2f', delimiter=' ')

# Ejemplo de uso:
if __name__ == "__main__":
    import xarray as xr
    ds = xr.open_dataset('/Users/daniela/Documents/swan/swan_experiments/run_peninsula_paracas/DATA/BATHY/bathy_data_hourly_2024-12-15_2025-01-15_-76.5_-76.3_-13.9_-13.7.nc')
    x = ds['longitude'].values
    y = ds['latitude'].values
    xx, yy = np.meshgrid(x, y)
    write_bottom_grid([
        '/Users/daniela/Documents/swan/swan_experiments/run_peninsula_paracas/DATA/BATHY/bathy_data_hourly_2024-12-15_2025-01-15_-76.5_-76.3_-13.9_-13.7.nc'
    ], 'bottom_test.txt', xx, yy) 