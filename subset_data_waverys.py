"""
This script is used to subset the data from the WAverys dataset.
"""


import sys
sys.path.append('functions/')
import xarray as xr
import pandas as pd
from pathlib import Path
from functions.mdatetime import datenum
from functions.data import generate_filename
import numpy as np


def subset_data_by_area(ds, start_date, end_date, lon_min, lon_max, lat_min, lat_max):
    """
    Subset the data from the Waverys dataset.
    """

    if 'valid_time' in ds.variables:
        time_var = 'valid_time'
        print('Variable valid_time found')
    elif 'time' in ds.variables:
        time_var = 'time'
        print('Variable time found')
    else:
        raise KeyError("No se encontró ni 'valid_time' ni 'time' en el archivo NetCDF.")
    
    units = ds[time_var].attrs.get("units", "")
    if not units:
        # Intentar con netCDF4
        try:
            from netCDF4 import Dataset
            # Obtener el path del archivo original
            nc_path = ds.encoding.get('source', None)
            if nc_path is None:
                raise ValueError("No se puede determinar el path del archivo NetCDF para netCDF4.")
            nc4ds = Dataset(nc_path)
            if time_var in nc4ds.variables:
                nc4var = nc4ds.variables[time_var]
                if 'units' in nc4var.ncattrs():
                    units = nc4var.getncattr('units')
                else:
                    raise ValueError(f"No se encontró el atributo 'units' en la variable '{time_var}' usando netCDF4.")
            else:
                raise ValueError(f"No se encontró la variable '{time_var}' en el archivo NetCDF usando netCDF4.")
            nc4ds.close()
        except ImportError:
            raise ImportError("No se pudo importar netCDF4. Instala el paquete para soporte extendido de atributos.")
    if "seconds since 1970" in units:
        t = ds[time_var][:].data.astype(float) + datenum((1970, 1, 1))
    elif "hours since 1900" in units:
        t = ds[time_var][:].data.astype(float) * 3600 + datenum((1900, 1, 1))
    else:
        raise ValueError(f"Unsupported or missing time units: '{units}' in variable '{time_var}'")
    
    return ds.sel({time_var: slice(start_date, end_date), "longitude": slice(lon_min, lon_max), "latitude": slice(lat_min, lat_max)})

def main(config):
    """
    Main function to run the script.
    """

    start_date = config['analysis']['dates']['start']
    end_date = config['analysis']['dates']['end']
    regional_grid = config['grids']['regional'] 
    lon_min = regional_grid['bounds']['lon_min']
    lon_max = regional_grid['bounds']['lon_max']
    lat_min = regional_grid['bounds']['lat_min']
    lat_max = regional_grid['bounds']['lat_max']

    source_path = Path(config['source']['path'])

    # --- WAVE ---
    ds_wave = xr.open_dataset(source_path / 'cmems_mod_glo_wav_myint_0.2deg_PT3H-i_1747134522272.nc')
    ds_wave_subset = subset_data_by_area(ds_wave, start_date, end_date, lon_min, lon_max, lat_min, lat_max)
    filename_wave = generate_filename(
        'wave_data', 'hourly', start_date, end_date,
        {'lon_min': lon_min, 'lon_max': lon_max, 'lat_min': lat_min, 'lat_max': lat_max}
    )
    output_dir_wave = Path(config['base']['path']) / config['output']['directory'] / config['output']['data']['wave']
    output_dir_wave.mkdir(parents=True, exist_ok=True)
    output_path_wave = output_dir_wave / filename_wave
    print(f"Saving subsetted wave data to {output_path_wave}")
    ds_wave_subset.to_netcdf(output_path_wave)

    # --- WIND ---
    ds_wind = xr.open_dataset(source_path / 'cmems_obs-wind_glo_phy_nrt_l4_0.125deg_PT1H_1747134907346.nc')
    ds_wind_subset = subset_data_by_area(ds_wind, start_date, end_date, lon_min, lon_max, lat_min, lat_max)
    filename_wind = generate_filename(
        'wind_data', 'hourly', start_date, end_date,
        {'lon_min': lon_min, 'lon_max': lon_max, 'lat_min': lat_min, 'lat_max': lat_max}
    )
    output_dir_wind = Path(config['base']['path']) / config['output']['directory'] / config['output']['data']['wind']
    output_dir_wind.mkdir(parents=True, exist_ok=True)
    output_path_wind = output_dir_wind / filename_wind
    print(f"Saving subsetted wind data to {output_path_wind}")
    encoding = {}
    for var in ds_wind_subset.data_vars:
        encoding[var] = {'_FillValue': np.float32(-32767)}

    # Check if the subset is empty before saving
    if any(dim == 0 for dim in ds_wind_subset.dims.values()):
        print("WARNING: Wind subset is empty for the specified lat/lon/time range. File will not be saved.")
    else:
        ds_wind_subset.to_netcdf(output_path_wind, encoding=encoding)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Run Waverys data subsetter.')
    parser.add_argument('--config', type=str, required=True, help='YAML config file for one run.')
    args = parser.parse_args()

    main(args.config)

