"""Request bathymetry data from ETOPO"""

import os
import sys
import yaml
import json
import numpy as np
from pathlib import Path
from erddapy import ERDDAP
import rasterio
from functions.data import generate_date_lists, generate_filename
import requests

    
def main(config):
    # Get parameters from config
    time_config = config['time']
    start_date = time_config['start']
    end_date = time_config['end']
    frequency = time_config['frequency']

    # Get bounds from regional grid
    regional_grid = config['grids']['regional']['bounds']

    # Base request parameters for bathymetry
    e = ERDDAP(
        server="https://oceanwatch.pifsc.noaa.gov/erddap",
        protocol="griddap",
    )

    dataset = "ETOPO_2022_v1_15s"
    e.dataset_id = dataset

    # Definir variable y slices para griddap
    variable = "z"
    lat_min = regional_grid['lat_min']
    lat_max = regional_grid['lat_max']
    lon_min = regional_grid['lon_min']
    lon_max = regional_grid['lon_max']

    # Funciones para convertir lat/lon a índices de grilla ETOPO 2022 v1 15s
    # Resolución: 1/240 grado = 0.0041666667 grados

    def lat_to_index(lat):
        return int(round((lat + 90) / 0.0041666667))

    def lon_to_index(lon):
        # Si lon < 0, conviértela a 0-360
        if lon < 0:
            lon = lon + 360
        return int(round(lon / 0.0041666667))

    # Usar los valores de la config para definir el área
    lat_min_idx = lat_to_index(lat_min)
    lat_max_idx = lat_to_index(lat_max)
    lon_min_idx = lon_to_index(lon_min)
    lon_max_idx = lon_to_index(lon_max)

    # Construir la query usando índices
    url = (
        f"https://oceanwatch.pifsc.noaa.gov/erddap/griddap/{dataset}.nc?"
        f"{variable}[{lat_min_idx}:{lat_max_idx}][{lon_min_idx}:{lon_max_idx}]"
    )

    base_path = Path(config['base']['path'])
    output_dir = base_path / config['output']['directory'] / config['output']['data']['bathy']
    output_dir.mkdir(parents=True, exist_ok=True)
    filename = output_dir / generate_filename('bathy_data', frequency, start_date, end_date, 
        {'lon_min': lon_min, 'lon_max': lon_max, 'lat_min': lat_min, 'lat_max': lat_max}
    )
    filename = filename.with_suffix('.nc')  # Guardar como .nc


    # Print request details
    print("\nETOPO Request Details:")
    print("===================")
    print(f"Dataset: {dataset}")
    print("\nVariables:", variable)
    print("Latitud:", lat_min, lat_max)
    print("Longitud:", lon_min, lon_max)
    print(f"\nOutput file: {filename}")
    print(f"\nDownload URL: {url}")
    
    # Descargar el archivo manualmente
    print("\nDescargando datos...")
    r = requests.get(url)
    with open(filename, "wb") as f:
        f.write(r.content)
    print("Descarga completada.")

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run SWAN config generator.')
    parser.add_argument('--config', type=str, required=True, help='Path to the experiment specifications file')
    args = parser.parse_args()

    main(args.config)