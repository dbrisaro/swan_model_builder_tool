import sys
sys.path.append('functions/')

import os
import inspect
import configparser

from pathlib import Path

import numpy as np
import xarray as xr
import geopandas as gpd
import rasterio
from rasterio.transform import from_origin

from scipy.io import loadmat
from functions.spatial import *
from functions.mdatetime import *




## 
def autoNestModels(models):
    """Detect nested models and apply default naming"""

    # get areas into an array of arrays
    areas = np.array([model.swanGrid.getBoundary() for model in models])

    # Get level of nesting for each model
    na = areas.shape[0]  # number of areas
    level = np.zeros((na,), dtype=np.int16)  # level of nesting
    inside = np.zeros((na, na), dtype=bool)  # in logical index
    for aa in range(na):
        for bb in range(na):
            if aa == bb:
                continue
            # is aa in bb?
            lgi = inPolygon(areas[aa], areas[bb])
            # if aa is in bb
            if np.all(lgi):
                level[aa] = level[aa] + 1
                inside[aa, bb] = True

    # Sort the areas based on level
    index = np.argsort(level)
    level = level[index]
    areas = areas[index]
    inside = inside[index, :][:, index]

    models = list(np.array(models)[index])

    # Shuffle models, rename and apply nesting
    for aa in range(na):
        model = models[aa]
        if level[aa] > 0:
            # If is nested find parent
            ii = np.max(np.where(inside[aa])[0])
            model.setParent(models[ii])

## 
from shapely.geometry import Polygon
import pandas as pd
import geopandas as gpd

'''
def create_rectangular_grid(lon_min, lon_max, lat_min, lat_max, x_len, y_len, name='REGIONAL', rotation=0.0):
    """
    Creates a single rectangle as a GeoDataFrame with metadata and a rotated grid of points.

    Parameters:
        lon_min (float): Minimum longitude
        lon_max (float): Maximum longitude
        lat_min (float): Minimum latitude
        lat_max (float): Maximum latitude
        x_len (float): X resolution in degrees
        y_len (float): Y resolution in degrees
        name (str): Name for the region
        rotation (float): Rotation angle in degrees

    Returns:
        gpd.GeoDataFrame: A GeoDataFrame with one polygon and metadata
    """
    from math import cos, sin, radians

    # Calculate center point
    center_lon = (lon_min + lon_max) / 2
    center_lat = (lat_min + lat_max) / 2

    # Create original grid points
    lons = np.arange(lon_min, lon_max + x_len, x_len)
    lats = np.arange(lat_min, lat_max + y_len, y_len)
    
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

    # Create polygon from rotated points
    polygon = Polygon(rotated_points)

    # Create DataFrame with metadata and grid points
    df = pd.DataFrame({
        "Name": [name],
        "Rotation": [rotation],
        "X Length": [x_len],
        "Y Length": [y_len],
        "geometry": [polygon],
        "grid_points": [rotated_points]  # Store the rotated grid points
    })

    gdf = gpd.GeoDataFrame(df, geometry='geometry', crs="EPSG:4326")

    return gdf
'''
def create_rectangular_grid(lon_min, lon_max, lat_min, lat_max, x_len, y_len, name='REGIONAL', rotation=0.0):
    """
    Creates a single rectangle as a GeoDataFrame with metadata and a rotated grid of points.

    Parameters:
        lon_min (float): Minimum longitude
        lon_max (float): Maximum longitude
        lat_min (float): Minimum latitude
        lat_max (float): Maximum latitude
        x_len (float): X resolution in degrees
        y_len (float): Y resolution in degrees
        name (str): Name for the region
        rotation (float): Rotation angle in degrees

    Returns:
        gpd.GeoDataFrame: A GeoDataFrame with one polygon and metadata
    """
    from math import cos, sin, radians

    # Calculate center point
    center_lon = (lon_min + lon_max) / 2
    center_lat = (lat_min + lat_max) / 2
    
    # Calculate n_points
    n_lon = int(np.ceil((lon_max-lon_min)/x_len))
    n_lat = int(np.ceil((lat_max-lat_min)/y_len))

    # Create original grid points
    lons = np.linspace(lon_min, lon_max, n_lon)
    lats = np.linspace(lat_min, lat_max, n_lat)
    
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

    # Create polygon from rotated points
    polygon = Polygon(rotated_points)

    # Create DataFrame with metadata and grid points
    df = pd.DataFrame({
        "Name": [name],
        "Rotation": [rotation],
        "X Length": [x_len],
        "Y Length": [y_len],
        "geometry": [polygon],
        "grid_points": [rotated_points]  # Store the rotated grid points
    })

    gdf = gpd.GeoDataFrame(df, geometry='geometry', crs="EPSG:4326")

    return gdf

##

def readGridsFromFile(file):

    """
    A function for building SwanGrid objects from a vector file with requisite attributes.

    PARAMETERS
    ----------
    file : string
        File path of vector file.

    RETURNS
    -------
    grids : list
        List of grid objects generated from vector layer feature attributes and geometry

    """

    # create dummy data source handle
    ds = None

    # find driver and get data source handle
    for ii in range(ogr.GetDriverCount()):
        driver = ogr.GetDriver(ii)
        
        try: 
            ds = driver.Open(file)
        except RuntimeError:
            ds = None

        # break loop when found
        if ds is not None:
            break

    # check data source is valid
    if ds is None:
        raise ValueError('Invalid file format')

    # get layer and srs handle
    layer = ds.GetLayer()
    srs = layer.GetSpatialRef()

    # check geometry type for polygon or polygon collection
    if (layer.GetGeomType() != 3) and (layer.GetGeomType() != 6):
        raise ValueError('Layer not polygon')

    # create empty list to store grids
    grids = list()

    # iterate over features
    for feature in layer:
        # get field data from the feature
        name = feature.GetField('Name')
        rotation = feature.GetField('Rotation')
        dx = feature.GetField('X Length')
        dy = feature.GetField('Y Length')

        # iterate down to single polygon
        geometry = feature.GetGeometryRef()
        while geometry.GetPointCount() == 0:
            geometry = geometry.GetGeometryRef(0)

        # get polygon as numpy array
        polygon = np.array(geometry.GetPoints())

        # create a new grid object and append to the list
        grids.append(SwanGrid(name, polygon, rotation, dx, dy, srs=srs))

    return grids

## 
# NEW FUNCTION WITH GEOPANDAS
def read_grids_from_file_gpd(file):
    # Read the vector file into a GeoDataFrame
    gdf = gpd.read_file(file)

    # Check geometry type
    if not all(gdf.geometry.geom_type.isin(['Polygon', 'MultiPolygon'])):
        raise ValueError("Layer not polygon")

    grids = []

    for _, row in gdf.iterrows():
        name = row['Name']
        rotation = row['Rotation']
        dx = row['X Length']
        dy = row['Y Length']
        
        # Convert geometry to single polygon and extract coordinates
        geom = row.geometry
        if geom.geom_type == 'MultiPolygon':
            geom = list(geom.geoms)[0]  # or use union/cascaded_union depending on your logic

        polygon = np.array(geom.exterior.coords)

        # Assuming SwanGrid is defined elsewhere
        grids.append(SwanGrid(name, polygon, rotation, dx, dy, srs=gdf.crs))

    return grids

##

def writeWindTsBc(inFile, outFile, swanModel):
    # read the wind time series
    with open(inFile, 'r') as f:
        # skip header
        f.readline()
        f.readline()

        # get the data
        names, fmts = ('TIME', 'U10', 'V10'), ('U19', 'f4', 'f4')
        data = np.loadtxt(f, delimiter=',', dtype={'names': names, 'formats': fmts})

        t, u, v = datenum(data['TIME']), data['U10'], data['V10']

        # get SWAN INPGRID and READINP parameters
        x = swanModel.swanGrid.getX()
        y = swanModel.swanGrid.getY()

        x1, x2 = x.min(), x.max()
        y1, y2 = y.min(), y.max()

        x0, nx, dx = x1, 2, x2 - x1
        y0, ny, dy = y1, 2, y2 - y1

        t1, t2 = datenum(swanModel.timeStart), datenum(swanModel.timeEnd)

        # find indices within model limits, if anything outside limits, return.
        try:
            it = np.arange(np.where(t <= t1)[0][-1], np.where(t2 <= t)[0][0] + 1)
        except IndexError:
            return None

        ts = datestr(t[0], '%Y%m%d.%H%M%S')
        tf = datestr(t[it[-1]], '%Y%m%d.%H%M%S')
        dt = np.median(np.diff(t))

        # write the x and y components of the wind into a text file based on read order
        with open(os.path.abspath(outFile), 'w') as f:
            template = \
                (
                        "INPGRID WIND REGULAR {x0:.6f} {y0:.6f} 0.00 {mx:d} {my:d} {dx:.6f} {dy:.6f} &\n" +
                        "NONSTATIONARY {ts} {dt:.2f} SEC {tf} \n" +
                        "READINP WIND +1 '{file}' 3 4 1 0 FORMAT '(({ni:}f7.2))'\n\n"
                )

            f.write(template.format(x0=x0, y0=y0, mx=(nx - 1), my=(ny - 1), dx=dx, dy=dy,
                                    ts=ts, dt=dt, tf=tf, file=outFile, ni=nx))

            for ii in it:
                # write the time stamp header
                f.write(datestr(t[ii], '%Y%m%d.%H%M%S') + '\n')

                # write the x-component
                np.savetxt(f, np.tile(u[ii], (2, 2)), '%7.2f', delimiter='')

                # write the y-component
                np.savetxt(f, np.tile(v[ii], (2, 2)), '%7.2f', delimiter='')


## 

def writeWindGridBc(inFile, outFile, swanModel):
    """
    Converts ERA5 wind data in netCDF4 format to SWAN input grid text file

    note: currently does not work for ERA5 grids with only 1 point. Must have 2 or more points in each dimension.
    """

    # get input netCDF4 file handle
    nc = Dataset(inFile, 'r')

    # get longitude, latitude and time arrays
    x = nc['longitude'][:].data.astype(float)
    y = nc['latitude'][:].data.astype(float)
    # flip the y-axis (otherwise descending)
    # Check if latitudes are descending (going from north to south)
    if np.all(np.diff(y) < 0):
        y = np.flip(y)

    nc_format = float(nc.getncattr('Conventions').split('-')[1])

    if 'valid_time' in nc.variables:
        time_var = 'valid_time'
        print('Variable valid_time found')
    elif 'time' in nc.variables:
        time_var = 'time'
        print('Variable time found')
    else:
        raise KeyError("No se encontró ni 'valid_time' ni 'time' en el archivo NetCDF.")
    
    units = nc[time_var].getncattr("units")
    if "seconds since 1970" in units:
        t = nc[time_var][:].data.astype(float) + datenum((1970, 1, 1))
    elif "hours since 1900" in units:
        t = nc[time_var][:].data.astype(float) * 3600 + datenum((1900, 1, 1))
    else:
        raise ValueError(f"Unsupported time units: {units}")

    # check if points need to be transformed
    srsA = swanModel.swanGrid.getSrs()

    srsB = osr.SpatialReference()

    result = srsB.ImportFromEPSG(4326)
    if result != 0:
        raise ValueError("No se pudo importar EPSG:4326 en srsB")

    transform = False
    if srsA is not None:
        if not isSrsEquivalent(srsA, srsB):
            transform = True

    # get the model spatial boundary
    mb = swanModel.swanGrid.getBoundary()

    # transform boundary to WGS84 if needed
    if transform:
        
        print('Aca estamos printeando el srsA')
        print(srsA)

        print('Aca estamos printeando el srsB')
        print(srsB)


        xt, yt = transformPoints(mb[:, 0], mb[:, 1], srsA, srsB)

        mb = np.column_stack((xt, yt))

    # get the model x, y, and t limits
    x1, x2 = mb[:, 0].min(), mb[:, 0].max()
    y1, y2 = mb[:, 1].min(), mb[:, 1].max()
    t1, t2 = swanModel.timeStart, swanModel.timeEnd

    # find interpolation indices for model domain
    ix = getInterpolationIndex(x, x1, x2)
    iy = getInterpolationIndex(y, y1, y2)

    it = getInterpolationIndex(t, t1, t2)

    # if nothing in model domain, return none
    if any([len(ix) == 0, len(iy) == 0, len(it) < 2]):
        return None

    # subset dimensions
    ts, xs, ys = t[it], x[ix], y[iy]

    # VALIDACIÓN para evitar errores de 0 puntos
    assert xs.size >= 2, "Not enough x-points for SWAN input."
    assert ys.size >= 2, "Not enough y-points for SWAN input."

    # transform points to model SRS, note: once transformed points will no longer be a grid
    # therefore average used to approximate grid. This will not exactly match original points.
    if transform:
        xs, ys = transformPoints(xs, ys, srsB, srsA, grid=True)
        xs, ys = np.mean(xs, axis=0), np.mean(ys, axis=1)

    # write a check file of points
    folder = os.path.split(outFile)[0]
    name, _ = os.path.splitext(os.path.split(outFile)[1])
    nametmp =  name + '.csv'
    checkFile = Path(folder) / nametmp

    with open(checkFile, 'w') as f:
        f.write('X,Y\n')

        xp, yp = np.meshgrid(xs, ys)
        xp, yp = xp.flatten(), yp.flatten()

        np.savetxt(f, np.column_stack((xp, yp)), delimiter=',')

    # get SWAN INPGRID and READINP parameters
    x0, nx, dx = xs[0], xs.size, np.median(np.diff(xs))
    y0, ny, dy = ys[0], ys.size, np.median(np.diff(ys))
    ts, tf, dt = datestr(ts[0], '%Y%m%d.%H%M%S'), \
                 datestr(ts[-1], '%Y%m%d.%H%M%S'), \
                 np.median(np.diff(ts))

    # write the x and y components of the wind into a text file based on read order
    with open(os.path.abspath(outFile), 'w') as f:
        template = \
            (
                    "INPGRID WIND REGULAR {x0:.6f} {y0:.6f} 0.00 {mx:d} {my:d} {dx:.6f} {dy:.6f} &\n" +
                    "NONSTATIONARY {ts} {dt:.2f} SEC {tf} \n" +
                    "READINP WIND +1 '{file}' 3 4 1 0 FORMAT '(({ni:}f7.2))'\n\n"
            )

        f.write(template.format(x0=x0, y0=y0, mx=(nx - 1), my=(ny - 1), dx=dx, dy=dy,
                                ts=ts, dt=dt, tf=tf, file=outFile, ni=nx))


        for ii in it:
            # intentar leer u10 o eastward_wind
            if 'u10' in nc.variables:
                u = np.flipud(nc['u10'][ii, :, :].data)
            elif 'eastward_wind' in nc.variables:
                u = np.flipud(nc['eastward_wind'][ii, :, :].data)
            else:
                raise KeyError("No se encontró ni 'u10' ni 'eastward_wind' en el archivo NetCDF.")

            # intentar leer v10 o northward_wind
            if 'v10' in nc.variables:
                v = np.flipud(nc['v10'][ii, :, :].data)
            elif 'northward_wind' in nc.variables:
                v = np.flipud(nc['northward_wind'][ii, :, :].data)
            else:
                raise KeyError("No se encontró ni 'v10' ni 'northward_wind' en el archivo NetCDF.")

            # truncate to grid size
            u = u[np.ix_(iy, ix)]
            v = v[np.ix_(iy, ix)]

            # write the time stamp header
            f.write(datestr(t[ii], '%Y%m%d.%H%M%S') + '\n')

            # write the x-component
            np.savetxt(f, u, '%7.2f', delimiter='')

            # write the y-component
            np.savetxt(f, v, '%7.2f', delimiter='')

    # remove netCDF4 file handle
    nc.close()

def writeBottomGridBc(inFiles, outFile, swanModel):
    """
    Inspect grid point elevations from many terrain models and write the output to a SWAN input grid text file
    """

    # get convenience handles
    x = swanModel.swanGrid.getX()
    y = swanModel.swanGrid.getY()
    dx = swanModel.swanGrid.getDx()
    dy = swanModel.swanGrid.getDy()
    nx = swanModel.swanGrid.getNx()
    ny = swanModel.swanGrid.getNy()
    r = swanModel.swanGrid.getRotation()
    srs = swanModel.swanGrid.getSrs()
    faces = swanModel.swanGrid.getFaces()

    # iteratively inspect the raster at points (x, y)
    zp0 = np.zeros((ny*nx,), dtype=np.float64)*np.nan
    
    for inFile in inFiles:
        if isMesh(inFile):
            zp1 = meshInspect(inFile, x.flatten(), y.flatten(), srs)
        else:
            import rasterio
            with rasterio.open(inFile) as src:
                bounds = src.bounds
                # Conversión robusta de longitudes
                if bounds.left > 0 and np.any(x.flatten() < 0):
                    x_to_use = x.flatten().copy()
                    x_to_use[x_to_use < 0] += 360
                elif bounds.left < 0 and np.any(x.flatten() > 180):
                    x_to_use = x.flatten().copy()
                    x_to_use[x_to_use > 180] -= 360
                else:
                    x_to_use = x.flatten()

                # Ahora consulta usando x_to_use
                vals = []
                for xi, yi in zip(x_to_use, y.flatten()):
                    for val in src.sample([(xi, yi)]):
                        vals.append(val[0])
                zp1 = np.array(vals)
                # Reemplaza fill values por np.nan
                FILL_VALUES = [-9969209968386869046778552952102584320.0, 9.96921e+36]
                for fv in FILL_VALUES:
                    zp1 = np.where(zp1 == fv, np.nan, zp1)
        # overwrite valid values
        zp0[np.isnan(zp1) == 0] = zp1[np.isnan(zp1) == 0]

    # write bathymetry grid to check file for visualization in GIS
    folder = os.path.split(outFile)[0]
    name, _ = os.path.splitext(os.path.split(outFile)[1])
    nametmp =  name + '.nc'
    checkFile = Path(folder) / nametmp

    meshWrite(checkFile, x.flatten(), y.flatten(), zp0, faces, swanModel.swanGrid.getSrs())

    # reshape to ny by nx array for SWAN format
    zp0 = np.reshape(zp0, (ny, nx))

    # open SWAN bottom grid output file
    with open(os.path.abspath(outFile), 'w') as f:


        # put the control file command in header (be sure to specify nhedf=3)
        template = \
            (
                    "INPGRID BOTTOM REGULAR {x0:.6f} {y0:.6f} {r:.2f} {mx:d} {my:d} {dx:.6f} {dy:.6f} &\n" +
                    "READINP BOTTOM +1 '{file}' 3 3 FORMAT '(({nj}f8.2))'\n\n"
            )



        f.write(template.format(x0=x[0, 0], y0=y[0, 0], r=r, mx=(nx - 1), my=(ny - 1), dx=dx, dy=dy, file=outFile, nj=nx))

        # write out grid point elevations with read order 3
        np.savetxt(f, -1*zp0, '%8.2f', delimiter='')
        print("Bathymetry file written")


def writeWaveGridBc(inFile, outFile, swanModel, spread=12):
    """ERA5 parametric wave boundary condition using TPAR files"""

    # get input wave data file netCDF4 handle
    nc = Dataset(inFile, 'r')

    # get longitude, latitude and time arrays
    x = nc['longitude'][:].data.astype(float)
    y = nc['latitude'][:].data.astype(float)
    # flip the y-axis (otherwise descending)
    # y = np.flip(y)

    # Check if latitudes are descending (going from north to south)
    if ~np.all(np.diff(y) < 0):
        y = np.flip(y)



    nc_format = float(nc.getncattr('Conventions').split('-')[1])

    if 'valid_time' in nc.variables:
        time_var = 'valid_time'
        print('Variable valid_time found')
    elif 'time' in nc.variables:
        time_var = 'time'
        print('Variable time found')
    else:
        raise KeyError("No se encontró ni 'valid_time' ni 'time' en el archivo NetCDF.")
    
    units = nc[time_var].getncattr("units")
    if "seconds since 1970" in units:
        t = nc[time_var][:].data.astype(float) + datenum((1970, 1, 1))
    elif "hours since 1900" in units:
        t = nc[time_var][:].data.astype(float) * 3600 + datenum((1900, 1, 1))
    else:
        raise ValueError(f"Unsupported time units: {units}")

    # check if points need to be transformed
    srsA = swanModel.swanGrid.getSrs()
    srsB = osr.SpatialReference()
    srsB.ImportFromEPSG(4326)

    transform = False
    if srsA is not None:
        if not isSrsEquivalent(srsA, srsB):
            transform = True

    # get the model spatial boundary
    mb = swanModel.swanGrid.getBoundary()

    # transform boundary to WGS84 if needed
    if transform:
        mb = np.column_stack(transformPoints(mb[:, 0], mb[:, 1], srsA, srsB))

    # get average resolution of ERA5 grid
    ds = np.mean(np.abs(np.concatenate((np.diff(x), np.diff(y)))))

    # make empty arrays for interpolation points
    xp, yp = np.array((), dtype=float), np.array((), dtype=float)

    # iterate over each side\stretch
    for aa in range(4):
        # get start & end coordinates
        x1, y1 = tuple(mb[aa])
        x2, y2 = tuple(mb[aa + 1])

        # get number of points to add
        dist = np.hypot(x2 - x1, y2 - y1)
        add = np.round(dist/ds) - 1

        # get interpolated points
        xp_aa = np.interp(np.arange(add + 1), [0, add + 1], [x1, x2])
        yp_aa = np.interp(np.arange(add + 1), [0, add + 1], [y1, y2])

        # stack the arrays
        xp = np.hstack((xp, xp_aa))
        yp = np.hstack((yp, yp_aa))

    # get time limits from model
    t1, t2 = swanModel.timeStart, swanModel.timeEnd

    # find interpolation indices for model domain
    it = getInterpolationIndex(t, t1, t2)

    # if nothing in model domain, return none
    if len(it) < 2:
        return None

    # pre-allocate data for points
    hs = np.zeros((it.size, xp.size))
    mwd = np.zeros((it.size, xp.size))
    pwp = np.zeros((it.size, xp.size))

    for aa in range(len(it)):
        # altura significativa
        if 'swh' in nc.variables:
            hs_ii = nc['swh'][it[aa], :, :].data.astype('float')
        elif 'VHM0' in nc.variables:
            hs_ii = nc['VHM0'][it[aa], :, :].data.astype('float')
        else:
            raise KeyError("No se encontró ni 'swh' ni 'VHM0' en el NetCDF.")

        # dirección media
        if 'mwd' in nc.variables:
            mwd_ii = nc['mwd'][it[aa], :, :].data.astype('float')
        elif 'VMDR' in nc.variables:
            mwd_ii = nc['VMDR'][it[aa], :, :].data.astype('float')
        else:
            raise KeyError("No se encontró ni 'mwd' ni 'VMDR' en el NetCDF.")

        # periodo pico
        if 'pp1d' in nc.variables:
            pwp_ii = nc['pp1d'][it[aa], :, :].data.astype('float')
        elif 'VTPK' in nc.variables:
            pwp_ii = nc['VTPK'][it[aa], :, :].data.astype('float')
        else:
            raise KeyError("No se encontró ni 'pp1d' ni 'VTPK' en el NetCDF.")

        #if nc_format <= 1.6:
        if "hours since 1900" in units:    
            # Fill with zero values
            bad = (hs_ii == -32767) | (mwd_ii == -32767) | (pwp_ii == -32767)
            hs_ii[bad] = np.nan
            mwd_ii[bad] = np.nan
            pwp_ii[bad] = np.nan

            # convert direction from nautical to cartesian
            mwd_ii = convertDirection(mwd_ii)

            # convert direction to unit vector
            xc_ii = np.cos(mwd_ii * np.pi / 180)
            yc_ii = np.sin(mwd_ii * np.pi / 180)

            # interpolate vector components to points
            xc_pp = linearGrid(x, y, xc_ii, xp, yp)
            yc_pp = linearGrid(x, y, yc_ii, xp, yp)
            pwp_pp = linearGrid(x, y, pwp_ii, xp, yp)
            hs_pp = linearGrid(x, y, hs_ii, xp, yp)

            # convert unit vector back into direction
            mwd_pp = np.arctan2(yc_pp, xc_pp) * 180 / np.pi
            mwd_pp[mwd_pp < 0] += 360

        else:
            # archivos CF-1.7 o superior: se asume que los datos están limpios

            hs_pp = linearGrid(x, y, hs_ii, xp, yp)
            mwd_pp = linearGrid(x, y, mwd_ii, xp, yp)
            pwp_pp = linearGrid(x, y, pwp_ii, xp, yp)

        # store the data (en ambos casos)
        hs[aa, :] = hs_pp
        mwd[aa, :] = mwd_pp
        pwp[aa, :] = pwp_pp

    # convert points to model SRS
    if transform:
        xp, yp = transformPoints(xp, yp, srsB, srsA)

    # write check file for points
    folder = os.path.split(outFile)[0]
    name, _ = os.path.splitext(os.path.split(outFile)[1])
    nametmp =  name + '.csv'
    checkFile = Path(folder) / nametmp

    with open(checkFile, 'w') as f:
        f.write('X,Y,NAME\n')
        for aa in range(len(xp)):
            f.write('{:.6f},{:.6f},P{:02d}\n'.format(xp[aa], yp[aa], aa))

    # start the control string
    #    nx, ny = swanModel.swanGrid.getNx(), swanModel.swanGrid.getNy()
    
    # TODO this is hardcoded, it should be removed and revised
    nx = swanModel.swanGrid.getNx() - 2
    ny = swanModel.swanGrid.getNy() - 2

    # 
    controlString = \
        "BOUND SHAPESPEC PM PEAK DSPR DEGREES\n" + \
        "BOUNDSPEC SEGMENT IJ 0,0 {mx},0 {mx},{my} 0,{my} 0,0 VARIABLE FILE &\n".format(mx=nx, my=ny)

    # start distance counter
    distance = 0.0

    # iterate over interpolated points
    for aa in range(len(xp)):
        # get TPAR output file path
        nametmp =  name + '_P{:02d}.tpar'.format(aa)
        tparFile = Path(folder) / nametmp

        # open TPAR file
        with open(os.path.abspath(tparFile), 'w') as f:
            # write the header
            f.write("TPAR\n")
            # format time in string array
            tString = datestr(t[it], '%Y%m%d.%H%M%S').astype(object)

            # write data to tpar
            data = np.column_stack((tString, hs[:, aa], pwp[:, aa], mwd[:, aa], spread*np.ones(it.shape)))

            np.savetxt(f, data, '%s,%.2f,%.2f,%.2f,%.2f')

        # get distance from start point
        if aa > 0:
            dx = xp[aa] - xp[aa - 1]
            dy = yp[aa] - yp[aa - 1]
            ds = np.hypot(dx, dy)
        else:
            ds = 0
        distance += ds

        # add statement to control string
        if np.all(~np.isnan(hs[:, aa])):
            controlString += "{:.6f} '{}' 1 &\n".format(distance, tparFile)

    # write the control file
    nametmp = name + '.txt'
    controlFile = Path(folder) / nametmp

    with open(os.path.abspath(controlFile), 'w') as f:
        f.write(controlString)

class SwanGrid:
    """
    A geometry class for building grids from a polygon. This class is unaware of the SWAN model and
    is supposed to be used only for generating the spatial component of the computational grid.

    PARAMETERS
    ----------
    name : string
        A unique name for identifying the grid
    polygon : array
        Polygon of model domain as (n, 2) array of (xp, yp)
    rotation : float
        Grid rotation, anticlockwise from positive x-axis
    dx : float
        x-cell size along x-axis of rotated grid
    dy : float
        y-cell size along y-axis of rotated grid
    nx : integer
        Number of x points in grid
    ny : integer
        Number of y points in grid
    nc : integer
        Number of total grid cells
    srs : osr.SpatialRef
        SRS object that defines spatial reference system

    ATTRIBUTES
    ----------
    w : float
        Width of grid along local x-axis
    h : float
        Height of grid along local y-axis
    x : array
        (ny, nx) array of grid point x-coordinates
    y : array
        (ny, nx) array of grid point y-coordinates
    """

    # constructor function
    def __init__(self, name, polygon=None, rotation=None, dx=None,
                 dy=None, nx=None, ny=None, nv=None, srs=None):

        # protected attributes
        self._name = name  # name of the grid
        self._polygon = polygon  # bounding polygon
        self._rotation = rotation  # rotation in degrees
        self._dx = dx  # cell size along x-axis
        self._dy = dy  # cell size along y-axis
        self._nx = nx  # number of points along x-axis
        self._ny = ny  # number of points along y-axis
        self._nv = nv  # total number of vertices (nx*ny)
        self._srs = srs  # osr.SpatialReference object

        self._w = None  # grid width
        self._h = None  # grid height
        self._x = None  # grid x coordinate
        self._y = None  # grid y coordinate

        self._isDefined = False  # flag to check if grid is defined

        # build the grid
        self._buildGrid()

    # protected member functions
    def _buildGrid(self):
        """Calculates grid properties"""

        # assume grid no longer defined
        self._isDefined = False

        # check polygon is defined
        if self._polygon is None:
            return

        # Get centroid for translation
        centre = np.mean(self._polygon, axis=0)

        # Get rotation in radians
        if self._rotation is not None:
            alpha = self._rotation * np.pi / 180
        else:
            alpha = 0

        # Get rotation matrix
        cos_a, sin_a = np.cos(alpha), np.sin(alpha)
        rm = np.array([[cos_a, sin_a], [-sin_a, cos_a]])

        # Rotate points into practical reference frame
        bounds_r = np.matmul(rm, (self._polygon - centre).transpose()).transpose()

        # Get rotated limits
        xlr = [np.min(bounds_r[:, 0]), np.max(bounds_r[:, 0])]
        ylr = [np.min(bounds_r[:, 1]), np.max(bounds_r[:, 1])]

        self._w = np.abs(np.diff(xlr))[0]
        self._h = np.abs(np.diff(ylr))[0]

        # solve for nx and ny, check all cases
        if self._nv is not None:  # CASE 1: nv defined
            # check for valid input
            if self._nv <= 0: return
            # calculate nx and ny from nv
            ds = np.sqrt(self._w * self._h / self._nv)
            self._nx = int(np.round(self._w / ds)) + 1
            self._ny = int(np.round(self._h / ds)) + 1
        elif self._dx is not None and self._dy is not None:  # CASE 2: dx and dy defined
            # check for valid input
            if (self._dx <= 0) or (self._dy <= 0): return
            # calculate nx and ny from dx and dy
            self._nx = int(np.round(self._w / self._dx)) + 1
            self._ny = int(np.round(self._h / self._dy)) + 1
        elif self._dx is not None and self._ny is not None:  # CASE 3: dx and ny defined
            # check for valid input
            if self._dx <= 0 or self._ny <= 1: return
            # calculate nx from dx
            self._nx = int(np.round(self._w / self._dx)) + 1
        elif self._nx is not None and self._dy is not None:  # CASE 4: nx and dy defined
            # check for valid input
            if self._dy <= 0 or self._nx <= 1: return
            # calculate ny from dy
            self._ny = int(np.round(self._h / self._dy)) + 1
        elif self._nx is not None and self._ny is not None:  # CASE 5: nx and ny defined
            # check for valid input
            if self._nx <= 1 or self._ny <= 1: return
        else:  # insufficient parameters to define grid
            return

        # recalculate other parameters
        self._dx = self._w / (self._nx - 1)
        self._dy = self._h / (self._ny - 1)
        self._nv = self._nx * self._ny

        # Get 1D rotated x & y
        xr = np.linspace(xlr[0], xlr[1], self._nx)
        yr = np.linspace(ylr[0], ylr[1], self._ny)

        # Mesh the points
        xr, yr = np.meshgrid(xr, yr)
        xr = xr.flatten()
        yr = yr.flatten()

        # Transform points back to original reference frame
        rm = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
        xy = np.matmul(rm, np.vstack((xr, yr))).transpose() + centre

        # Store gridded coordinates
        self._x = np.reshape(xy[:, 0], (self._ny, self._nx))
        self._y = np.reshape(xy[:, 1], (self._ny, self._nx))

        # build complete, everything is ok
        self._isDefined = True

    # public member functions
    def getName(self):
        return self._name

    def setName(self, name):
        self._name = name

    def getPolygon(self):
        return self._polygon

    def setPolygon(self, polygon):
        # reset variables
        self._polygon = polygon
        self._nx = None
        self._ny = None
        self._nv = None
        self._w = None
        self._h = None

        # rebuild grid
        self._buildGrid()

    def getRotation(self):
        if self._rotation is None:
            return 0
        else:
            return self._rotation

    def setRotation(self, rotation):
        # reset variables
        self._rotation = rotation
        self._nx = None
        self._ny = None
        self._nv = None
        self._w = None
        self._h = None

        # rebuild grid
        self._buildGrid()

    def getDx(self):
        return self._dx

    def setDx(self, dx):
        # reset variables
        self._dx = dx
        self._nx = None
        self._nv = None

        # rebuild grid
        self._buildGrid()

    def getDy(self):
        return self._dy

    def setDy(self, dy):
        # reset variables
        self._dy = dy
        self._ny = None
        self._nv = None

        # rebuild grid
        self._buildGrid()

    def getNx(self):
        return self._nx

    def setNx(self, nx):
        # reset variables
        self._nx = nx
        self._dx = None
        self._nv = None

        # rebuild grid
        self._buildGrid()

    def getNy(self):
        return self._ny

    def setNy(self, ny):
        # reset variables
        self._ny = ny
        self._dy = None
        self._nv = None

        # rebuild grid
        self._buildGrid()

    def getNv(self):
        return self._nv

    def setNv(self, nv):
        # reset variables
        self._nv = nv
        self._nx = None
        self._ny = None
        self._dx = None
        self._dy = None

        # rebuild grid
        self._buildGrid()

    def getSrs(self):
        return self._srs

    def setSrs(self, srs):
        self._srs = srs

    def getW(self):
        return self._w

    def getH(self):
        return self._h

    def getX(self):
        return self._x

    def getY(self):
        return self._y

    def isDefined(self):
        return self._isDefined

    def isSpherical(self):
        return isSpherical(self._srs)

    def getFaces(self):
        ni, nj = self._ny, self._nx

        ii = np.arange(ni - 1)
        jj = np.arange(nj - 1)

        jj, ii = np.meshgrid(jj, ii)
        bl = np.ravel(jj + nj * ii)
        tl = np.ravel(jj + nj * (ii + 1))

        faces = np.vstack((bl, bl + 1, tl + 1, tl)).transpose()

        return faces

    def getNodes(self):
        nodeX = self._x.flatten()
        nodeY = self._y.flatten()

        nodes = np.column_stack((nodeX, nodeY))

        return nodes

    def getBoundary(self):
        index = [[0, 0], [0, -1], [-1, -1], [-1, 0], [0, 0]]

        x = [self._x[ii, jj] for (ii, jj) in index]
        y = [self._y[ii, jj] for (ii, jj) in index]

        return np.column_stack((x, y))

    def getBoundingBox(self):
        boundary = self.getBoundary()
        x1 = boundary[:, 0].min()
        y1 = boundary[:, 1].min()
        x2 = boundary[:, 0].max()
        y2 = boundary[:, 1].max()

        return x1, y1, x2, y2


class SwanConfig:
    """
    A class for configuring SWAN models. This class prevents overloading the SwanModel class and also allows for
    a single configuration across many runs. This class is responsible for holding information for the SWAN
    control file SET, MODE, CGRID (spectra domain only) and NUMERIC commands. It is designed to hold all parameters
    which are likely to be constant across many runs.

    PARAMETERS
    ----------
    level : float
        Constant water level in SET command
    dirNorth : float
        Direction of North in SET command
    depthMin : float
        Minimum depth in SET command
    maxMsessages : int
        Maximum number of messages in SET command
    maxError : int
        Maximum error level
    """

    def __init__(self, level=0, dirNorth=90, depthMin=0.05, maxMessages=200, maxError=1, gravity=9.81, rho=1025,
                 maxDrag=99999., stationary=False, spherical=True, numDir=36, freqMin=0.04, freqMax=1.00, numFreq=34,
                 dAbsolute=0.005, dRelative=0.01, curvature=0.005, numPoints=99.5, maxIterations=5, limiter=0.01,
                 outputVars=None, outputPoints=None, outputType='.mat', timeStep=3600, timeUnit='SEC'):

        # store all parameter names in order, useful for iterating
        self.parameters = inspect.getfullargspec(self.__init__).args[1:]

        # public attributes for SET command
        self.level = level
        self.dirNorth = dirNorth
        self.depthMin = depthMin
        self.maxMessages = maxMessages
        self.maxError = maxError
        self.gravity = gravity
        self.rho = rho
        self.maxDrag = maxDrag

        # public attributes for the MODE command
        self.stationary = stationary
        self.spherical = spherical

        # public attributes for the CGRID command
        self.numDir = numDir
        self.freqMin = freqMin
        self.freqMax = freqMax
        self.numFreq = numFreq

        # public attributes for the NUMERIC command
        self.dAbsolute = dAbsolute
        self.dRelative = dRelative
        self.curvature = curvature
        self.numPoints = numPoints
        self.maxIterations = maxIterations
        self.limiter = limiter

        # public attributes for the output command
        self.outputVars = outputVars
        self.outputPoints = outputPoints
        self.outputType = outputType

        if self.outputVars is None:
            self.outputVars = 'XP YP HSIGN TPS PDIR DIR UBOT TMBOT FORCE DEPTH'

        # public attributes for temporal commands
        self.timeStep = timeStep
        self.timeUnit = timeUnit

    @staticmethod
    def isFloat(string):
        try:
            float(string)
            return True
        except ValueError:
            return False

    @staticmethod
    def isInteger(string):
        try:
            int(string)
            return True
        except ValueError:
            return False

    @staticmethod
    def convertType(string):
        if SwanConfig.isInteger(string):
            return int(string)
        elif SwanConfig.isFloat(string):
            return float(string)
        elif string == 'None':
            return None
        elif string == 'False':
            return False
        elif string == 'True':
            return True
        else:
            return string

    def read(self, file):
        """Simple wrapper for ConfigParser.read() method"""

        # instantiate configparser object
        cp = configparser.ConfigParser()
        cp.optionxform = str

        # read configuration from file
        cp.read(file)

        # pass data from dictionary to object
        for name in self.parameters:
            if name in cp['SWAN CONFIG']:
                setattr(self, name, SwanConfig.convertType(cp['SWAN CONFIG'][name]))

        return self

    def write(self, file):
        """Simple wrapper for ConfigParser.write() method"""

        # instantiate configparser object
        cp = configparser.ConfigParser()
        cp.optionxform = str

        # pass data to the parser dictionary
        cp['SWAN CONFIG'] = {name: str(getattr(self, name)) for name in self.parameters}

        # save it to specified destination
        with open(file, 'w') as f:
            cp.write(f)

class SwanModel:
    """
    A class for a single SWAN run. This class is a simple wrapper, it makes no assumptions about the project
    setup and related file paths. It is made aware of nesting by manually setting the parent attribute. It strictly
    mimics a SWAN control file.

    PARAMETERS
    ----------
    modelName : string
        A unique name for identifying the model
    swanGrid : SwanGrid
        Grid object that defines spatial component of computational domain
    swanConfig : SwanConfig
        Configuration object for model, sets parameters across many models
    modelParent : SwanModel
        Parent model if model is nested (used to configure nest BCs)
    timeStart : float
        Model start time as python time stamp
    timeEnd : float
        Model end time as python time stamp
    controlFile : string
        File path of model control file (writes to this location)
    bottomFile : string
        File path of SWAN bottom grid input text file
    windFile : string
        File path of SWAN wind grid input text file
    specFile : string
        File path of SWAN spectral input text file
    outputFile : string
        File path of SWAN netCDF or .mat output files
    templateFile : string
        File path of template, default: swangis\template.txt
    """

    def __init__(self, modelName, swanGrid, swanConfig=None, modelParent=None, timeStart=None, timeEnd=None,
                 controlFile=None, bottomFile=None, windFile=None, specFile=None, outputFile=None, templateFile=None):

        # protected attributes
        self._parent = None
        self._children = []

        # set protected attributes
        self.setParent(modelParent)

        # set public attributes
        self.modelName = modelName
        self.swanGrid = swanGrid
        self.swanConfig = swanConfig

        self.timeStart = timeStart
        self.timeEnd = timeEnd

        self.controlFile = controlFile
        self.bottomFile = bottomFile
        self.windFile = windFile
        self.specFile = specFile
        self.outputFile = outputFile
        self.templateFile = templateFile

        # set some defaults
        if templateFile is None:  # set default
            self.templateFile = Path(os.path.dirname(__file__)) / 'template'
        if swanConfig is None:
            self.swanConfig = SwanConfig()

    def __call__(self, *args, **kwargs):
        return self.getControlString()

    def getParent(self):
        return self._parent

    def setParent(self, parent):
        # clear connection to current parent
        if self._parent is not None:
            self._parent._children.remove(self)

        # set new parent
        self._parent = parent
        if parent is not None:
            parent._children.append(self)

    def isStationary(self):
        return (self.timeStart is None) or (self.timeEnd is None)

    def getStartUp(self):
        """Returns control file start-up commands"""

        params = \
            {
                'level': self.swanConfig.level,
                'north': self.swanConfig.dirNorth,
                'depmin': self.swanConfig.depthMin,
                'maxmes': self.swanConfig.maxMessages,
                'maxerr': self.swanConfig.maxError,
                'gravity': self.swanConfig.gravity,
                'rho': self.swanConfig.rho,
                'cdcap': self.swanConfig.maxDrag
            }

        # specify the mode (stationary or non-stationary)
        if self.swanConfig.stationary:
            params['mode'] = 'STATIONARY'
        else:
            params['mode'] = 'NONSTATIONARY'

        # specify the CRS type
        if self.swanConfig.spherical:
            params['crs'] = 'SPHERICAL'
        else:
            params['crs'] = 'CARTESIAN'

        # overwrite the CRS type
        if self.swanGrid.isSpherical():
            params['crs'] = 'SPHERICAL'
        else:
            params['crs'] = 'CARTESIAN'

        template = \
            "PROJECT '' ''\n" \
            "SET {level:.2f} {north:.2f} {depmin:.2f} {maxmes:d} {maxerr:d} {gravity:.2f} {rho:.2f} {cdcap:.3f}\n" \
            "MODE {mode:} TWODIMENSIONAL\n" \
            "COORDINATES {crs:}\n"

        startUp = template.format(**params)

        return startUp

    def getCGridInput(self):
        """Returns control file command for computational grid specification."""

        params = \
            {
                'x0': self.swanGrid.getX()[0, 0],
                'y0': self.swanGrid.getY()[0, 0],
                'alpha': self.swanGrid.getRotation(),
                'w': self.swanGrid.getW(),
                'h': self.swanGrid.getH(),
                'mx': self.swanGrid.getNx() - 1,
                'my': self.swanGrid.getNy() - 1,
                'nd': self.swanConfig.numDir,
                'fmin': self.swanConfig.freqMin,
                'fmax': self.swanConfig.freqMax,
                'nf': self.swanConfig.numFreq,
            }

        template = \
            "CGRID REGULAR {x0:.6f} {y0:.6f} {alpha:.2f} {w:.6f} {h:.6f} {mx:d} {my:d} &\n" \
            "CIRCLE {nd:d} {fmin:.3f} {fmax:.3f} {nf:d}\n"

        if self.swanGrid.isDefined():
            cGridInput = template.format(**params)
        else:
            cGridInput = '$ N/A NO CGRID INPUT'

        return cGridInput

    def getBottomInput(self):
        """Returns control file input command for bathymetry."""

        if self.bottomFile is not None:
            with open(os.path.abspath(self.bottomFile), 'r') as f:
                bottomInput = f.readline() + f.readline()
        else:
            bottomInput = '$ N/A NO BOTTOM INPUT'

        return bottomInput

    def getWindInput(self):
        """Returns control file input command for wind"""

        if self.windFile is not None:
            with open(os.path.abspath(self.windFile), 'r') as f:
                windInput = f.readline() + f.readline() + f.readline()
        else:
            windInput = '$ N/A NO WIND INPUT'

        return windInput

    def getSurfaceInput(self):
        return '$ N/A NO SURFACE INPUT'

    def getCurrentInput(self):
        return '$ N/A NO CURRENT INPUT'

    def getSpecInput(self):
        """Returns control file input command for spectral boundary condition"""

        if self.specFile is not None:
            with open(os.path.abspath(self.specFile), 'r') as f:
                specInput = f.read()
        else:
            specInput = '$ N/A NO SPECTRAL INPUT'

        return specInput

    def getNestInput(self):
        """Returns control file input command for parent model boundary condition if model has parent."""

        if self.getParent() is not None:
            # find nest file in output folder of parent
            folder = os.path.split(self.getParent().outputFile)[0]
            name, _ = os.path.splitext(os.path.split(self.controlFile)[1])
            nametmp = name + '.nest'

            nestFile = Path(folder) / nametmp

            nestInput = "BOUNDNEST1 NEST '{}' CLOSED".format(nestFile)
        else:
            nestInput = '$ N/A NO NEST INPUT'

        return nestInput

    def getBoundaryInput(self):
        """Returns control file input command for boundary input."""

        if self.getParent() is None:
            boundaryInput = self.getSpecInput()
        else:
            boundaryInput = self.getNestInput()

        return boundaryInput

    def getNestOutput(self):
        """Returns control file output command for nested child model boundary conditions"""

        template = \
            "NGRID '{name}' {x0:.6f} {y0:6f} {alpha:.2f} {w:.6f} {h:.6f}\n" \
            "NESTOUT '{name}' '{file}' OUTPUT {timeStart} {timeStep} {timeUnit}\n"

        nestOutput = ''
        for aa, child in enumerate(self._children):
            folder = os.path.split(self.outputFile)[0]
            name, _ = os.path.splitext(os.path.split(child.controlFile)[1])
            nametmp = name +'.nest'

            nestFile = Path(folder) / nametmp

            grid = child.swanGrid

            nestOutput += template.format(x0=grid.getX()[0, 0], y0=grid.getY()[0, 0], alpha=grid.getRotation(),
                                          w=grid.getW(), h=grid.getH(), name='NEST{:02d}'.format(aa), file=nestFile,
                                          timeStart=datestr(child.timeStart, '%Y%m%d.%H%M%S'),
                                          timeStep=child.swanConfig.timeStep, timeUnit=child.swanConfig.timeUnit)

        return nestOutput

    def getGridOutput(self):
        """Returns control file output command for netCDF4 block output"""

        template = \
            "BLOCK 'COMPGRID' NOHEADER '{outFile}' &\n" \
            "LAY-OUT 3 {variables} OUTPUT {timeStart} {timeStep} {timeUnit}\n"

        folder = os.path.split(self.outputFile)[0]
        name, _ = os.path.splitext(os.path.split(self.outputFile)[1])
        nametmp = name + self.swanConfig.outputType

        outFile = Path(folder) / nametmp


        gridOutput = template.format(outFile=outFile, variables=self.swanConfig.outputVars,
                                   timeStart=datestr(self.timeStart, '%Y%m%d.%H%M%S'),
                                   timeStep=self.swanConfig.timeStep, timeUnit=self.swanConfig.timeUnit)

        return gridOutput

    def getSpecOutput(self):
        """Returns control file output command for point spectral output"""

        if self.swanConfig.outputPoints is None or self.swanConfig.outputPoints == 'None':
            return '$ N/A No spectral output'

        template = \
            "POINTS '{desc}' {x:.6f} {y:.6f}\n" \
            "SPECOUT  '{desc}' SPEC2D ABS '{file}' &\n" \
            "OUTPUT {timeStart} {timeStep} {timeUnit}\n"

        folder = os.path.split(self.outputFile)[0]
        name, _ = os.path.splitext(os.path.split(self.outputFile)[1])

        specOutput = ''
        for desc, (x, y) in self.swanConfig.outputPoints.items():
            nametmp = name + '_' + desc + '.spec'
            specFile = Path(folder) / nametmp

            specOutput += template.format(desc=desc, x=x, y=y, file=specFile,
                                           timeStart=datestr(self.timeStart, '%Y%m%d.%H%M%S'),
                                           timeStep=self.swanConfig.timeStep, timeUnit=self.swanConfig.timeUnit)

        return specOutput

    def getOutput(self):
        return self.getNestOutput() + '\n' + self.getGridOutput() + '\n' + self.getSpecOutput()

    def getNumeric(self):
        """Returns control file numeric command"""

        params = \
            {
                'dabs': self.swanConfig.dAbsolute,
                'drel': self.swanConfig.dRelative,
                'curvat': self.swanConfig.curvature,
                'npnts': self.swanConfig.numPoints,
                'maxitr': self.swanConfig.maxIterations,
                'limiter': self.swanConfig.limiter,
            }

        # specify the mode (stationary or non-stationary)
        if self.swanConfig.stationary:
            params['mode'] = 'STATIONARY'
        else:
            params['mode'] = 'NONSTATIONARY'

        template = \
            "NUMERIC ACCUR {dabs:.3f} {drel:.3f} {curvat:.3f} {npnts:.2f} {mode:} {maxitr:d} {limiter:.2f}"

        numeric = template.format(**params)

        return numeric

    def getComputeLockUp(self):

        template = \
            "INITIAL DEFAULT\n" \
            "COMPUTE NONSTATIONARY {timeStart} {timeStep} {timeUnit} {timeEnd}\n" \
            "STOP\n"

        return template.format(timeStart=datestr(self.timeStart, '%Y%m%d.%H%M%S'),
                               timeStep=self.swanConfig.timeStep, timeUnit=self.swanConfig.timeUnit,
                               timeEnd=datestr(self.timeEnd, '%Y%m%d.%H%M%S'))

    def getControlString(self):
        """Returns contents of control file as string, perhaps this should be returned on call?"""

        # read control file template
        with open(self.templateFile, 'r') as f:
            template = f.read()

        # create command hash table
        commands = \
            {
                'START_UP': self.getStartUp(),
                'CGRID_INPUT': self.getCGridInput(),
                'BOTTOM_INPUT': self.getBottomInput(),
                'WIND_INPUT': self.getWindInput(),
                'SURFACE_INPUT': self.getSurfaceInput(),
                'CURRENT_INPUT': self.getCurrentInput(),
                'BOUNDARY_INPUT': self.getBoundaryInput(),
                'OUTPUT': self.getOutput(),
                'NUMERIC': self.getNumeric(),
                'COMPUTE_LOCKUP': self.getComputeLockUp()
            }

        # format the template with commands
        controlString = template.format(**commands)

        # return the control file contents
        return controlString

    def writeControlFile(self):
        """Writes the control file using set model parameters"""
        with open(self.controlFile, 'w') as f:
            f.write(self.getControlString())


class SwanBuilder:
    """
    The SwanBuilder is used for automating nesting, run generation and file management. In future it would be good
    for people to make their own builder classes from a template to allow to easy customisation. Current style is
    rigid for QA purposes.

    PARAMETERS
    ----------
    rootFolder : string
        Path of SWAN folder. Default is current directory './'.
    templateSource : string
        Path to template file used to generate control files. Default used if not specified.
    configSource : string
        Path to configuration file used to set model parameters. Default used if not specified.
    gridSource : string
        Path to vector file containing SWAN grid data as polygon features with fields.
    bottomSource : list
        List of file paths to SWAN bathymetry sources in descending priority.
    windSource : string
        File path to wind source (ERA5 .nc file or .csv time series)
    waveSource : string
        File path to wave source (ERA5 .nc file or .csv time series)
    """

    def __init__(self, rootFolder='./', templateSource=None, configSource=None,
                 gridSource=None, bottomSource=None, windSource=None, waveSource=None):


        # protected attributes
        self._rootFolder = None
        self._geoFolder = None
        self._bcFolder = None
        self._simFolder = None
        self._resFolder = None

        # public attributes
        self.templateSource = templateSource
        self.configSource = configSource
        self.gridSource = gridSource
        self.bottomSource = bottomSource
        self.windSource = windSource
        self.waveSource = waveSource

        # set the root folder
        self.setRootFolder(rootFolder)

    def getRootFolder(self):
        return self._rootFolder

    def setRootFolder(self, rootFolder):
        # set root folder attribute
        self._rootFolder = os.path.abspath(rootFolder)

        # use default names for respective project sub-folders
        self._geoFolder = os.path.join(rootFolder, '01_geometry')
        self._bcFolder = os.path.join(rootFolder, '02_bc_dbase')
        self._simFolder = os.path.join(rootFolder, '03_simulation')
        self._resFolder = os.path.join(rootFolder, '04_results')

    def buildRun(self, timeStart, timeEnd, prefix='', suffix=''):
        # make the root folder if it doesn't exist
        if not os.path.exists(self._rootFolder):
            os.mkdir(self._rootFolder)

        # make the geometry folder if it doesn't exist
        if not os.path.exists(self._geoFolder):
            os.mkdir(self._geoFolder)

        # make the bc_dbase folder if it doesn't exist
        if not os.path.exists(self._bcFolder):
            os.mkdir(self._bcFolder)

        # make the simulation folder if it doesn't exist
        if not os.path.exists(self._simFolder):
            os.mkdir(self._simFolder)

        # make the results folder if it doesn't exist
        if not os.path.exists(self._resFolder):
            os.mkdir(self._resFolder)


        # change working directory to simulation folder
        os.chdir(self._simFolder)

        # create an empty list to store models
        models = list()

        # read grids from vector source
        grids = ()
        if self.gridSource is not None:
            grids = readGridsFromFile(self.gridSource)

        # read configuration from source
        config = SwanConfig()
        if self.configSource is not None:
            config.read(self.configSource)

        # convert time range to strings for naming runs
        tsString = datestr(timeStart, '%Y%m%d')
        teString = datestr(timeEnd, '%Y%m%d')

        # iterate over grids
        for grid in grids:
            # automatically name the model using the grid name and time range
            name = prefix + grid.getName() + '_' + tsString + '_' + teString + suffix

            # create a model for each grid
            model = SwanModel(modelName=name, swanGrid=grid, swanConfig=config, templateFile=self.templateSource)

            # set the model start and end times
            model.timeStart, model.timeEnd = timeStart, timeEnd

            # get control file and output file paths based on the model name
            nametmp = name + '.swn'
            model.controlFile = Path(os.path.relpath(self._simFolder)) / nametmp
            nametmp = name + '.nc'            
            model.outputFile = Path(os.path.relpath(self._resFolder)) / nametmp

            # append to list of models
            models.append(model)

        # automatically nest models in list
        autoNestModels(models)

        # write out bottom grid BC for each model
        for model in models:
            if all([os.path.isfile(f) for f in self.bottomSource]):
                print(model.timeStart)
                print(type(model.timeStart))
                nametmp = model.modelName + '_BOTTOM.txt'
                model.bottomFile = Path(os.path.relpath(self._geoFolder)) / nametmp
                writeBottomGridBc(self.bottomSource, model.bottomFile, model)

        if self.windSource is not None:
            if os.path.isfile(self.windSource):

                _, ext = os.path.splitext(os.path.split(self.windSource)[1])
                if ext == '.csv':
                    writer = writeWindTsBc
                elif ext == '.nc':                    
                    writer = writeWindGridBc

                for model in models:
                    nametmp = model.modelName + '_WIND.txt'
                    model.windFile = Path(os.path.relpath(self._bcFolder)) / nametmp
                    writer(self.windSource, model.windFile, model)

                    if not os.path.isfile(model.windFile):
                        model.windFile = None

        # write out wave spectra BC for outer most models
        if self.waveSource is not None:
            if os.path.isfile(self.windSource):
                for model in models:
                    if model.getParent() is None:
                        nametmp = model.modelName + '_SPECTRA.txt'
                        model.specFile = Path(os.path.relpath(self._bcFolder)) / nametmp
                        writeWaveGridBc(self.waveSource, model.specFile, model)

                        if not os.path.isfile(model.specFile):
                            model.specFile = None

        # write out control file for each model
        for model in models:
            model.writeControlFile()

        # write out build file
        self.write(Path(self._rootFolder) / 'BUILD.ini')

    def read(self, file):
        """Simple wrapper for ConfigParser.read() method"""

        # instantiate configparser object
        cp = configparser.ConfigParser()
        cp.optionxform = str

        # read build paths from file
        cp.read(file)

        parameters = ['_rootFolder', 'templateSource', 'configSource', 'gridSource',
                      'bottomSource', 'windSource', 'waveSource']

        for name in parameters:
            if name.replace('_', '') in cp['SWAN BUILD']:
                setattr(self, name, cp['SWAN BUILD'][name.replace('_', '')])

        if self.bottomSource is not None:
            self.bottomSource = self.bottomSource.split('\n')[1:]

        return self

    def write(self, file):
        """Simple wrapper for ConfigParser.write() method"""
        
        # instantiate ConfigParser object
        cp = configparser.ConfigParser()
        cp.optionxform = str

        # pass the build data to ConfigParser
        parameters = ['_rootFolder', 'templateSource', 'configSource', 'gridSource',
                      'bottomSource', 'windSource', 'waveSource']

        cp['SWAN BUILD'] = {}
        for name in parameters:
            value = getattr(self, name)
            if value is not None:
                cp['SWAN BUILD'][name.replace('_', '')] = str(value)

        cp['SWAN BUILD']['bottomSource'] = '\n\t' + '\n\t'.join(self.bottomSource)

        # save it to specified destination
        with open(file, 'w') as f:
            cp.write(f)
