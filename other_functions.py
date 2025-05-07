#!/usr/bin/env python3
"""Functions for processing SWAN results"""

import os
import sys
import numpy as np
from netCDF4 import Dataset
from scipy.io import loadmat
from scipy.io.matlab import mio5_params
from functions.mdatetime import *
import warnings

def convertMat2Nc(matFile, ncFile, isRotated, isNonSpherical):
    """Converts a SWAN result file from .mat to .nc format"""

    # Store some dictionaries for attribute names
    attr_dict = {'Hsig_units': 'm',
                 'Hsig_long_name': 'significant wave height',
                 'Hsig_standard_name': 'sea surface_wave_significant_height',
                 'Dir_units': 'degrees',
                 'Dir_long_name': 'mean wave direction (cartesian going to)',
                 'Dir_standard_name': 'sea_surface_wave_cartesian_mean_to_direction',
                 'PkDir_units': 'degrees',
                 'PkDir_long_name': 'peak wave direction (cartesian going to)',
                 'PkDir_standard_name': 'sea_surface_wave_cartesian_peak_to_direction',
                 'TPsmoo_units': 's',
                 'TPsmoo_long_name': 'peak wave period',
                 'TPsmoo_standard_name': 'sea_surface_wave_smoothed_peak_period',
                 'Ubot_units': 'm s-1',
                 'Ubot_long_name': 'near bottom wave orbital velocity',
                 'Ubot_standard_name': 'sea_surface_wave_near_bottom_orbital_velocity',
                 'TmBot_units': 's',
                 'TmBot_long_name': 'near bottom wave period',
                 'TmBot_standard_name': 'sea_surface_wave_bottom_wave_period',
                 'WForce_x_units': 'N m-2',
                 'WForce_x_long_name': 'x-component of wave force',
                 'WForce_x_standard_name': 'sea_surface_wave_x_force_component',
                 'WForce_y_units': 'N m-2',
                 'WForce_y_long_name': 'y-component of wave force',
                 'WForce_y_standard_name': 'sea_surface_wave_y_force_component',
                 'Windv_x_units': 'm s-1',
                 'Windv_x_long_name': 'x-component of 10 metre wind',
                 'Windv_x_standard_name': '10m_u_wind_component',
                 'Windv_y_units': 'm s-1',
                 'Windv_y_long_name': 'y-component of 10 metre wind',
                 'Windv_y_standard_name': '10m_v_wind_component',
                 'Depth_units': 'm',
                 'Depth_long_name': 'depth',
                 'Depth_standard_name' : 'depth'
                 }

    # Suppress MatReadWarning about duplicate variables
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=UserWarning, module='scipy.io.matlab')
        inData = loadmat(matFile)

    # pull list of keys from .mat file
    keys = list(inData.keys())  # Convert to list to avoid view changes during iteration
    
    # Filter out metadata keys
    keys = [k for k in keys if not k.startswith('__')]

    # get the grid geometry
    # Check if Xp and Yp are in the keys before trying to pop them
    if 'Xp' in keys and 'Yp' in keys:
        x = inData.pop('Xp')
        y = inData.pop('Yp')
        # Remove from keys list to avoid processing them again
        keys.remove('Xp')
        keys.remove('Yp')
    else:
        # Try to find them in the time series data
        x_keys = [k for k in keys if k.startswith('Xp_')]
        y_keys = [k for k in keys if k.startswith('Yp_')]
        if x_keys and y_keys:
            x = inData[x_keys[0]]
            y = inData[y_keys[0]]
        else:
            raise ValueError("Could not find grid coordinates (Xp, Yp) in the .mat file")

    # if grid is 1D, convert to 2D array for consistency
    if (x.ndim == 1) and (y.ndim == 1):
        x, y = np.meshgrid(x, y)

    # create output netCDF4 file
    ncOut = Dataset(ncFile, 'w')

    # set output netCDF4 file attributes
    ncOut.setncattr('History', 'SWAN GIS Tools Version 0.0.2')

    # create dimensions
    ncOut.createDimension('Xp', x.shape[1])
    ncOut.createDimension('Yp', y.shape[0])
    ncOut.createDimension('time', 0)

    # create dimension variables
    if isRotated:
        ncOut.createVariable('Xp', np.float32, ('Yp', 'Xp',), fill_value=-999999)
        ncOut.createVariable('Yp', np.float32, ('Yp', 'Xp',), fill_value=-999999)
        ncOut['Xp'][:] = x
        ncOut['Yp'][:] = y
    else:
        ncOut.createVariable('Xp', np.float32, ('Xp',), fill_value=-999999)
        ncOut.createVariable('Yp', np.float32, ('Yp',), fill_value=-999999)
        ncOut['Xp'][:] = np.linspace(x[0, 0], x[0, -1], x.shape[1])
        ncOut['Yp'][:] = np.linspace(y[0, 0], y[-1, 0], y.shape[0])

    ncOut.createVariable('time', np.float64, ('time',), fill_value=-999999)

    if isNonSpherical:
        ncOut['Xp'].setncattr('units', 'm')
        ncOut['Xp'].setncattr('long_name', 'x coordinate of projection')
        ncOut['Yp'].setncattr('units', 'm')
        ncOut['Yp'].setncattr('long_name', 'y coordinate of projection')
    else:
        ncOut['Xp'].setncattr('units', 'degrees_east')
        ncOut['Xp'].setncattr('long_name', 'longitude')
        ncOut['Yp'].setncattr('units', 'degrees_north')
        ncOut['Yp'].setncattr('long_name', 'latitude')

    # Time attributes
    ncOut['time'].setncattr('units', 'seconds since 1970-01-01 00:00:00.0')
    ncOut['time'].setncattr('long_name', 'time')
    ncOut['time'].setncattr('calendar', 'gregorian')

    # write data to netCDF4
    ii = 0
    for key in keys:
        # Skip coordinate variables
        if key.startswith('Xp_') or key.startswith('Yp_'):
            continue

        # get parts of key
        parts = key.split('_')

        # first part is variable name
        name = '_'.join(parts[0:-2])

        # create new variable in nc if not yet created
        if name not in ncOut.variables:
            ncOut.createVariable(name, np.float32, ('time', 'Yp', 'Xp',), fill_value=-999999)

            # Set units and long names
            try:
                ncOut[name].setncattr('units', attr_dict[name + '_units'])
                ncOut[name].setncattr('long_name', attr_dict[name + '_long_name'])
            except KeyError:
                print('Can not find key: ' + name + '_long_name')

        # second and third is date and time
        timeString = parts[-2] + parts[-1]
        timeStamp = datenum(timeString, '%Y%m%d%H%M%S')

        # update time index if needed
        if (ii == 0) and (ncOut['time'].size == 0):
            ncOut['time'][ii] = timeStamp
        elif timeStamp > ncOut['time'][-1].data:
            ncOut['time'][ii + 1] = timeStamp
            ii += 1
        else:
            pass

        # pass data to nc file
        data = inData[key]
        data[np.isnan(data)] = -999999
        ncOut[name][ii] = data

    ncOut.close()




def extractSwanTs(inFile, outFile, pointsFile):
    """Extracts a time series from a .nc SWAN result file at specified points"""

    # get input netCDF4 file handle
    nc = Dataset(inFile, 'r')

    # get longitude, latitude and time arrays
    x, y = None, None
    if 'Xp' in nc.variables and 'Yp' in nc.variables:
        x = nc['Xp'][:].data.astype(float)
        y = nc['Yp'][:].data.astype(float)
    elif 'longitude' in nc.variables and 'latitude' in nc.variables:
        x = nc['longitude'][:].data.astype(float)
        y = nc['latitude'][:].data.astype(float)
    t = nc['time'][:].data.astype(float)

    # read the points file
    structure = {'names': ['NAME','X','Y'], 'formats': ['U50','f8','f8']}

    try:
        points = np.loadtxt(pointsFile, delimiter=',', dtype=structure, skiprows=1)
    except:
        raise ValueError('File format not supported or empty. Try header: NAME, X, Y')

    # Allow for no, single or multiple points
    if points.size == 1:
        points = np.repeat(points, repeats=2, axis=0)
        n = 1
    elif points.size == 0:
        raise ValueError('No points found in file or format not supported. Try header: NAME, X, Y')
    else:
        n = points.shape[0]

    # get index of nearest point
    i = np.zeros((n,), np.int16)
    j = np.zeros((n,), np.int16)
    xp = np.zeros((n,), np.float64)
    yp = np.zeros((n,), np.float64)
    for pp in range(n):
        if x.ndim == 1 and y.ndim == 1:
            j[pp] = np.argmin(np.abs(x - points[pp][1]))
            i[pp] = np.argmin(np.abs(y - points[pp][2]))
            xp[pp], yp[pp] = x[j[pp]], y[i[pp]]
        elif x.ndim == 2 and y.ndim == 2:
            jj = np.arange(x.shape[1])
            ii = np.arange(y.shape[0])
            jj, ii = np.meshgrid(jj, ii)
            jj, ii = jj.flatten(), ii.flatten()

            dx = x.flatten() - points[pp][1]
            dy = y.flatten() - points[pp][2]
            ds = np.hypot(dx, dy)

            nearest = np.argmin(ds)

            i[pp], j[pp] = ii[nearest], jj[nearest]
            xp[pp], yp[pp] = x[i[pp], j[pp]], y[i[pp], j[pp]]

    # pull out all time series of all time varying variables
    header, data, sub = 'TIME', [], np.arange(n)
    for var in nc.variables:
        if nc[var].ndim == 3:
            header = header + ',' + var.upper()
            data.append(nc[var][:, i, j].data[:, sub, sub])

    # remove netCDF4 file handle
    nc.close()

    # stack the data for each variable
    data = np.dstack(data)

    # get time as string
    ts = datestr(t).astype(object)

    # iterate over points
    for pp in range(n):
        # get output for point
        output = data[:, pp, :]

        # append the time string array
        output = np.column_stack((ts, output))

        # get output format for np.savetxt
        fmt = '%s' + ',' + ','.join(['%.2f' for _ in range(output.shape[1] - 1)])

        # append point name tag to outFile
        outFilePoint = outFile.replace('.csv', '_' + points[pp][0] + '.csv')

        # write CSV file
        with open(outFilePoint, 'w') as f:
            f.write('{:.6f},{:.6f}\n'.format(xp[pp], yp[pp]))
            f.write(header + '\n')
            np.savetxt(f, output, fmt)