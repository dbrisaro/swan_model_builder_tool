import os
import sys
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, LineString
from datetime import datetime
import yaml
import configparser
from pathlib import Path

# Add the parent directory to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def datestr(date_str, format_str):
    """Convert date string or datetime object to another format.
    
    Args:
        date_str: String in format 'DD/MM/YYYY HH:MM:SS' or datetime/date object
        format_str: Output format string (e.g. '%Y%m%d')
    
    Returns:
        Formatted date string
    """
    if isinstance(date_str, str):
        try:
            dt = datetime.strptime(date_str, '%d/%m/%Y %H:%M:%S')
        except ValueError:
            try:
                dt = datetime.strptime(date_str, '%Y-%m-%d')
            except ValueError:
                dt = datetime.strptime(date_str, '%Y%m%d')
    else:
        if hasattr(date_str, 'hour'):  # datetime object
            dt = date_str
        else:  # date object
            dt = datetime.combine(date_str, datetime.min.time())
    return dt.strftime(format_str)

class SwanGrid:
    def __init__(self, name, polygon=None, rotation=None, dx=None, dy=None, nx=None, ny=None, nv=None, srs=None):
        self.name = name
        self.polygon = polygon
        self.rotation = rotation
        self.dx = dx
        self.dy = dy
        self.nx = nx
        self.ny = ny
        self.nv = nv
        self.srs = srs

    def getName(self):
        return self.name

    def setName(self, name):
        self.name = name

    def getPolygon(self):
        return self.polygon

    def setPolygon(self, polygon):
        self.polygon = polygon

    def getRotation(self):
        return self.rotation

    def setRotation(self, rotation):
        self.rotation = rotation

    def getDx(self):
        return self.dx

    def setDx(self, dx):
        self.dx = dx

    def getDy(self):
        return self.dy

    def setDy(self, dy):
        self.dy = dy

    def getNx(self):
        return self.nx

    def setNx(self, nx):
        self.nx = nx

    def getNy(self):
        return self.ny

    def setNy(self, ny):
        self.ny = ny

    def getNv(self):
        return self.nv

    def setNv(self, nv):
        self.nv = nv

    def getSrs(self):
        return self.srs

    def setSrs(self, srs):
        self.srs = srs

    def getW(self):
        return self.nx * self.dx

    def getH(self):
        return self.ny * self.dy

    def getX(self):
        return np.linspace(0, self.getW(), self.nx + 1)

    def getY(self):
        return np.linspace(0, self.getH(), self.ny + 1)

    def isDefined(self):
        return all([self.nx, self.ny, self.dx, self.dy])

    def isSpherical(self):
        return self.srs == 'EPSG:4326'

    def getFaces(self):
        if not self.isDefined():
            return None
        x = self.getX()
        y = self.getY()
        faces = []
        for i in range(self.nx):
            for j in range(self.ny):
                faces.append([
                    [x[i], y[j]],
                    [x[i+1], y[j]],
                    [x[i+1], y[j+1]],
                    [x[i], y[j+1]]
                ])
        return np.array(faces)

    def getNodes(self):
        if not self.isDefined():
            return None
        x = self.getX()
        y = self.getY()
        nodes = []
        for i in range(self.nx + 1):
            for j in range(self.ny + 1):
                nodes.append([x[i], y[j]])
        return np.array(nodes)

    def getBoundary(self):
        if not self.isDefined():
            return None
        x = self.getX()
        y = self.getY()
        return np.array([
            [x[0], y[0]],
            [x[-1], y[0]],
            [x[-1], y[-1]],
            [x[0], y[-1]]
        ])

    def getBoundingBox(self):
        if not self.isDefined():
            return None
        x = self.getX()
        y = self.getY()
        return {
            'xmin': x[0],
            'xmax': x[-1],
            'ymin': y[0],
            'ymax': y[-1]
        }

class SwanConfig:
    def __init__(self, level=0, dirNorth=90, depthMin=0.05, maxMessages=200, maxError=1, gravity=9.81, rho=1025,
                 maxDrag=0.0025, stationary=False, spherical=True, numDir=36, freqMin=0.04, freqMax=1.00, numFreq=34,
                 dAbsolute=0.005, dRelative=0.01, curvature=0.005, numPoints=99.5, maxIterations=5, limiter=0.01,
                 outputVars=None, outputPoints=None, outputType='.nc', timeStep=30, timeUnit='DAY'):
        self.level = level
        self.dirNorth = dirNorth
        self.depthMin = depthMin
        self.maxMessages = maxMessages
        self.maxError = maxError
        self.gravity = gravity
        self.rho = rho
        self.maxDrag = maxDrag
        self.stationary = stationary
        self.spherical = spherical
        self.numDir = numDir
        self.freqMin = freqMin
        self.freqMax = freqMax
        self.numFreq = numFreq
        self.dAbsolute = dAbsolute
        self.dRelative = dRelative
        self.curvature = curvature
        self.numPoints = numPoints
        self.maxIterations = maxIterations
        self.limiter = limiter
        self.outputVars = outputVars if outputVars else "XP YP HSIGN TPS PDIR DIR UBOT TMBOT FORCE DEPTH"
        self.outputPoints = outputPoints
        self.outputType = outputType
        self.timeStep = timeStep
        self.timeUnit = timeUnit

    def write(self, file):
        with open(file, 'w') as f:
            f.write('[SWAN CONFIG]\n')
            for attr, value in self.__dict__.items():
                if value is not None:
                    f.write(f"{attr} = {value}\n")

class SwanModel:
    def __init__(self, modelName, swanGrid, swanConfig=None, modelParent=None, timeStart=None, timeEnd=None,
                 controlFile=None, bottomFile=None, windFile=None, specFile=None, outputFile=None, templateFile=None):
        self.modelName = modelName
        self.swanGrid = swanGrid
        self.swanConfig = swanConfig
        self.modelParent = modelParent
        self.timeStart = timeStart
        self.timeEnd = timeEnd
        self.controlFile = controlFile
        self.bottomFile = bottomFile
        self.windFile = windFile
        self.specFile = specFile
        self.outputFile = outputFile
        self.templateFile = templateFile

    def getParent(self):
        return self.modelParent

    def setParent(self, parent):
        self.modelParent = parent

    def isStationary(self):
        return self.swanConfig.stationary if self.swanConfig else False

    def getStartUp(self):
        if self.modelParent:
            return f"INPGRID BOTTOM REGULAR {self.swanGrid.getNx()} {self.swanGrid.getNy()} {self.swanGrid.getDx()} {self.swanGrid.getDy()} {self.swanGrid.getBoundingBox()['xmin']} {self.swanGrid.getBoundingBox()['ymin']} {self.swanGrid.getRotation()}"
        return ""

    def getCGridInput(self):
        if self.swanConfig and self.swanConfig.spherical:
            return f"CGRID SPHERICAL {self.swanGrid.getNx()} {self.swanGrid.getNy()} {self.swanGrid.getDx()} {self.swanGrid.getDy()} {self.swanGrid.getBoundingBox()['xmin']} {self.swanGrid.getBoundingBox()['ymin']} {self.swanGrid.getRotation()}"
        return f"CGRID REGULAR {self.swanGrid.getNx()} {self.swanGrid.getNy()} {self.swanGrid.getDx()} {self.swanGrid.getDy()} {self.swanGrid.getBoundingBox()['xmin']} {self.swanGrid.getBoundingBox()['ymin']} {self.swanGrid.getRotation()}"

    def getBottomInput(self):
        if self.bottomFile:
            return f"INPGRID BOTTOM REGULAR {self.swanGrid.getNx()} {self.swanGrid.getNy()} {self.swanGrid.getDx()} {self.swanGrid.getDy()} {self.swanGrid.getBoundingBox()['xmin']} {self.swanGrid.getBoundingBox()['ymin']} {self.swanGrid.getRotation()}\nREADINP BOTTOM 1 '{self.bottomFile}' 3 0 FREE"
        return ""

    def getWindInput(self):
        if self.windFile:
            return f"INPGRID WIND REGULAR {self.swanGrid.getNx()} {self.swanGrid.getNy()} {self.swanGrid.getDx()} {self.swanGrid.getDy()} {self.swanGrid.getBoundingBox()['xmin']} {self.swanGrid.getBoundingBox()['ymin']} {self.swanGrid.getRotation()}\nREADINP WIND 1 '{self.windFile}' 3 0 FREE"
        return ""

    def getSurfaceInput(self):
        return ""

    def getCurrentInput(self):
        return ""

    def getSpecInput(self):
        if self.specFile:
            return f"BOUNDSPEC SEGMENT IJ 1 {self.swanGrid.getNy()} {self.swanGrid.getNx()} 1 {self.swanGrid.getNy()} 1 '{self.specFile}' 3 0 FREE"
        return ""

    def getNestInput(self):
        if self.modelParent:
            return f"NESTGRID {self.swanGrid.getNx()} {self.swanGrid.getNy()} {self.swanGrid.getDx()} {self.swanGrid.getDy()} {self.swanGrid.getBoundingBox()['xmin']} {self.swanGrid.getBoundingBox()['ymin']} {self.swanGrid.getRotation()}"
        return ""

    def getBoundaryInput(self):
        if self.modelParent:
            return "BOUNDNEST1 NEST 1"
        return ""

    def getNestOutput(self):
        if self.modelParent:
            return f"SAVENEST NEST 1 '{self.outputFile}'"
        return ""

    def getGridOutput(self):
        if self.outputFile:
            return f"OUTPUT OPTIONS BLOCK 1\nBLOCK 'COMPGRID' NOHEAD '{self.outputFile}' LAY 3 XP YP HSIGN TPS PDIR DIR UBOT TMBOT FORCE DEPTH"
        return ""

    def getSpecOutput(self):
        if self.outputFile:
            return f"OUTPUT OPTIONS BLOCK 1\nBLOCK 'COMPGRID' NOHEAD '{self.outputFile}' LAY 3 XP YP HSIGN TPS PDIR DIR UBOT TMBOT FORCE DEPTH"
        return ""

    def getOutput(self):
        if self.modelParent:
            return self.getNestOutput()
        return self.getGridOutput()

    def getNumeric(self):
        if self.swanConfig:
            return f"NUMERIC {self.swanConfig.dAbsolute} {self.swanConfig.dRelative} {self.swanConfig.curvature} {self.swanConfig.numPoints} {self.swanConfig.maxIterations} {self.swanConfig.limiter}"
        return ""

    def getComputeLockUp(self):
        if self.swanConfig and self.swanConfig.stationary:
            return "COMPUTE"
        return f"COMPUTE {self.timeStart} {self.timeEnd} {self.swanConfig.timeStep} {self.swanConfig.timeUnit}"

    def getControlString(self):
        control = []
        
        # Add startup if needed
        startup = self.getStartUp()
        if startup:
            control.append(startup)
        
        # Add CGRID
        control.append(self.getCGridInput())
        
        # Add bottom input
        bottom = self.getBottomInput()
        if bottom:
            control.append(bottom)
        
        # Add wind input
        wind = self.getWindInput()
        if wind:
            control.append(wind)
        
        # Add spec input
        spec = self.getSpecInput()
        if spec:
            control.append(spec)
        
        # Add nest input
        nest = self.getNestInput()
        if nest:
            control.append(nest)
        
        # Add boundary input
        boundary = self.getBoundaryInput()
        if boundary:
            control.append(boundary)
        
        # Add numeric
        numeric = self.getNumeric()
        if numeric:
            control.append(numeric)
        
        # Add compute
        control.append(self.getComputeLockUp())
        
        # Add output
        output = self.getOutput()
        if output:
            control.append(output)
        
        return "\n".join(control)

    def writeControlFile(self):
        if self.controlFile:
            with open(self.controlFile, 'w') as f:
                f.write(self.getControlString())

class SwanBuilder:
    def __init__(self, rootFolder='./', templateSource=None, configSource=None,
                 gridSource=None, bottomSource=None, windSource=None, waveSource=None):
        self.rootFolder = rootFolder
        self.templateSource = templateSource
        self.configSource = configSource
        self.gridSource = gridSource
        self.bottomSource = bottomSource
        self.windSource = windSource
        self.waveSource = waveSource

    def getRootFolder(self):
        return self.rootFolder

    def setRootFolder(self, rootFolder):
        self.rootFolder = rootFolder

    def buildRun(self, timeStart, timeEnd, prefix='', suffix=''):
        # Create output directory
        os.makedirs(self.rootFolder, exist_ok=True)
        
        # Create grid
        self.generate_grids()
        
        # Create config
        self.generate_config()
        
        # Request wave data
        self.request_wave_data()
        
        # Request wind data
        self.request_wind_data()
        
        # Create model
        model = SwanModel(
            modelName=f"{prefix}swan{suffix}",
            swanGrid=self.gridSource,
            swanConfig=self.configSource,
            timeStart=timeStart,
            timeEnd=timeEnd,
            controlFile=os.path.join(self.rootFolder, f"{prefix}swan{suffix}.swn"),
            bottomFile=os.path.join(self.rootFolder, "bottom.txt"),
            windFile=os.path.join(self.rootFolder, "wind.txt"),
            specFile=os.path.join(self.rootFolder, "spec.txt"),
            outputFile=os.path.join(self.rootFolder, f"{prefix}swan{suffix}.nc")
        )
        
        # Write control file
        model.writeControlFile()
        
        return model

    def write(self, file):
        with open(file, 'w') as f:
            f.write(f"rootFolder={self.rootFolder}\n")
            f.write(f"templateSource={self.templateSource}\n")
            f.write(f"configSource={self.configSource}\n")
            f.write(f"gridSource={self.gridSource}\n")
            f.write(f"bottomSource={self.bottomSource}\n")
            f.write(f"windSource={self.windSource}\n")
            f.write(f"waveSource={self.waveSource}\n")

def create_rectangular_grid(lon_min, lon_max, lat_min, lat_max, x_len, y_len, name='REGIONAL', rotation=0.0):
    """Create a rectangular grid for SWAN model.
    
    Args:
        lon_min: Minimum longitude
        lon_max: Maximum longitude
        lat_min: Minimum latitude
        lat_max: Maximum latitude
        x_len: Number of cells in x direction
        y_len: Number of cells in y direction
        name: Grid name
        rotation: Grid rotation in degrees
    
    Returns:
        SwanGrid object
    """
    # Create polygon
    polygon = Polygon([
        (lon_min, lat_min),
        (lon_max, lat_min),
        (lon_max, lat_max),
        (lon_min, lat_max)
    ])
    
    # Calculate grid dimensions
    dx = (lon_max - lon_min) / x_len
    dy = (lat_max - lat_min) / y_len
    
    # Create grid
    grid = SwanGrid(
        name=name,
        polygon=polygon,
        rotation=rotation,
        dx=dx,
        dy=dy,
        nx=x_len,
        ny=y_len,
        srs='EPSG:4326'
    )
    
    return grid


def convertMat2Nc(matFile, ncFile, isRotated, isNonSpherical):
    """Converts a SWAN result file from .mat to .nc format"""

    if not scipy_loaded:
        return 'scipy not loaded - if using QGIS, please check QGIS installation includes the "scipy" Python library via OSGeo4W installer'

    # we want to be able to visualize SWAN results in QGIS (compatible mesh file) but also be able to use
    # the results on TUFLOWFV in one-way nesting. This may not be possible.


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

    inData = loadmat(matFile)

    # pull list of keys from .mat file
    keys = inData.keys()

    # get the grid geometry
    x = inData.pop('Xp')
    y = inData.pop('Yp')

    # if grid is 1D, convert to 2D array for consistency
    if (x.ndim == 1) and (y.ndim == 1):
        x, y = np.meshgrid(x, y)

    # create output netCDF4 file
    ncOut = Dataset(ncFile, 'w')

    # set output netCDF4 file attributes
    ncOut.setncattr('History', 'SWAN GIS Tools Version 0.0.2')
    # MJS Consider adding cf convention global attr

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
        # skip magic keys
        if '__' in key:
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

    # close the netCDF4 file handle
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

from shapely.geometry import Polygon
import pandas as pd
import geopandas as gpd

def create_rectangular_grid(lon_min, lon_max, lat_min, lat_max, x_len, y_len, name='REGIONAL', rotation=0.0):
    """
    Creates a single rectangle as a GeoDataFrame with metadata.

    Parameters:
        lon_min (float): Minimum longitude
        lon_max (float): Maximum longitude
        lat_min (float): Minimum latitude
        lat_max (float): Maximum latitude
        x_len (float): X resolution in degrees (metadata only)
        y_len (float): Y resolution in degrees (metadata only)
        name (str): Name for the region
        rotation (float): Rotation (metadata only)

    Returns:
        gpd.GeoDataFrame: A GeoDataFrame with one polygon and metadata

    """

    # Define the 5-point closed polygon (ring)
    points = [
        (lon_min, lat_min),
        (lon_min, lat_max),
        (lon_max, lat_max),
        (lon_max, lat_min),
        (lon_min, lat_min)  # Close the polygon
    ]

    polygon = Polygon(points)

    df = pd.DataFrame({
        "Name": [name],
        "Rotation": [rotation],
        "X Length": [x_len],
        "Y Length": [y_len],
        "geometry": [polygon]
    })

    gdf = gpd.GeoDataFrame(df, geometry='geometry', crs="EPSG:4326")

    return gdf


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
    y = np.flip(y)

    nc_format = float(nc.getncattr('Conventions').split('-')[1])

    if nc_format > 1.8:
        t = nc['time'][:].data.astype(float)
        t = t + datenum((1970, 1, 1))

    else:
        t = nc['time'][:].data.astype(float)
        t = t * 3600 + datenum((1900, 1, 1))

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
        xt, yt = transformPoints(mb[:, 0], mb[:, 1], srsA, srsB)
        mb = np.column_stack((xt, yt))

    # get the model x, y, and t limits
    x1, x2 = mb[:, 0].min(), mb[:, 0].max()
    y1, y2 = mb[:, 1].min(), mb[:, 1].max()
    t1, t2 = datenum(swanModel.timeStart), datenum(swanModel.timeEnd)

    # find interpolation indices for model domain
    ix = getInterpolationIndex(x, x1, x2)
    iy = getInterpolationIndex(y, y1, y2)
    it = getInterpolationIndex(t, t1, t2)

    # if nothing in model domain, return none
    if any([len(ix) == 0, len(iy) == 0, len(it) < 2]):
        return None

    # subset dimensions
    ts, xs, ys = t[it], x[ix], y[iy]

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
            # read x & y components of wind vector
            u = np.flipud(nc['u10'][ii, :, :].data)
            v = np.flipud(nc['v10'][ii, :, :].data)

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
            zp1 = rasterInspect(inFile, x.flatten(), y.flatten(), srs)

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

def writeWaveGridBc(inFile, outFile, swanModel, spread=12):
    """ERA5 parametric wave boundary condition using TPAR files"""

    # get input wave data file netCDF4 handle
    nc = Dataset(inFile, 'r')

    # get longitude, latitude and time arrays
    x = nc['longitude'][:].data.astype(float)
    y = nc['latitude'][:].data.astype(float)

    nc_format = float(nc.getncattr('Conventions').split('-')[1])

    if nc_format > 1.8:
        t = nc['time'][:].data.astype(float)
        t = t + datenum((1970, 1, 1))

    else:
        t = nc['time'][:].data.astype(float)
        t = t * 3600 + datenum((1900, 1, 1))

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
    t1, t2 = datenum(swanModel.timeStart), datenum(swanModel.timeEnd)

    # find interpolation indices for model domain
    it = getInterpolationIndex(t, t1, t2)

    # if nothing in model domain, return none
    if len(it) < 2:
        return None

    # pre-allocate data for points
    hs = np.zeros((it.size, xp.size))
    mwd = np.zeros((it.size, xp.size))
    pwp = np.zeros((it.size, xp.size))

    # iterate over time indices
    for aa in range(len(it)):
        # get wave parameter data
        hs_ii = nc['swh'][it[aa], :, :].data.astype('float')
        mwd_ii = nc['mwd'][it[aa], :, :].data.astype('float')
        pwp_ii = nc['pp1d'][it[aa], :, :].data.astype('float')

    if nc_format <= 1.6:
        # Fill with zero values
        bad = (hs_ii == -32767) | \
              (mwd_ii == -32767) | \
              (pwp_ii == -32767)

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
        mwd_pp[mwd_pp < 0] = 360 + mwd_pp[mwd_pp < 0]

        # store the data
        hs[aa :] = hs_pp
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
    nx, ny = swanModel.swanGrid.getNx(), swanModel.swanGrid.getNy()

    controlString = \
        "BOUND SHAPESPEC PM PEAK DSPR DEGREES\n" + \
        "BOUNDSPEC SEGMENT IJ 0,0 {mx:},0 {mx:},{my:} 0,{my:} 0,0 VARIABLE FILE &\n".format(mx=nx - 1, my=ny - 1)

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
