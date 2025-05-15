"""Functions for data handling"""

import pandas as pd
from datetime import datetime, timedelta

def generate_date_lists(start_date, end_date, frequency):
    """Generate lists of dates for data requests
    
    Parameters
    ----------
    start_date : str or datetime.date
        Start date in format 'YYYY-MM-DD' or datetime.date object
    end_date : str or datetime.date
        End date in format 'YYYY-MM-DD' or datetime.date object
    frequency : str
        Frequency of data ('hourly', 'daily', 'monthly', etc.)
    
    Returns
    -------
    list
        List of dates
    """
    # Convert to datetime objects if they are strings
    if isinstance(start_date, str):
        start = datetime.strptime(start_date, '%Y-%m-%d')
    else:
        start = datetime.combine(start_date, datetime.min.time())
    
    if isinstance(end_date, str):
        end = datetime.strptime(end_date, '%Y-%m-%d')
    else:
        end = datetime.combine(end_date, datetime.min.time())
    
    if frequency == 'hourly':
        delta = timedelta(hours=1)
    elif frequency == 'daily':
        delta = timedelta(days=1)
    elif frequency == 'monthly':
        delta = timedelta(days=30)
    else:
        raise ValueError(f"Unsupported frequency: {frequency}")
    
    dates = []
    current = start
    while current <= end:
        if frequency == 'hourly':
            dates.append(current.strftime('%Y-%m-%d %H:%M:%S'))
        else:
            dates.append(current.strftime('%Y-%m-%d'))
        current += delta
    
    return dates

def generate_filename(variable, frequency, start_date, end_date, bounds):
    """Generate filename for data files
    
    Parameters
    ----------
    variable : str
        Variable name
    frequency : str
        Frequency of data
    start_date : str
        Start date
    end_date : str
        End date
    bounds : dict
        Bounding box coordinates
    
    Returns
    -------
    str
        Generated filename
    """
    # Asegurarse de que las fechas estén en formato YYYY-MM-DD y sin '/'
    start = str(start_date).replace('/', '-')
    end = str(end_date).replace('/', '-')
    # Formatear los límites con decimales
    lon_min = f"{bounds['lon_min']}"
    lon_max = f"{bounds['lon_max']}"
    lat_min = f"{bounds['lat_min']}"
    lat_max = f"{bounds['lat_max']}"
    return f"{variable}_{frequency}_{start}_{end}_{lon_min}_{lon_max}_{lat_min}_{lat_max}.nc" 