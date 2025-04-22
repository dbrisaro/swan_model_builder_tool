"""Functions for data handling"""

import pandas as pd
from datetime import datetime, timedelta

def generate_date_lists(start_date, end_date, frequency):
    """Generate lists of dates for data requests
    
    Parameters
    ----------
    start_date : str
        Start date in format 'YYYY-MM-DD'
    end_date : str
        End date in format 'YYYY-MM-DD'
    frequency : str
        Frequency of data ('daily', 'monthly', etc.)
    
    Returns
    -------
    list
        List of dates
    """
    start = datetime.strptime(start_date, '%Y-%m-%d')
    end = datetime.strptime(end_date, '%Y-%m-%d')
    
    if frequency == 'daily':
        delta = timedelta(days=1)
    elif frequency == 'monthly':
        delta = timedelta(days=30)
    else:
        raise ValueError(f"Unsupported frequency: {frequency}")
    
    dates = []
    current = start
    while current <= end:
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
    return f"{variable}_{frequency}_{start_date}_{end_date}_{bounds['lon_min']}_{bounds['lon_max']}_{bounds['lat_min']}_{bounds['lat_max']}.nc" 