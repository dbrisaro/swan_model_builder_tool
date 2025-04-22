"""Functions for configuration handling"""

import yaml
import os

def generate_config_ini(config):
    """Generate SWAN configuration file
    
    Parameters
    ----------
    config : dict
        Configuration parameters
    
    Returns
    -------
    str
        Configuration file content
    """
    config_str = f"""PROJECT '{config['project_name']}' '{config['case_name']}'
MODE STAT
SET level=0.0
CGRID {config['grid']['x_len']} {config['grid']['y_len']} {config['grid']['dx']} {config['grid']['dy']} {config['grid']['rotation']}
INPGRID BOTTOM {config['grid']['x_len']} {config['grid']['y_len']} {config['grid']['dx']} {config['grid']['dy']} {config['grid']['rotation']}
READINP BOTTOM 1 '{config['bathymetry_file']}' 3 0 FREE
"""
    return config_str

def create_directory_structure(base_dir):
    """Create directory structure for SWAN model
    
    Parameters
    ----------
    base_dir : str
        Base directory path
    
    Returns
    -------
    dict
        Dictionary with paths to created directories
    """
    dirs = {
        'input': os.path.join(base_dir, 'input'),
        'output': os.path.join(base_dir, 'output'),
        'grids': os.path.join(base_dir, 'grids'),
        'config': os.path.join(base_dir, 'config')
    }
    
    for dir_path in dirs.values():
        os.makedirs(dir_path, exist_ok=True)
    
    return dirs 