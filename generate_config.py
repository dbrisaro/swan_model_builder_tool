"""Generate SWAN configuration files"""

from pathlib import Path
import yaml

def generate_swan_config():
    """Generate SWAN configuration file in the correct format
    
    Returns
    -------
    str
        Configuration file content in the correct format
    """
    config_str = """[SWAN CONFIG]
level = 0
dirNorth = 90
depthMin = 0.05
maxMessages = 200
maxError = 1
gravity = 9.81
rho = 1025
maxDrag = 0.0025
stationary = False
spherical = True
numDir = 36
freqMin = 0.04
freqMax = 1.0
numFreq = 34
dAbsolute = 0.005
dRelative = 0.01
curvature = 0.005
numPoints = 99.5
maxIterations = 5
limiter = 0.01
outputVars = XP YP HSIGN TPS PDIR DIR UBOT TMBOT FORCE DEPTH
outputType = .mat
timeStep = 3600
timeUnit = SEC
"""
    return config_str

def main(config):

    base_path = Path(config['base']['path'])
    output_dir = base_path / config['output']['directory'] / 'SWAN'
    
    # Generate SWAN configuration
    config_str = generate_swan_config()
    
    # Write configuration to file
    config_file = output_dir / 'CONFIG.ini'
    with open(config_file, 'w') as f:
        f.write(config_str)
    
    print(f"Configuration file generated: {config_file}")

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run SWAN config generator.')
    parser.add_argument('--config', type=str, required=True, help='Path to the experiment specifications file')
    args = parser.parse_args()

    main(args.config)