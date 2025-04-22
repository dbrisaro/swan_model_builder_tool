"""Generate SWAN configuration files"""

import os
import yaml
from functions.config import generate_config_ini, create_directory_structure

def main():
    # Load configuration
    with open('experiments_specs.txt', 'r') as f:
        config = yaml.safe_load(f)
    
    # Create directory structure
    base_dir = f"/Users/daniela/Documents/swan/swan_experiments/{config['output']['directory']}"
    dirs = create_directory_structure(base_dir)
    
    # Generate SWAN configuration
    config_str = generate_config_ini(config)
    
    # Write configuration to file
    config_file = os.path.join(dirs['config'], 'swan_config.ini')
    with open(config_file, 'w') as f:
        f.write(config_str)
    
    print(f"Configuration file generated: {config_file}")

if __name__ == '__main__':
    main() 