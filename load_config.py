import yaml

def load_run_config(config_path, run_name):
    """Load the config for a specific run by name from the YAML file."""
    with open(config_path, "r") as f:
        all_config = yaml.safe_load(f)
    for run in all_config['runs']:
        if run['name'] == run_name:
            return run
    raise ValueError(f"Run '{run_name}' not found in config.")

def main(config_path, run_name):
    config = load_run_config(config_path, run_name)
    return config