#!/bin/bash

# Set up environment
echo "Setting up environment..."
source .venv/bin/activate

# Step 1: Generate configuration
echo "Step 1: Generating configuration..."
python generate_config.py

# Step 2: Generate grids
echo "Step 2: Generating grids..."
python create_grid.py

# Step 3: Plot grids
echo "Step 3: Plotting grids..."
python plot_grids.py

# Step 4: Request data if needed
echo "Step 4: Requesting data..."
read -p "Do you need to request wind data? (y/n): " request_wind
if [ "$request_wind" = "y" ]; then
    python request_wind_data.py
fi

read -p "Do you need to request wave data? (y/n): " request_wave
if [ "$request_wave" = "y" ]; then
    python request_wave_data.py
fi

# Step 5: Build and run SWAN model
echo "Step 5: Building and running SWAN model..."
python build_and_run.py

# Step 6: Run SWAN simulations
echo "Step 6: Running SWAN simulations..."

# Read experiment configuration using Python
CONFIG_FILE="/Users/daniela/Documents/swan/swan_experiments/experiments_specs.txt"
OUTPUT_DIR="/Users/daniela/Documents/swan/swan_experiments/$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['output']['directory'])")"
SIM_DIR="$OUTPUT_DIR/SWAN/03_simulation"

# Change to simulation directory
cd "$SIM_DIR"

# Get grid names in order from config file
GRID_NAMES=($(python3 -c "import yaml; config = yaml.safe_load(open('$CONFIG_FILE')); print(' '.join([config['grids'][grid]['name'] for grid in config['grids']]))"))

echo "Found grid names in config: ${GRID_NAMES[@]}"

# First, print the files that will be executed
echo "Will execute the following .swn files:"
for grid_name in "${GRID_NAMES[@]}"; do
    swn_file=$(ls *"$grid_name"*.swn 2>/dev/null)
    if [ -n "$swn_file" ]; then
        echo "- $swn_file"
    else
        echo "Warning: No .swn file found for grid $grid_name"
    fi
done

# Ask for confirmation
read -p "Continue with execution? (y/n): " confirm
if [ "$confirm" != "y" ]; then
    echo "Execution cancelled"
    exit 0
fi

# Run simulations in order
for grid_name in "${GRID_NAMES[@]}"; do
    # Find the .swn file that contains this grid name
    swn_file=$(ls *"$grid_name"*.swn 2>/dev/null)
    if [ -n "$swn_file" ]; then
        echo "Running SWAN simulation for $swn_file..."
        swanrun -input "$swn_file"
    else
        echo "Warning: No .swn file found for grid $grid_name"
    fi
done

# Step 7: Process results
echo "Step 7: Processing results..."
cd ..
python process_results.py

echo "Workflow completed!" 