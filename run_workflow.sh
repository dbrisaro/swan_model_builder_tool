#!/bin/bash

# Set up environment
echo "Setting up environment..."
source .venv/bin/activate

# Create directory structure
echo "Creating directory structure..."
CONFIG_FILE="/Users/daniela/Documents/swan/swan_experiments/experiments_specs.txt"
OUTPUT_DIR="/Users/daniela/Documents/swan/swan_experiments/$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['output']['directory'])")"

# Create main directories
mkdir -p "$OUTPUT_DIR/DATA/BATHY"
mkdir -p "$OUTPUT_DIR/DATA/WIND"
mkdir -p "$OUTPUT_DIR/DATA/WAVE"
mkdir -p "$OUTPUT_DIR/QGIS"
mkdir -p "$OUTPUT_DIR/SWAN/01_geometry"
mkdir -p "$OUTPUT_DIR/SWAN/02_bc_dbase"
mkdir -p "$OUTPUT_DIR/SWAN/03_simulation"
mkdir -p "$OUTPUT_DIR/SWAN/04_results"

echo "Directory structure created successfully!"

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


# TODO: for gabo to be solved. 
# read -p "Do you need to request bathymetry data? (y/n): " request_bathymetry
# if [ "$request_bathymetry" = "y" ]; then
#     python request_bathymetry_data.py
# fi

# Step 5: Build and run SWAN model
echo "Step 5: Building and running SWAN model..."
python build_and_run.py

# Step 6: Run SWAN simulations
echo "Step 6: Running SWAN simulations..."
python run_simulations.py

# Step 7: Process results
echo "Step 7: Processing results..."
cd ..
python process_results.py

echo "Workflow completed!" 