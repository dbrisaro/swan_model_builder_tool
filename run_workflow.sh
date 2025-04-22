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

echo "Workflow completed!" 