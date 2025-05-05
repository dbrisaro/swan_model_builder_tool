# SWAN Model Builder Tool

This tool helps automate the setup and execution of SWAN wave model simulations. It provides a structured workflow for generating configurations, creating grids, and running simulations.

## Directory Structure

The tool follows this directory structure:
```
swan_experiments/
├── run_climatology_southern_peru/  # Experiment directory
│   ├── DATA/
│   │   ├── BATHY/     # Bathymetry data
│   │   ├── WIND/      # Wind data
│   │   └── WAVE/      # Wave data
│   ├── QGIS/          # Grid files and plots
│   └── SWAN/          # SWAN model files
│       ├── 01_geometry/  # Grid geometry files
│       ├── 02_bc_dbase/  # Boundary condition files
│       ├── 03_simulation/ # SWAN input files
│       └── 04_results/   # Output files
```

## Configuration

The tool uses a YAML configuration file (`experiments_specs.txt`) that specifies:
- Grid definitions (regional and transition grids)
- Time settings (start/end dates, frequency)
- Output settings (directory, subdirectories)

## Workflow Steps

1. **Environment Setup**
   ```bash
   source .venv/bin/activate
   ```

2. **Generate Configuration**
   ```bash
   python generate_config.py
   ```
   Creates SWAN configuration files with model parameters.

3. **Generate Grids**
   ```bash
   python create_grid.py
   ```
   Creates the computational grids based on the configuration.

4. **Plot Grids**
   ```bash
   python plot_grids.py
   ```
   Generates visualizations of the grids.

5. **Request Data (Optional)**
   ```bash
   python request_wind_data.py
   python request_wave_data.py
   ```
   Downloads wind and wave data if needed.

6. **Build and Run SWAN Model**
   ```bash
   python build_and_run.py
   ```
   Generates SWAN input files and runs the simulations.

7. **Process Results**
   ```bash
   python process_results.py
   ```
   Processes SWAN results and generates visualizations:
   - Converts .mat files to NetCDF format
   - Generates daily mean significant wave height (Hsig) maps
   - Extracts and plots time series at specific points
   - Supports both regional and transition grids

## Running the Workflow

To run the entire workflow:
```bash
bash run_workflow.sh
```

The script will:
1. Set up the environment
2. Generate configurations
3. Create and plot grids
4. Optionally request wind and wave data
5. Build and run the SWAN model
6. Process the results

## Input Data Requirements

- **Bathymetry**: GEBCO data in .tif format
- **Wind Data**: ERA5 data in .nc format
- **Wave Data**: ERA5 data in .nc format
- **Grid Definition**: Shapefile with grid polygons
- **Points of Interest**: CSV file with coordinates (points.csv)
- **Analysis Period**: CSV file with dates (dates.csv)

## Output Files

The tool generates:
- Grid files in the QGIS directory
- SWAN input files (.swn) in the simulation directory
- NetCDF output files in the results directory
- Processed results in the processed_results directory:
  - Daily mean Hsig maps (4-column grid layout)
  - Time series plots for points of interest
  - Statistics and analysis results

## Results Processing

The results processor (`process_results.py`) requires:

1. `experiments_specs.txt`: YAML configuration file containing:
   - Base path for results
   - Grid configurations (bounds, resolution)
   - Output directory structure

2. `points.csv`: CSV file with points of interest:
   ```
   NAME,LON,LAT
   Paita,-81.1,-5.1
   Chicama,-79.5,-7.7
   Callao,-77.1,-12.0
   Ilo,-71.3,-17.6
   ```

3. `dates.csv`: CSV file with analysis period:
   ```
   START_DATE,END_DATE
   24/12/2024,03/01/2025
   ```

## Requirements

- Python 3.x
- Required packages:
  - numpy
  - matplotlib
  - xarray
  - pandas
  - cartopy
  - pyyaml
  - scipy

## Notes

- The tool supports various simulation types (climatology, storm events, etc.)
- Grids are processed in order (regional followed by transition)
- All paths are relative to the experiment directory
- The tool checks for required input files before running
- Results processing includes:
  - Automatic coordinate conversion (meters to lat/lon)
  - Consistent color scales across maps
  - Time series extraction at specified points
  - Daily mean calculations
- Dates should be in DD/MM/YYYY format
- Coordinates should be in decimal degrees

# SWAN Results Processor

This script processes SWAN model results, converting .mat files to NetCDF format and generating visualizations.

## Features

- Converts SWAN .mat files to NetCDF format
- Generates daily mean significant wave height (Hsig) maps
- Extracts and plots time series at specific points
- Supports both regional and transition grids
- Configurable through YAML and CSV files

## Requirements

- Python 3.x
- Required packages:
  - numpy
  - matplotlib
  - xarray
  - pandas
  - cartopy
  - pyyaml
  - scipy

## Input Files

1. `experiments_specs.txt`: YAML configuration file containing:
   - Base path for results
   - Grid configurations (bounds, resolution)
   - Output directory structure

2. `points.csv`: CSV file with points of interest:
   ```
   NAME,LON,LAT
   Paita,-81.1,-5.1
   Chicama,-79.5,-7.7
   Callao,-77.1,-12.0
   Ilo,-71.3,-17.6
   ```

3. `dates.csv`: CSV file with analysis period:
   ```
   START_DATE,END_DATE
   24/12/2024,03/01/2025
   ```

## Output

The script generates:

1. NetCDF files in the `processed_results` directory
2. Daily mean Hsig maps in the `processed_results/maps` directory
   - Maps are arranged in a 4-column grid
   - Single colorbar for all subplots
   - Latitude/longitude labels
   - Coastline and borders

3. Time series plots in the `processed_results/time_series` directory
   - Hsig time series for each point
   - Legend with point names
   - Grid for better readability

## Usage

1. Ensure all input files are in place:
   - `experiments_specs.txt`
   - `points.csv`
   - `dates.csv`

2. Run the script:
   ```bash
   python process_results.py
   ```

## Notes

- The script assumes SWAN results are in the `SWAN/04_results` directory
- All dates should be in DD/MM/YYYY format
- Coordinates in `points.csv` should be in decimal degrees
- The script automatically handles coordinate conversion from meters to lat/lon if needed 