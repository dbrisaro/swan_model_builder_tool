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

7. **Process Results (TODO)**
   ```bash
   python process_results.py
   ```
   TODO: Implement results processing script that will:
   - [ ] Read output NetCDF files
   - [ ] Calculate statistics and climatologies
   - [ ] Generate plots and visualizations
   - [ ] Create summary reports
   - [ ] Integrate with workflow steps

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
6. Process the results (once implemented)

## Input Data Requirements

- **Bathymetry**: GEBCO data in .tif format
- **Wind Data**: ERA5 data in .nc format
- **Wave Data**: ERA5 data in .nc format
- **Grid Definition**: Shapefile with grid polygons

## Output Files

The tool generates:
- Grid files in the QGIS directory
- SWAN input files (.swn) in the simulation directory
- NetCDF output files in the results directory
- Plots and visualizations
- Processed results and statistics (once implemented)

## Configuration Parameters

Key parameters in the SWAN configuration:
- `timeStep`: Time step for output (configurable)
- `timeUnit`: Time unit for output (configurable)
- `outputVars`: Variables to output (HSIGN, TPS, PDIR, etc.)
- `outputType`: Output file format (default: .nc)

## Notes

- The tool supports various simulation types (climatology, storm events, etc.)
- Grids are processed in order (regional followed by transition)
- All paths are relative to the experiment directory
- The tool checks for required input files before running
- Results processing will be integrated into the workflow once implemented
- Time step and frequency can be adjusted based on simulation needs 