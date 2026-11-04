# PhysiCell Parameter Calibration from Experimental Cell Trajectories

MATLAB code to calibrate the parameters of a [PhysiCell](http://physicell.org/) cell-migration model against experimental cell trajectories. The workflow compares the time evolution of the mean position and the covariance of a cell population (experiment vs. simulation) and fits the model parameters with a non-linear least-squares optimizer (`lsqnonlin`).

> **Status:** research code. An experimental-computational framework to characterize immune cell migration in the glioblastoma microenvironment. https://doi.org/10.1016/j.cmpb.2026.109690

## How it works

1. **Experimental data processing** (`experimental_data_processor.m`)
   Reads the tracked cell trajectories, optionally generates a bootstrap sample (resampling cells with replacement), and computes at each time point the number of cells, the mean position (`Xmean`, `Ymean`) and the covariance (`Sxx`, `Syy`, `Sxy`). It also writes the initial cell positions as CSV files that PhysiCell uses as initial conditions.
2. **Quantity of interest (QoI)** (`MAIN_PhysiCell.m`)
   Builds a time-dependent weighted vector `[pmx.*Xmean; pmy.*Ymean; psx.*Sxx; psy.*Syy]`. The weights depend on the number of cells and on the variance at each time point.
3. **Optimization** (`MAIN_PhysiCell.m`, `loss_function_PhysiCell.m`)
   `lsqnonlin` (Levenberg-Marquardt) minimizes the residual between the experimental and simulated QoI over three parameters: `b`, `s` and `T` (see research paper https://doi.org/10.1016/j.cmpb.2026.109690, where 'b' is the migration bias, 's' the migration speed and 'T' the persistance time)
   At each iteration the loss function:
   - runs `Nreplicas` PhysiCell simulations through a bash script (`data_generator_PhysiCell.m`),
   - reads the simulated cell positions (`load_data_PhysiCell.m`),
   - computes the same statistics and weighted QoI for the simulations (`data_mean_processor_PhysiCell.m`, `quantity_of_interest_PhysiCell.m`),
   - averages over replicates and returns the residual.

## Repository structure

| File | Description |
|------|-------------|
| `MAIN_PhysiCell.m` | Main script: loads data, sets up and runs the optimization |
| `experimental_data_processor.m` | Processes experimental trajectories, bootstrap resampling, writes PhysiCell input CSVs |
| `loss_function_PhysiCell.m` | Objective function passed to `lsqnonlin` |
| `data_generator_PhysiCell.m` | Launches the PhysiCell simulations via a bash script |
| `load_data_PhysiCell.m` | Reads PhysiCell output files and builds the trajectory arrays |
| `data_mean_processor_PhysiCell.m` | Computes mean position and covariance per time point |
| `quantity_of_interest_PhysiCell.m` | Assembles the QoI vector from the statistics |

## Requirements

- MATLAB (R2025b) with the **Optimization Toolbox** (`lsqnonlin`)
- A working **PhysiCell** installation (see the [PhysiCell repository](https://github.com/MathCancer/PhysiCell))
- A Unix-like system (Linux or macOS), since simulations are launched with MATLAB's `unix()` command

## Installation

**The first thing to do is to install PhysiCell and everything it requires.** Only after that can the MATLAB code in this repository be used.

### 1. Install PhysiCell

PhysiCell is an open-source, C++ agent-based framework for multicellular simulations. Follow the official installation instructions:

- Repository and releases: <https://github.com/MathCancer/PhysiCell>
- Homepage: <http://physicell.org/>
- User Guide: `documentation/UserGuide.pdf` inside the PhysiCell folder

In short, you need a C++ compiler with OpenMP support (e.g. `g++`) and `make`. On macOS, you must also define the `PHYSICELL_CPP` environment variable pointing to an OpenMP-enabled `g++` (see the PhysiCell Quick Start for details).

Download or clone PhysiCell so that the folder is called `PhysiCell-master` and sits next to the Fitt_param_from_exp folder.

If you use PhysiCell, please cite it (see [Citation](#citation)).

### 2. Load and compile the user project `PBMConChip`

The PhysiCell model used here is the user project **`PBMConChip`**. Included in this repository under PhysiCell-master/user_projects/PBMConChip".

Place the project folder in `PhysiCell-master/user_projects/PBMConChip`, then, from inside the `PhysiCell-master` folder, load it and compile it:

```
make load PROJ=PBMConChip
make
```

`make load PROJ=PBMConChip` copies the project (its `Makefile`, `main.cpp`, `config/` and `custom_modules/`) into the PhysiCell working folders, and `make` compiles it and creates the `project` executable that the calibration code runs. If the project is distributed as a zip file, place `PBMConChip.zip` in `PhysiCell-master/user_projects/` and run `make unpack PROJ=PBMConChip` before loading it (available in PhysiCell 1.13.1 or later).

Recompile with `make` whenever you change the C++ code of the project.

### 3. Check that the project reads the correct experimental cells

The project must read the cell positions generated from **your** experimental trajectories, and the file names depend on the experimental condition (`cond`) and day (`D`) set in `MAIN_PhysiCell.m`. `experimental_data_processor.m` writes two files into `PhysiCell-master/config/`, named after those two values (e.g. `cond = 4`, `D = 1`):

| File written by MATLAB | Contents | Read by |
|------------------------|----------|---------|
| `cells_opt_c4_d1.csv` | `[t, x, y, z]` of every tracked cell, at the time it appears | `custom_modules/custom.cpp` |
| `initial_cells_opt_c4_d1.csv` | `[x, y, z, type]` of the cells present at `t = 0` | the XML configuration |

Before running anything, make sure that **all three places agree** on the same `cond` and `D`:

1. **`MAIN_PhysiCell.m`**: the values of `cond` and `D`.
2. **`custom_modules/custom.cpp`**: the line in `setup_tissue()` that loads the cell trajectories:
   ```cpp
   t_cell_positions = load_file_vector("./config/cells_opt_c4_d1.csv");
   ```
3. **The XML configuration** (`config/template5.xml`): the `cell_positions` block with the initial positions:
   ```xml
   <cell_positions type="csv" enabled="true">
       <folder>./config</folder>
       <filename>initial_cells_opt_c4_d1.csv</filename>
   </cell_positions>
   ```

If you change the line in `custom.cpp`, **recompile** the project with `make` (from `PhysiCell-master`). Changes in the XML do not require recompiling.

> `config/template5.xml` is a *template*, not a runnable configuration: it contains the placeholders `SESGO` (migration bias), `VELOCIDAD` (speed), `TPERS` (persistence time) and `SEMILLA` (random seed). `OPTphysicell.sh` replaces them with the values of each optimization step and replicate, and generates the `config/template.xml` that is actually passed to PhysiCell. Always edit `template5.xml`, never the generated `template*.xml` files.

### 4. Set up the simulation script

Copy `BashScripting/OPTphysicell.sh` from this repository into `PhysiCell-master/BashScripting/` and make sure it is executable. The script detects the PhysiCell directory from its own location (or from the `PHYSICELL_DIR` environment variable if set):

```bash
chmod +x PhysiCell-master/BashScripting/OPTphysicell.sh
```

## Expected folder layout

The code uses relative paths and assumes the following layout, where the MATLAB scripts are run from the Fitt_param_from_exp folder:

```
repository/
├── PhysiCell-master/
│   ├── BashScripting/
│   │   └── OPTphysicell.sh      # copy it from this repo's BashScripting/ folder
│   ├── config/                  # initial cell CSVs are written here
│   └── Results/                 # simulation outputs (OPT1, OPT2, ...)
└── Fitt_param_from_exp/
    ├── MAIN_PhysiCell.m
    ├── ...
    └── EXP_RESULTS/
        └── trajectories_collagen<cond>_day<D>/
            └── all_trajectories.mat
```

If your layout differs, edit the paths in `data_generator_PhysiCell.m`, `loss_function_PhysiCell.m`, `experimental_data_processor.m` and `MAIN_PhysiCell.m`.

Copy `BashScripting/OPTphysicell.sh` from this repository into `PhysiCell-master/BashScripting/` and make sure it is executable. The script detects the PhysiCell directory from its own location (or from the `PHYSICELL_DIR` environment variable if set):

```bash
chmod +x ../PhysiCell-master/BashScripting/OPTphysicell.sh
```

## Input data format

`all_trajectories.mat` must contain a variable `M`: a cell array with one entry per tracked cell. Each entry is a matrix with one row per time point and columns `[time, x, y]`.

## Usage

1. Set the experimental condition at the top of `MAIN_PhysiCell.m`:
   ```matlab
   cond = 4;  % collagen concentration: 2 or 4 mg/mL
   D    = 0;  % day 0 or day 1
   B    = 1;  % number of bootstrap samples
   ```
2. **First run:** uncomment the call to `experimental_data_processor(M,cond,D,bst)` inside the bootstrap loop in `MAIN_PhysiCell.m`. This generates `experimental_data_B###.mat` and the PhysiCell initial-condition CSVs. It can be commented out again on later runs to reuse the processed data.
3. Adjust the optimization settings if needed (initial guess, bounds, number of replicates `Nreplicas`, tolerances).
4. Run:
   ```matlab
   MAIN_PhysiCell
   ```

## Outputs

Generated in the working directory (and ignored by git):

- `experimental_data_B###.mat`: processed experimental statistics
- `PHYSICELL_DATA.mat`, `PHY_DATA.mat`, `OPT.mat`: simulated data and residuals from the last evaluation
- `results_c<cond>_d<D>_B###.mat`: full workspace saved after each bootstrap sample

## License

Released under the Apache License 2.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE).

## Citation

If you use this code, please cite: https://doi.org/10.1016/j.cmpb.2026.109690

This work relies on **PhysiCell**. Please also cite:

> Ghaffarizadeh, A., Heiland, R., Friedman, S. H., Mumenthaler, S. M., & Macklin, P. (2018). PhysiCell: An open source physics-based cell simulator for 3-D multicellular systems. *PLoS Computational Biology*, 14(2), e1005991. <https://doi.org/10.1371/journal.pcbi.1005991>

```bibtex
@article{ghaffarizadeh2018physicell,
  title   = {{PhysiCell}: An open source physics-based cell simulator for 3-D multicellular systems},
  author  = {Ghaffarizadeh, Ahmadreza and Heiland, Randy and Friedman, Samuel H. and Mumenthaler, Shannon M. and Macklin, Paul},
  journal = {PLoS Computational Biology},
  volume  = {14},
  number  = {2},
  pages   = {e1005991},
  year    = {2018},
  doi     = {10.1371/journal.pcbi.1005991}
}
```

PhysiCell is developed by the PhysiCell project and is distributed under its own license (see the PhysiCell repository).
