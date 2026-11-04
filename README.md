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
4. **Visualization**
   `MAIN_PhysiCell.m` plots experimental vs. computational mean positions and standard deviations over time.

## Repository structure

| File | Description |
|------|-------------|
| `MAIN_PhysiCell.m` | Main script: loads data, sets up and runs the optimization, plots results |
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

https://doi.org/10.1016/j.cmpb.2026.109690
