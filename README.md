# Non-natural base pairs in TMSD

## Project Overview

This repository contains the analysis environment for the _Systematic incorporation of non-natural base pairs for increased toehold-mediated strand displacement kinetics in random DNA pools_. It houses the curated time-course datasets, Jupyter notebooks used for fitting/reporting and generating figures that appear in the manuscript. All workflows revolve around reproducing the cascade kinetics comparison between UW and TUM experimental contexts.

## Repository Map

- `Data/` – CSV exports of calibrated plate-reader data (UW_* and TUM_* groups)
- `src/` – Reusable Python modules (`modeling.py`, `statistics.py`, `utilities.py`)
- `*.ipynb` – Notebook workflows (model fitting, figure generation, parameter sweeps)
- `README.md` – You are here

## Raw Data

The `Data/` directory is organized by site and experiment. Experiments on the 5' toehold systems are prepended with TUM, and those on the 3' toehold systems, multi isoC, and long overhang systems are prepended with UW. The contents of indidivual datasets for each tag (TUM/UW):

### TUM Data
- `CascadeLeak.csv` - Used for plotting the 5' cacade system in the absence of trigger
- `Cascades_CascIncub.csv` - Kinetics for the 5' cascade system with the cascade pre-incubated in random background
- `Cascades_NoIncub.csv` - Kinetics for the 5' cascade system with the cascade directly added to random background
- `Cascades_TrigIncub.csv` - Kinetics for the 5' cascade system with the invader pre-incubated in random background
- `controls.csv` - Baseline measurements and  kinetics 5' one-step system in the absence of background
- `NN_Mismatches_InvIsoC.csv` - 5' one-step system with mismatched between natural bases and isoC in the invader
- `NN_Mismatches_RepIsoG.csv` - 5' one-step system with mismatched between natural bases and isoG in the reporter
- `no.csv` - Kinetics for the 5' one-step system with the system directly added to random background
- `rep.csv` - Kinetics for the 5' one-step system with the reporter pre-incubated in random background
- `trig.csv` - Kinetics for the 5' one-step system with the invader pre-incubated in random background

### UW Data
- `backgr.csv` - Background readings for the 3' one-step reporter without trigger and different combinations of background
- `CascadeModel_Natural_Trial1.csv` One replicate trial for the natural 3' toehold cascade used to fit a detailed kinetic model
- `CascadeModel_Natural_Trial2.csv` A second replicate trial for the natural 3' toehold cascade used to fit a detailed kinetic model
- `CascadeModel_NN.csv` - Data for fitting a detailed kinetic model to the non-natural 3' toehold cascade
- `Cascades.csv` - Kinetics for the 3' toehold cascades across all incubation conditions
- `LongOH_NoInc` - Kinetics for the long overhang systems added directly to background
- `LongOH_RepInc1` - Kinetics for the long overhang systems with the reporter and gate pre-incubated in background
- `Multiiso_NN2_close_NoInc.csv` - Kinetics for the multiple isoC experiments with 2 isoC in the close spacing with the system added directly to the background
- `Multiiso_NN2_close_TrigInc.csv` - Kinetics for the multiple isoC experiments with 2 isoC in the close spacing with the invader preincubated in background
- `Multiiso_NoInc.csv` - Kinetics for the multiple isoC experiments with either 1 isoC or 2 isoC in the far spacing with the system added directly to the background
- `Multiiso_TrigInc.csv` - Kinetics for the multiple isoC experiments with either 1 isoC or 2 isoC in the far spacing with the invader preincubated in background
- `no.csv` - Kinetics for the 3' one-step system with the system directly added to the background
- `trig.csv` - Kinetics for the 3' one-step system with the invader preincubated in background
- `rep.csv` - Kinetics for the 3' one-step system with the reporter preincubated in background
- `RepNat_Kinetics.csv` - Data for fitting a precise rate constant to the natural 3' toehold reporter
- `RepNN__MM_Kinetics.csv` - Data for fitting a precise rate constant to the non-natural 3' toehold reporter when the mismatching natural invader was used
- `RepNN__NoMM_Kinetics.csv` - Data for fitting a precise rate constant to the non-natural 3' toehold reporter with the complementary non-natural invader

## Analysis Notebooks

- `cacade_modeling.ipynb` - Fits second order rate models to all cascade data
- `cascade_leak.ipynb` – Analyzes leak for both the 5' and 3' teohold cascade systems
- `long_overhang_figure.ipynb` – Analyzes data from the long overhang experiment and generates Figure 8B
- `multi_iso_figure.ipynb` – Analyzes data from the multi-isoC experiment and generates Figure 8A
- `mismatch_ddg_calculations.ipynb` - Calculates ∆∆G for mismatches in all trinucleotide contexts
-  `nonnatural_mismatches.ipynb` – Explores kinetics of single natural-non-natural mismatches
- `multi_parameter_model.ipynb` – Fits a secondary supportive mechanistic model to the 5' toehold one-step system no incubation data
- `nat_vs_nn_detailed.ipynb` – Compares natural vs non-natural single step reactions and cascades with detailed kinetic modeling.
- `one_step_modeling.ipynb` – Fits second order rate models to all one-step system data

Each notebook expects the CSVs in `Data/` and imports helpers from `src/`. Run them top-to-bottom to regenerate intermediate tables and SVGs.

## Helper Python Modules (`src/`)

- `modeling.py` – Defines ODE-based kinetic models (`model_one_step`, `model_two_step`, `model_cascade`, etc.) plus fitting utilities (`residuals`, `fit_model`).
- `statistics.py` – Provides statistical tests (Wald Z, Holm correction, ratio Z-tests) for reporting confidence in fitted parameters.
- `utilities.py` – Shared helpers for time parsing, calibration, averaging, uncertainty propagation, and plotting conveniences.
