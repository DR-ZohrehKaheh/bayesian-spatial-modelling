# Bayesian Spatial Modelling with C++ and Python

A native C++ sampling workflow with Python orchestration and diagnostics for Bayesian spatial modelling. The project demonstrates simulation-based computation, numerical integrity checks and repeatable execution across multiple datasets and chains.

**Project supervisor: Dr. Zohreh Kaheh.**  
**Undergraduate collaborator: Jingxuan Hu.**

Dr. Kaheh provided instruction in the project concepts and methods, supervised implementation and guided improvements to the student’s initial code.

## Technical focus

- Native C++ implementation of a Bayesian sampling workflow.
- Synthetic data generation and repeated simulation runs.
- Python orchestration, trace summaries and diagnostics.
- Checks for trace completeness, parameter support and cached calculations.

## Run

The native workflow was verified on Linux with GCC 13.3.0 and Python 3.12.14. The wrapper expects `g++` on PATH; use `--compiler` to select a compatible compiler. Other operating systems and compilers have not been verified. Install `requirements.txt` in an isolated environment, then run from this directory:

```bash
python -m pip install -r requirements.txt
python run_pipeline.py --output outputs
```

The wrapper compiles and executes the native checks, simulates input data, runs the fixed native sampling workflow, summarises the generated traces and checks their completeness and numeric integrity. It uses six synthetic datasets, two chains per dataset and 96,000 iterations per chain, with the supplied burn-in and thinning settings. The verification run used this complete native configuration.

Source files:

- `bd_secr_core.hpp`: sampler and input validation.
- `test_bd_secr_core.cpp`: synthetic native checks.
- `run_multidataset_recovery.cpp`: data simulation and sampling.
- `summarize_recovery.py`: summaries and diagnostics of generated outputs.
- `run_pipeline.py`: compilation, execution and integrity checks.

## Scope

R/Rcpp files are excluded because that interface was not executable in the verification environment. The larger extended research experiment is not part of this release and is not claimed to have been rerun. Passing execution checks is not proof of identifiability, convergence for arbitrary data or scientific validity. No recorded research traces or paper results are distributed.

## Verification and contributions

See [VERIFICATION.md](VERIFICATION.md) for the tested execution scope and [ATTRIBUTION.md](ATTRIBUTION.md) for supervision and student contributions. This repository contains code and technical documentation only.
