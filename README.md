# Bayesian Spatial Modelling with C++ and Python

[Run and inputs](#run) · [Scope](#scope) · [Citation](#citation) · [Contributions](ATTRIBUTION.md) · [Reuse status](LICENSE_STATUS.md)

A native C++ sampling workflow with Python orchestration and diagnostics for Bayesian spatial modelling. The project demonstrates simulation-based computation, numerical integrity checks and repeatable execution across multiple datasets and chains.

**Availability:** Public source code for technical review and citation. See [publication and reuse terms](LICENSE_STATUS.md).

**Project design, supervision and subsequent development: Dr. Zohreh Kaheh.**  
**Initial implementation (undergraduate collaborator): Jingxuan Hu.**

The project builds on Dr Kaheh’s earlier research and supervised student collaboration. See [ATTRIBUTION.md](ATTRIBUTION.md) for the contribution details.

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

## Citation

Please cite the archived software and credit both contributors:

Kaheh, Z., & Hu, J. (2026). *Bayesian Spatial Modelling with C++ and Python* (Version snapshot-717b3bc) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23312496

Download [CITATION.cff](CITATION.cff) or [BibTeX](CITATION.bib). The [archived version](https://zenodo.org/records/23312496) is `snapshot-717b3bc`, from source commit [`717b3bc`](https://github.com/DR-ZohrehKaheh/bayesian-spatial-modelling/tree/717b3bc1e143587203edcdd56f2aeb8ae257c6f1). Cite the version or commit you used; later GitHub documentation may differ from the fixed archive. See [ZENODO_DEPOSIT.md](ZENODO_DEPOSIT.md) for the deposit record and checksums.

For reuse permissions, see [LICENSE_STATUS.md](LICENSE_STATUS.md).

## Questions and corrections

For technical questions or documentation corrections, [open an issue](https://github.com/DR-ZohrehKaheh/bayesian-spatial-modelling/issues) with the version or commit, environment, command and a minimal error report. Do not include private datasets, credentials or unpublished results. Discuss the scope and permissions of proposed code contributions with the maintainer first; [CONTRIBUTOR_RIGHTS.md](CONTRIBUTOR_RIGHTS.md) describes the required record.

## Continuous checks

[GitHub Actions](https://github.com/DR-ZohrehKaheh/bayesian-spatial-modelling/actions/workflows/checks.yml) runs the existing dataset-independent checks on Python 3.12 and Ubuntu. This covers utility tests or synthetic native/model checks, not the complete data-dependent research pipeline or scientific conclusions. The workflow does not publish generated research outputs.

