# Execution verification

**Project supervisor: Dr. Zohreh Kaheh.**

The released `run_pipeline.py` compiled and ran the synthetic native checks and the full six-dataset/two-chain native simulation. The supplied 96,000 iterations per chain were retained. Trace completeness, finite values, parameter support and intensity-cache checks passed, and the summary script completed.

Execution receipt: private log `native_final.log`, SHA-256 `8934c96f8ea30aa34487c8c05fc2dee87ae3e10d651441374d0d23c7d1359cc7`. No traces or generated plots are distributed.

The published verification scope covers the documented execution path and supplied inputs. Independent certification of research findings and validation across other inputs are outside that scope.

R/Rcpp interfaces and extended research experiments are outside this published edition and its verification scope.

## Environment

Python 3.12.14 on Linux. Native code was compiled with GCC 13.3.0 using C++17. Verified Python dependencies are pinned in `requirements.txt`: `matplotlib==3.10.8`, `numpy==2.3.5`, `pandas==2.2.3`, `scipy==1.17.0`.

## Version and archive

The execution record above applies to the documented software and inputs. The archived DOI version, source commit and checksums are recorded in [ZENODO_DEPOSIT.md](ZENODO_DEPOSIT.md). Later editorial updates do not rerun or extend the scientific validation.
