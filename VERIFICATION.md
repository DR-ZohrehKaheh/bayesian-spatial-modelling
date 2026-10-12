# Execution verification

**Project supervisor: Dr. Zohreh Kaheh.**

The released `run_pipeline.py` compiled and ran the synthetic native checks and the full six-dataset/two-chain native simulation. The supplied 96,000 iterations per chain were retained. Trace completeness, finite values, parameter support and intensity-cache checks passed, and the summary script completed.

Execution receipt: private log `native_final.log`, SHA-256 `8934c96f8ea30aa34487c8c05fc2dee87ae3e10d651441374d0d23c7d1359cc7`. No traces or generated plots are distributed.

Verification covers the documented execution path and supplied inputs; it is not certification of research findings or all possible inputs.

R/Rcpp and extended research experiments are excluded and were not verified.

## Environment

Python 3.12.14 on Linux. Native code was compiled with GCC 13.3.0 using C++17. Verified Python dependencies are pinned in `requirements.txt`: `matplotlib==3.10.8`, `numpy==2.3.5`, `pandas==2.2.3`, `scipy==1.17.0`.

## Documentation review

Repository documentation was reviewed on 11 October 2026 for consistency with the included source files and recorded execution scope. This editorial review did not change executable code or rerun training. `SHA256SUMS.txt` records the current file contents.

## Rights and citation documentation — 12 October 2026

A documentation-only review added the reuse-status notice, contributor-rights register and software citation, and clarified the contribution account. Executable sources, dependencies and configuration are unchanged from commit `2c651dda46214d344f99fa4c9ae2923fe8c6d0b8`. No training or scientific validation was rerun for this update. Existing execution receipts do not establish publication or licensing authority. `SHA256SUMS.txt` was refreshed for the current files, excluding itself.

The subsequent clarification on 12 October 2026 records publication-only permission and the maintainer's acceptance of standard notice retention with a scholarly citation request. MIT remains conditional and was not added. This follow-up changes documentation and its checksums only.

DOI preparation added a BibTeX export and a Zenodo metadata/rights record on 12 October 2026. At that stage, these were preparation documents and no deposit or DOI had been created. Source code and the existing CFF citation were unchanged; documentation checksums were refreshed.

## DOI publication — 12 October 2026

Archived snapshot `717b3bc1e143587203edcdd56f2aeb8ae257c6f1` was published as `snapshot-717b3bc` in [Zenodo](https://zenodo.org/records/23312496), DOI `10.5281/zenodo.23312496`. The public record, both creator names, archive checksum and custom rights notice were checked. Citation files and README now identify this archived version. The original archive retains its pre-DOI citation files. Source code, dependency pins and scientific execution scope are unchanged; no training was rerun. See [ZENODO_DEPOSIT.md](ZENODO_DEPOSIT.md). Documentation checksums were refreshed.

## Portfolio presentation review — 12 October 2026

The README navigation, contribution summary, citation guidance and issue-reporting instructions were reviewed for clarity and consistency with the existing rights record. Documentation links and citation identifiers were checked. Source code, dependencies, configuration and the archived DOI version are unchanged; no training or scientific validation was rerun. `SHA256SUMS.txt` was refreshed for the current files.
