#!/usr/bin/env python3
# Project supervisor: Dr. Zohreh Kaheh.
# Undergraduate project collaborator: Jingxuan Hu.
# Initial student coding was developed with supervisory teaching and improvement guidance.
"""Summarise and plot the fixed-design multi-dataset recovery experiment."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata


PARAMETERS = ("N", "sigma", "beta0", "beta1")
LABELS = {"N": r"$N$", "sigma": r"$\sigma$", "beta0": r"$\beta_0$", "beta1": r"$\beta_1$"}


def split_chains(values: np.ndarray) -> np.ndarray:
    """Split m x n chains into 2m x floor(n/2) chains."""
    m, n = values.shape
    half = n // 2
    return np.concatenate((values[:, :half], values[:, -half:]), axis=0)


def basic_rhat(chains: np.ndarray) -> float:
    m, n = chains.shape
    chain_means = chains.mean(axis=1)
    within = np.mean(np.var(chains, axis=1, ddof=1))
    between = n * np.var(chain_means, ddof=1)
    if within == 0:
        return 1.0 if between == 0 else float("inf")
    variance = (n - 1.0) / n * within + between / n
    return float(np.sqrt(variance / within))


def rank_normalize(chains: np.ndarray) -> np.ndarray:
    flat = chains.ravel()
    ranks = rankdata(flat, method="average")
    transformed = norm.ppf((ranks - 0.375) / (len(flat) + 0.25))
    return transformed.reshape(chains.shape)


def rank_split_rhat(chains: np.ndarray) -> float:
    split = split_chains(chains)
    bulk = basic_rhat(rank_normalize(split))
    folded = np.abs(split - np.median(split))
    tail = basic_rhat(rank_normalize(folded))
    return max(bulk, tail)


def autocovariance_fft(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    n = len(x)
    size = 1 << (2 * n - 1).bit_length()
    transform = np.fft.rfft(x, size)
    acov = np.fft.irfft(transform * np.conjugate(transform), size)[:n]
    return acov / np.arange(n, 0, -1)


def bulk_ess(chains: np.ndarray) -> float:
    """Rank-normalised multi-chain ESS with Geyer's positive sequence."""
    chains = rank_normalize(split_chains(chains))
    m, n = chains.shape
    chain_var = np.var(chains, axis=1, ddof=1)
    within = chain_var.mean()
    between = n * np.var(chains.mean(axis=1), ddof=1)
    var_plus = (n - 1.0) / n * within + between / n
    if var_plus <= 0:
        return float(m * n)
    acov = np.array([autocovariance_fft(chain) for chain in chains])
    rho = np.ones(n)
    for lag in range(1, n):
        rho[lag] = 1.0 - (within - acov[:, lag].mean()) / var_plus
    pair_sums = []
    for lag in range(1, n - 1, 2):
        pair = rho[lag] + rho[lag + 1]
        if pair < 0:
            break
        pair_sums.append(pair)
    # Enforce Geyer's initial monotone sequence.
    for index in range(1, len(pair_sums)):
        pair_sums[index] = min(pair_sums[index], pair_sums[index - 1])
    tau = max(1.0, -1.0 + 2.0 * (1.0 + sum(pair_sums)))
    return float(min(m * n, m * n / tau))


def chain_matrix(group: pd.DataFrame, parameter: str) -> np.ndarray:
    chains = [part[parameter].to_numpy() for _, part in group.groupby("chain")]
    lengths = {len(chain) for chain in chains}
    if len(chains) < 2 or len(lengths) != 1:
        raise ValueError("diagnostics require at least two equal-length chains")
    return np.stack(chains)


def main() -> None:
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    results = root / "results" / "full_recovery"
    figures = root / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    traces = pd.read_csv(results / "recovery_traces.csv")
    acceptance = pd.read_csv(results / "recovery_acceptance.csv")
    design_frame = pd.read_csv(results / "simulation_design.csv")
    design = dict(zip(design_frame["quantity"], design_frame["value"]))
    truth = {
        "N": float(design["true_N"]),
        "sigma": float(design["true_sigma"]),
        "beta0": float(design["true_beta0"]),
        "beta1": float(design["true_beta1"]),
    }

    rows = []
    for dataset, group in traces.groupby("dataset"):
        for parameter in PARAMETERS:
            values = group[parameter].to_numpy()
            matrix = chain_matrix(group, parameter)
            q025, q975 = np.quantile(values, [0.025, 0.975])
            rows.append(
                {
                    "dataset": int(dataset),
                    "parameter": parameter,
                    "truth": truth[parameter],
                    "mean": values.mean(),
                    "sd": values.std(ddof=1),
                    "q025": q025,
                    "q975": q975,
                    "covered": int(q025 <= truth[parameter] <= q975),
                    "rank_split_rhat": rank_split_rhat(matrix),
                    "bulk_ess": bulk_ess(matrix),
                }
            )
    by_dataset = pd.DataFrame(rows)
    by_dataset["error"] = by_dataset["mean"] - by_dataset["truth"]
    by_dataset["interval_width"] = by_dataset["q975"] - by_dataset["q025"]
    by_dataset.to_csv(results / "recovery_by_dataset.csv", index=False)

    summary_rows = []
    for parameter, group in by_dataset.groupby("parameter", sort=False):
        summary_rows.append(
            {
                "parameter": parameter,
                "truth": truth[parameter],
                "mean_posterior_mean": group["mean"].mean(),
                "bias": group["error"].mean(),
                "rmse": np.sqrt(np.mean(group["error"] ** 2)),
                "coverage_95": group["covered"].mean(),
                "mean_interval_width": group["interval_width"].mean(),
                "median_rhat": group["rank_split_rhat"].median(),
                "max_rhat": group["rank_split_rhat"].max(),
                "median_bulk_ess": group["bulk_ess"].median(),
                "min_bulk_ess": group["bulk_ess"].min(),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(results / "recovery_summary.csv", index=False)

    diagnostics = {
        "datasets": int(traces["dataset"].nunique()),
        "chains_per_dataset": int(traces["chain"].nunique()),
        "retained_draws_per_chain": int(traces.groupby(["dataset", "chain"]).size().iloc[0]),
        "all_truths_covered": bool(by_dataset["covered"].all()),
        "overall_coverage": float(by_dataset["covered"].mean()),
        "max_rank_split_rhat": float(by_dataset["rank_split_rhat"].max()),
        "min_bulk_ess": float(by_dataset["bulk_ess"].min()),
        "max_cache_abs_error": float(acceptance["cache_max_abs_error"].max()),
        "acceptance_rate_ranges": {
            column: [float(acceptance[column].min()), float(acceptance[column].max())]
            for column in ("sigma", "beta0", "beta1", "location", "birth_death")
        },
        "diagnostic_definitions": {
            "rhat": "rank-normalised split R-hat, maximum of bulk and folded variants",
            "ess": "rank-normalised bulk ESS using Geyer's initial positive and monotone sequence",
        },
    }
    (results / "recovery_diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2) + "\n", encoding="utf-8"
    )

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6), constrained_layout=True)
    for axis, parameter in zip(axes.ravel(), PARAMETERS):
        subset = by_dataset[by_dataset["parameter"] == parameter]
        y = subset["dataset"].to_numpy()
        mean = subset["mean"].to_numpy()
        lower = mean - subset["q025"].to_numpy()
        upper = subset["q975"].to_numpy() - mean
        colors = ["#2471A3" if value else "#C0392B" for value in subset["covered"]]
        for yi, mi, lo, hi, color in zip(y, mean, lower, upper, colors):
            axis.errorbar(mi, yi, xerr=np.array([[lo], [hi]]), fmt="o", color=color,
                          capsize=2.5, markersize=4)
        axis.axvline(truth[parameter], color="black", linestyle="--", linewidth=1)
        axis.set_title(LABELS[parameter])
        axis.set_ylabel("Simulated dataset")
        axis.set_yticks(range(1, 7))
        axis.invert_yaxis()
    fig.suptitle("Posterior means and 95% credible intervals; dashed lines show generating values")
    fig.savefig(figures / "fig4_parameter_recovery.pdf", bbox_inches="tight")
    fig.savefig(figures / "fig4_parameter_recovery.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    representative = traces[traces["dataset"] == 1]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6), constrained_layout=True)
    for axis, parameter in zip(axes.ravel(), PARAMETERS):
        for chain, part in representative.groupby("chain"):
            axis.plot(part["draw"], part[parameter], linewidth=0.55, alpha=0.8,
                      label=f"Chain {chain}")
        axis.axhline(truth[parameter], color="black", linestyle="--", linewidth=1)
        diagnostics_row = by_dataset[(by_dataset["dataset"] == 1) &
                                     (by_dataset["parameter"] == parameter)].iloc[0]
        axis.set_title(f"{LABELS[parameter]}  R-hat={diagnostics_row['rank_split_rhat']:.3f}")
        axis.set_xlabel("Retained draw")
    axes[0, 0].legend(frameon=False, ncol=2)
    fig.suptitle("Trace diagnostics for representative simulated dataset 1")
    fig.savefig(figures / "fig5_trace_diagnostics.pdf", bbox_inches="tight")
    fig.savefig(figures / "fig5_trace_diagnostics.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    print(summary.to_string(index=False, float_format=lambda value: f"{value:.4f}"))
    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()
