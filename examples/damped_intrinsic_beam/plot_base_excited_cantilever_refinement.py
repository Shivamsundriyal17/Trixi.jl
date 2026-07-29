#!/usr/bin/env python3
"""Visualize response and work sensitivity to spatial refinement."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    example_directory = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=(
            example_directory
            / "reference"
            / "base_excited_cantilever_refinement_comparison.csv"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            example_directory
            / "results"
            / "base_excited_cantilever_refinement_diagnostics.png"
        ),
    )
    arguments = parser.parse_args()

    rows = sorted(
        (
            row
            for row in read_rows(arguments.input)
            if row["campaign"] == "05g_upper"
        ),
        key=lambda row: float(row["normalized_frequency"]),
    )
    frequencies = np.array(
        [float(row["normalized_frequency"]) for row in rows]
    )

    figure, axes = plt.subplots(
        1, 3, figsize=(11.2, 3.5), constrained_layout=True
    )
    response_series = (
        ("transverse_change_fraction", r"$w_{\mathrm{tip}}$", "#2864a5"),
        (
            "longitudinal_change_fraction",
            r"$u_{\mathrm{tip}}$",
            "#c43c39",
        ),
        ("rotation_change_fraction", r"$\psi_{\mathrm{tip}}$", "#2f8f5b"),
    )
    for key, label, color in response_series:
        axes[0].plot(
            frequencies,
            100.0 * np.abs([float(row[key]) for row in rows]),
            marker="o",
            markersize=4,
            linewidth=1.5,
            color=color,
            label=label,
        )
    axes[0].set(
        xlabel=r"$f/f_1$",
        ylabel="absolute refinement change [%]",
        title="(a) Response sensitivity",
    )
    axes[0].grid(True, linewidth=0.4, alpha=0.3)
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].axhline(0.0, color="0.6", linewidth=0.8)
    axes[1].plot(
        frequencies,
        100.0
        * np.array(
            [float(row["baseline_numerical_fraction"]) for row in rows]
        ),
        color="#c43c39",
        marker="o",
        markersize=4,
        linewidth=1.5,
        label="baseline numerical",
    )
    axes[1].plot(
        frequencies,
        100.0
        * np.array(
            [float(row["refined_numerical_fraction"]) for row in rows]
        ),
        color="#2864a5",
        marker="D",
        markersize=4,
        linewidth=1.5,
        label="refined numerical",
    )
    axes[1].plot(
        frequencies,
        100.0
        * np.array(
            [float(row["baseline_material_fraction"]) for row in rows]
        ),
        color="#c43c39",
        linestyle=":",
        linewidth=1.2,
        label="baseline Kelvin–Voigt",
    )
    axes[1].plot(
        frequencies,
        100.0
        * np.array(
            [float(row["refined_material_fraction"]) for row in rows]
        ),
        color="#2864a5",
        linestyle=":",
        linewidth=1.2,
        label="refined Kelvin–Voigt",
    )
    axes[1].set(
        xlabel=r"$f/f_1$",
        ylabel="signed fraction of root work [%]",
        title="(b) Work-partition sensitivity",
    )
    axes[1].grid(True, linewidth=0.4, alpha=0.3)
    axes[1].legend(frameon=False, fontsize=7)

    peak = rows[-1]
    labels = ["baseline", "refined"]
    material = 100.0 * np.array(
        [
            float(peak["baseline_material_fraction"]),
            float(peak["refined_material_fraction"]),
        ]
    )
    interface = 100.0 * np.array(
        [
            float(peak["baseline_jump_fraction"]),
            float(peak["refined_jump_fraction"]),
        ]
    )
    boundary = 100.0 * np.array(
        [
            float(peak["baseline_net_boundary_fraction"]),
            float(peak["refined_net_boundary_fraction"]),
        ]
    )
    positions = np.arange(2)
    axes[2].bar(
        positions,
        material,
        color="#2f8f5b",
        width=0.62,
        label="Kelvin–Voigt",
    )
    axes[2].bar(
        positions,
        interface,
        bottom=material,
        color="#e49c3f",
        width=0.62,
        label="interface",
    )
    axes[2].bar(
        positions,
        boundary,
        bottom=material + interface,
        color="#735aa8",
        width=0.62,
        label="net boundary",
    )
    numerical = interface + boundary
    for index, value in enumerate(numerical):
        axes[2].text(
            index,
            material[index] + 0.5 * value,
            f"{value:.1f}%\nnumerical",
            ha="center",
            va="center",
            fontsize=8,
            color="white",
            fontweight="bold",
        )
    axes[2].set(
        xticks=positions,
        xticklabels=labels,
        ylabel="fraction of root work [%]",
        ylim=(0.0, 105.0),
        title=rf"(c) Peak split at $f/f_1={frequencies[-1]:.5f}$",
    )
    axes[2].grid(True, axis="y", linewidth=0.4, alpha=0.3)
    axes[2].legend(frameon=False, fontsize=7, loc="lower left")

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(arguments.output, dpi=240)
    figure.savefig(arguments.output.with_suffix(".pdf"))
    plt.close(figure)
    print(f"Wrote {arguments.output}")
    print(f"Wrote {arguments.output.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
