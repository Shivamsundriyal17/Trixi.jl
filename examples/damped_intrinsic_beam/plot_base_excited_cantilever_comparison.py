#!/usr/bin/env python3
"""Plot the archived Farokhi branch comparison for inspection."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    example_directory = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--experimental",
        type=Path,
        default=(
            example_directory
            / "reference"
            / "farokhi_2022_experimental_frequency_response.csv"
        ),
    )
    parser.add_argument(
        "--comparison",
        type=Path,
        default=(
            example_directory
            / "reference"
            / "base_excited_cantilever_branch_comparison.csv"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            example_directory
            / "results"
            / "base_excited_cantilever_branch_comparison.png"
        ),
    )
    arguments = parser.parse_args()

    experimental = read_rows(arguments.experimental)
    comparison = read_rows(arguments.comparison)
    figure, axes = plt.subplots(
        2, 2, figsize=(8.0, 5.4), sharex=True, constrained_layout=True
    )

    for column, acceleration in enumerate((0.2, 0.5)):
        experiment = [
            row
            for row in experimental
            if float(row["acceleration_rms_g"]) == acceleration
        ]
        for component, row_index in (("transverse", 0), ("longitudinal", 1)):
            axis = axes[row_index, column]
            experimental_key = (
                "transverse_peak"
                if component == "transverse"
                else "longitudinal_minimum"
            )
            experimental_values = [
                abs(float(row[experimental_key])) for row in experiment
            ]
            axis.scatter(
                [float(row["normalized_frequency"]) for row in experiment],
                experimental_values,
                s=22,
                facecolors="none",
                edgecolors="#c43c39",
                linewidths=1.0,
                label="experiment",
                zorder=3,
            )

            for branch, linestyle in (("upper", "-"), ("lower", "--")):
                campaign_rows = sorted(
                    (
                        row
                        for row in comparison
                        if float(row["acceleration_rms_g"]) == acceleration
                        and row["branch"] == branch
                    ),
                    key=lambda row: float(row["normalized_frequency"]),
                )
                numerical_key = f"numerical_{component}"
                axis.plot(
                    [
                        float(row["normalized_frequency"])
                        for row in campaign_rows
                    ],
                    [abs(float(row[numerical_key])) for row in campaign_rows],
                    linestyle,
                    color="#2864a5",
                    linewidth=1.5,
                    marker=".",
                    markersize=4,
                    label=f"Trixi {branch}",
                    zorder=2,
                )

            axis.grid(True, linewidth=0.4, alpha=0.35)
            if row_index == 0:
                axis.set_title(f"{acceleration:.1f}g RMS")
            if column == 0:
                label = (
                    r"$\max |w_{\mathrm{tip}}|/L$"
                    if component == "transverse"
                    else r"$|\min u_{\mathrm{tip}}|/L$"
                )
                axis.set_ylabel(label)
            if row_index == 1:
                axis.set_xlabel(r"$f/f_1$")

    axes[0, 0].legend(frameon=False, fontsize=8)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(arguments.output, dpi=220)
    print(f"Wrote {arguments.output}")


if __name__ == "__main__":
    main()
