#!/usr/bin/env python3
"""Compare the 0.5g cantilever response with and without dead gravity."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


PERIODICITY_TOLERANCE = 5.0e-4


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def read_key_value(path: Path, key: str) -> float:
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            name, separator, value = line.partition("=")
            if separator and name.strip() == key:
                return float(value.strip())
    raise ValueError(f"{key!r} not found in {path}")


def read_numerical_check(path: Path, quantity: str) -> float:
    for row in read_rows(path):
        if row["check"] == quantity:
            return float(row["value"])
    raise ValueError(f"{quantity!r} not found in {path}")


def matching_rows(
    gravity_rows: list[dict[str, str]],
    zero_gravity_rows: list[dict[str, str]],
) -> list[tuple[dict[str, str], dict[str, str]]]:
    pairs = []
    for zero_row in zero_gravity_rows:
        frequency = float(zero_row["normalized_frequency"])
        candidates = [
            row
            for row in gravity_rows
            if abs(float(row["normalized_frequency"]) - frequency) < 1.0e-12
        ]
        if len(candidates) != 1:
            raise ValueError(
                f"expected one gravity-on row at normalized frequency "
                f"{frequency}, found {len(candidates)}"
            )
        pairs.append((candidates[0], zero_row))
    return pairs


def first_nonperiodic_bracket(
    rows: list[dict[str, str]],
) -> tuple[dict[str, str], dict[str, str]]:
    for index, row in enumerate(rows):
        if float(row["periodicity_error"]) > PERIODICITY_TOLERANCE:
            if index == 0:
                raise ValueError("first response point is not periodic")
            previous = rows[index - 1]
            if (
                float(previous["periodicity_error"])
                > PERIODICITY_TOLERANCE
            ):
                raise ValueError("nonperiodic response has no accepted predecessor")
            return previous, row
    raise ValueError("no nonperiodic transition point found")


def write_comparison(
    path: Path,
    pairs: list[tuple[dict[str, str], dict[str, str]]],
) -> None:
    fieldnames = [
        "normalized_frequency",
        "gravity_frequency_hz",
        "zero_gravity_frequency_hz",
        "gravity_transverse_peak",
        "zero_gravity_transverse_peak",
        "gravity_longitudinal_minimum",
        "zero_gravity_longitudinal_minimum",
        "gravity_rotation_peak",
        "zero_gravity_rotation_peak",
        "gravity_periodicity_error",
        "zero_gravity_periodicity_error",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for gravity, zero in pairs:
            writer.writerow(
                {
                    "normalized_frequency": gravity["normalized_frequency"],
                    "gravity_frequency_hz": gravity["frequency_hz"],
                    "zero_gravity_frequency_hz": zero["frequency_hz"],
                    "gravity_transverse_peak": gravity["transverse_peak"],
                    "zero_gravity_transverse_peak": zero["transverse_peak"],
                    "gravity_longitudinal_minimum": gravity[
                        "longitudinal_minimum"
                    ],
                    "zero_gravity_longitudinal_minimum": zero[
                        "longitudinal_minimum"
                    ],
                    "gravity_rotation_peak": gravity["rotation_peak"],
                    "zero_gravity_rotation_peak": zero["rotation_peak"],
                    "gravity_periodicity_error": gravity[
                        "periodicity_error"
                    ],
                    "zero_gravity_periodicity_error": zero[
                        "periodicity_error"
                    ],
                }
            )


def plot_comparison(
    path: Path,
    gravity_rows: list[dict[str, str]],
    zero_gravity_rows: list[dict[str, str]],
) -> None:
    figure, axes = plt.subplots(
        1, 3, figsize=(10.5, 3.35), constrained_layout=True
    )
    quantities = (
        ("transverse_peak", r"$\max |w_{\mathrm{tip}}|/L$", "(a) Transverse"),
        (
            "longitudinal_minimum",
            r"$|\min u_{\mathrm{tip}}|/L$",
            "(b) Longitudinal",
        ),
        ("rotation_peak", r"$\max |\psi_{\mathrm{tip}}|$", "(c) Rotation"),
    )
    series = (
        (gravity_rows, "gravity", "#2864a5", "o"),
        (zero_gravity_rows, "zero gravity", "#c43c39", "D"),
    )
    for axis, (key, ylabel, title) in zip(axes, quantities):
        for rows, label, color, marker in series:
            accepted = [
                row
                for row in rows
                if float(row["periodicity_error"]) <= PERIODICITY_TOLERANCE
            ]
            rejected = [
                row
                for row in rows
                if float(row["periodicity_error"]) > PERIODICITY_TOLERANCE
            ]
            values = [float(row[key]) for row in accepted]
            if key == "longitudinal_minimum":
                values = [abs(value) for value in values]
            axis.plot(
                [float(row["normalized_frequency"]) for row in accepted],
                values,
                color=color,
                marker=marker,
                markersize=4,
                linewidth=1.35,
                label=label,
            )
            rejected_values = [float(row[key]) for row in rejected]
            if key == "longitudinal_minimum":
                rejected_values = [abs(value) for value in rejected_values]
            axis.scatter(
                [float(row["normalized_frequency"]) for row in rejected],
                rejected_values,
                color=color,
                marker="x",
                s=32,
                linewidths=1.2,
            )
        axis.set(xlabel=r"$f/f_1$", ylabel=ylabel, title=title)
        axis.grid(True, linewidth=0.4, alpha=0.3)
    axes[0].legend(frameon=False, fontsize=8)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=220)
    plt.close(figure)


def main() -> None:
    example_directory = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gravity",
        type=Path,
        default=(
            example_directory
            / "reference"
            / "farokhi_05g_k4n2_up_sweep.csv"
        ),
    )
    parser.add_argument(
        "--zero-gravity",
        type=Path,
        default=(
            example_directory
            / "results"
            / "gravity_sensitivity_05g_zero_k4n2.csv"
        ),
    )
    parser.add_argument(
        "--zero-gravity-metadata",
        type=Path,
        default=(
            example_directory
            / "results"
            / "gravity_sensitivity_05g_zero_k4n2_metadata.txt"
        ),
    )
    parser.add_argument(
        "--comparison-output",
        type=Path,
        default=(
            example_directory
            / "results"
            / "base_excited_cantilever_gravity_comparison.csv"
        ),
    )
    parser.add_argument(
        "--figure-output",
        type=Path,
        default=(
            example_directory
            / "results"
            / "base_excited_cantilever_gravity_comparison.png"
        ),
    )
    arguments = parser.parse_args()

    gravity_rows = read_rows(arguments.gravity)
    zero_gravity_rows = read_rows(arguments.zero_gravity)
    pairs = matching_rows(gravity_rows, zero_gravity_rows)
    gravity_matching = [pair[0] for pair in pairs]

    numerical_checks = (
        example_directory
        / "reference"
        / "base_excited_cantilever_numerical_checks.csv"
    )
    gravity_frequency = read_numerical_check(
        numerical_checks, "linear_frequency_hz_k4_n2"
    )
    zero_gravity_frequency = read_key_value(
        arguments.zero_gravity_metadata, "linear_frequency_hz"
    )
    gravity_last, gravity_transition = first_nonperiodic_bracket(
        gravity_matching
    )
    zero_last, zero_transition = first_nonperiodic_bracket(zero_gravity_rows)

    write_comparison(arguments.comparison_output, pairs)
    plot_comparison(arguments.figure_output, gravity_matching, zero_gravity_rows)

    linear_reduction = 1.0 - gravity_frequency / zero_gravity_frequency
    normalized_shift = float(gravity_last["normalized_frequency"]) - float(
        zero_last["normalized_frequency"]
    )
    physical_shift = float(gravity_last["frequency_hz"]) - float(
        zero_last["frequency_hz"]
    )
    print(
        "gravity lowers the discrete linear frequency by "
        f"{100.0 * linear_reduction:.3f}% "
        f"({zero_gravity_frequency:.6f} to {gravity_frequency:.6f} Hz)"
    )
    print(
        "last accepted upper point, zero/full gravity: "
        f"{float(zero_last['normalized_frequency']):.9f} / "
        f"{float(gravity_last['normalized_frequency']):.9f}"
    )
    print(
        "normalized last-upper shift from gravity: "
        f"{normalized_shift:+.9f}"
    )
    print(
        "physical-frequency last-upper shift from gravity: "
        f"{physical_shift:+.6f} Hz"
    )
    print(
        "first nonperiodic transition, zero/full gravity: "
        f"{float(zero_transition['normalized_frequency']):.9f} / "
        f"{float(gravity_transition['normalized_frequency']):.9f}"
    )
    print(f"wrote {arguments.comparison_output}")
    print(f"wrote {arguments.figure_output}")


if __name__ == "__main__":
    main()
