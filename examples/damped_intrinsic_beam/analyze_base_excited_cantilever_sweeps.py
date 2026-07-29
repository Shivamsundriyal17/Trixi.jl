#!/usr/bin/env python3
"""Compare converged Trixi cantilever branches with Farokhi et al. markers."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Campaign:
    label: str
    branch: str
    filename: str
    maximum_normalized_frequency: float | None = None
    maximum_periodicity_error: float | None = None


CAMPAIGNS = (
    Campaign(
        "02g_upper",
        "upper",
        "farokhi_02g_k4n2_sweep.csv",
        maximum_normalized_frequency=1.0206086,
        maximum_periodicity_error=0.002,
    ),
    Campaign(
        "02g_lower",
        "lower",
        "farokhi_02g_k4n2_down_sweep.csv",
        maximum_periodicity_error=0.0002,
    ),
    Campaign(
        "05g_upper",
        "upper",
        "farokhi_05g_k4n2_up_sweep.csv",
        maximum_normalized_frequency=1.0411751,
        maximum_periodicity_error=0.0005,
    ),
    Campaign(
        "05g_lower",
        "lower",
        "farokhi_05g_k4n2_down_sweep.csv",
        maximum_periodicity_error=0.0002,
    ),
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def experimental_marker(
    rows: list[dict[str, str]],
    acceleration: float,
    normalized_frequency: float,
    branch: str,
) -> dict[str, str]:
    candidates = [
        row
        for row in rows
        if abs(float(row["acceleration_rms_g"]) - acceleration) < 1.0e-12
        and abs(float(row["normalized_frequency"]) - normalized_frequency)
        < 2.0e-6
    ]
    if not candidates:
        raise ValueError(
            f"no experimental marker at {acceleration}g, "
            f"frequency {normalized_frequency}"
        )
    selector = max if branch == "upper" else min
    return selector(candidates, key=lambda row: float(row["transverse_peak"]))


def error_metrics(
    comparisons: list[dict[str, float | str]], numerical_key: str,
    experimental_key: str
) -> dict[str, float]:
    errors = [
        float(row[numerical_key]) - float(row[experimental_key])
        for row in comparisons
    ]
    experimental_peak = max(
        abs(float(row[experimental_key])) for row in comparisons
    )
    rmse = math.sqrt(sum(error * error for error in errors) / len(errors))
    return {
        "rmse": rmse,
        "peak_normalized_rmse": rmse / experimental_peak,
        "mae": sum(abs(error) for error in errors) / len(errors),
        "bias": sum(errors) / len(errors),
        "maximum_absolute_error": max(abs(error) for error in errors),
    }


def analyze_campaign(
    campaign: Campaign,
    results_directory: Path,
    experimental_rows: list[dict[str, str]],
) -> tuple[
    list[dict[str, float | str]],
    list[dict[str, float | str]],
    dict[str, float | str],
]:
    numerical_rows = read_rows(results_directory / campaign.filename)
    comparisons: list[dict[str, float | str]] = []
    for numerical in numerical_rows:
        frequency = float(numerical["normalized_frequency"])
        periodicity_error = float(numerical["periodicity_error"])
        if (
            campaign.maximum_normalized_frequency is not None
            and frequency > campaign.maximum_normalized_frequency
        ):
            continue
        if (
            campaign.maximum_periodicity_error is not None
            and periodicity_error > campaign.maximum_periodicity_error
        ):
            continue

        acceleration = float(numerical["acceleration_rms_g"])
        experimental = experimental_marker(
            experimental_rows, acceleration, frequency, campaign.branch
        )
        numerical_transverse = float(numerical["transverse_peak"])
        experimental_transverse = float(experimental["transverse_peak"])
        numerical_longitudinal = float(numerical["longitudinal_minimum"])
        experimental_longitudinal = float(
            experimental["longitudinal_minimum"]
        )
        physical_root_work = float(numerical["physical_root_work"])
        left_boundary_dissipation = float(
            numerical["left_boundary_dissipation"]
        )
        right_boundary_dissipation = float(
            numerical["right_boundary_dissipation"]
        )
        sat_data_work = float(numerical["sat_data_work"])
        net_boundary_dissipation = (
            left_boundary_dissipation
            + right_boundary_dissipation
            - sat_data_work
        )
        comparisons.append(
            {
                "campaign": campaign.label,
                "branch": campaign.branch,
                "acceleration_rms_g": acceleration,
                "normalized_frequency": frequency,
                "frequency_hz": float(numerical["frequency_hz"]),
                "numerical_transverse": numerical_transverse,
                "experimental_transverse": experimental_transverse,
                "transverse_error": (
                    numerical_transverse - experimental_transverse
                ),
                "numerical_longitudinal": numerical_longitudinal,
                "experimental_longitudinal": experimental_longitudinal,
                "longitudinal_error": (
                    numerical_longitudinal - experimental_longitudinal
                ),
                "periodicity_error": periodicity_error,
                "total_energy_change": float(
                    numerical["total_energy_change"]
                ),
                "physical_root_work": physical_root_work,
                "material_dissipation": float(
                    numerical["material_dissipation"]
                ),
                "jump_dissipation": float(numerical["jump_dissipation"]),
                "left_boundary_dissipation": left_boundary_dissipation,
                "right_boundary_dissipation": right_boundary_dissipation,
                "sat_data_work": sat_data_work,
                "net_boundary_dissipation": net_boundary_dissipation,
                "relative_ledger_residual": float(
                    numerical["relative_ledger_residual"]
                ),
            }
        )

    if not comparisons:
        raise ValueError(f"campaign {campaign.label} has no accepted rows")

    summary: list[dict[str, float | str]] = []
    for component, numerical_key, experimental_key in (
        (
            "transverse",
            "numerical_transverse",
            "experimental_transverse",
        ),
        (
            "longitudinal",
            "numerical_longitudinal",
            "experimental_longitudinal",
        ),
    ):
        metrics = error_metrics(comparisons, numerical_key, experimental_key)
        summary.append(
            {
                "campaign": campaign.label,
                "branch": campaign.branch,
                "component": component,
                "points": len(comparisons),
                **metrics,
                "maximum_periodicity_error": max(
                    float(row["periodicity_error"]) for row in comparisons
                ),
                "maximum_relative_ledger_residual": max(
                    float(row["relative_ledger_residual"])
                    for row in comparisons
                ),
            }
        )

    peak_row = max(
        comparisons, key=lambda row: float(row["numerical_transverse"])
    )
    root_work = float(peak_row["physical_root_work"])
    if root_work <= 0.0:
        raise ValueError(
            f"campaign {campaign.label} has non-positive root work "
            "at its maximum accepted response"
        )
    compact_ledger_residual = (
        float(peak_row["total_energy_change"])
        + float(peak_row["material_dissipation"])
        + float(peak_row["jump_dissipation"])
        + float(peak_row["net_boundary_dissipation"])
        - root_work
    )
    energy_summary: dict[str, float | str] = {
        "campaign": campaign.label,
        "branch": campaign.branch,
        "normalized_frequency": float(peak_row["normalized_frequency"]),
        "numerical_transverse": float(peak_row["numerical_transverse"]),
        "periodicity_error": float(peak_row["periodicity_error"]),
        "total_energy_change": float(peak_row["total_energy_change"]),
        "physical_root_work": root_work,
        "material_dissipation": float(peak_row["material_dissipation"]),
        "jump_dissipation": float(peak_row["jump_dissipation"]),
        "net_boundary_dissipation": float(
            peak_row["net_boundary_dissipation"]
        ),
        "energy_change_fraction": (
            float(peak_row["total_energy_change"]) / root_work
        ),
        "material_fraction": (
            float(peak_row["material_dissipation"]) / root_work
        ),
        "jump_fraction": (
            float(peak_row["jump_dissipation"]) / root_work
        ),
        "net_boundary_fraction": (
            float(peak_row["net_boundary_dissipation"]) / root_work
        ),
        "compact_ledger_residual": compact_ledger_residual,
        "compact_relative_ledger_residual": (
            abs(compact_ledger_residual) / root_work
        ),
        "relative_ledger_residual": float(
            peak_row["relative_ledger_residual"]
        ),
    }
    return comparisons, summary, energy_summary


def write_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


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
        "--results-directory",
        type=Path,
        default=example_directory / "results",
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=example_directory / "results",
    )
    arguments = parser.parse_args()

    experimental_rows = read_rows(arguments.experimental)
    comparisons: list[dict[str, float | str]] = []
    summaries: list[dict[str, float | str]] = []
    energy_summaries: list[dict[str, float | str]] = []
    for campaign in CAMPAIGNS:
        (
            campaign_comparisons,
            campaign_summary,
            campaign_energy_summary,
        ) = analyze_campaign(
            campaign, arguments.results_directory, experimental_rows
        )
        comparisons.extend(campaign_comparisons)
        summaries.extend(campaign_summary)
        energy_summaries.append(campaign_energy_summary)

    comparison_path = (
        arguments.output_directory
        / "base_excited_cantilever_branch_comparison.csv"
    )
    summary_path = (
        arguments.output_directory
        / "base_excited_cantilever_branch_summary.csv"
    )
    energy_summary_path = (
        arguments.output_directory
        / "base_excited_cantilever_energy_summary.csv"
    )
    write_csv(comparison_path, comparisons)
    write_csv(summary_path, summaries)
    write_csv(energy_summary_path, energy_summaries)
    print(f"Wrote {comparison_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {energy_summary_path}")


if __name__ == "__main__":
    main()
