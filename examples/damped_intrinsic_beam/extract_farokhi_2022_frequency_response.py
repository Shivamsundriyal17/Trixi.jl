#!/usr/bin/env python3
"""Digitize the experimental markers in Farokhi et al. (2022).

The publisher PDF stores the red circular markers as vector paths. This script
extracts their centers rather than sampling raster pixels. It intentionally
does not redistribute the article: pass a legally obtained PDF as input.

Reference:
  H. Farokhi, Y. Xia, and A. Erturk,
  Nonlinear Dynamics 107 (2022), 457--475.
  https://doi.org/10.1007/s11071-021-07023-9
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path


MARKER_COLOR = "#ff1d09"
# The marker path spans x in [-2.374, 0], so its center is at -1.187.
MARKER_CENTER_OFFSET_X = -1.187


@dataclass(frozen=True)
class Axes:
    x_svg_min: float
    x_svg_max: float
    y_svg_at_data_min: float
    y_svg_at_data_max: float
    x_data_min: float
    x_data_max: float
    y_data_min: float
    y_data_max: float

    def transform(self, point: tuple[float, float]) -> tuple[float, float]:
        x_svg, y_svg = point
        x_data = self.x_data_min + (
            (x_svg - self.x_svg_min)
            / (self.x_svg_max - self.x_svg_min)
            * (self.x_data_max - self.x_data_min)
        )
        y_data = self.y_data_min + (
            (y_svg - self.y_svg_at_data_min)
            / (self.y_svg_at_data_max - self.y_svg_at_data_min)
            * (self.y_data_max - self.y_data_min)
        )
        return x_data, y_data


FIGURE_4_TRANSVERSE = Axes(
    71.378, 262.454, 477.300, 295.323, 0.8, 1.2, 0.0, 0.9
)
FIGURE_4_LONGITUDINAL = Axes(
    303.028, 491.598, 296.594, 476.185, 0.8, 1.2, 0.0, -0.7
)
FIGURE_7_TRANSVERSE = Axes(
    71.378, 262.454, 250.777, 68.800, 0.8, 1.2, 0.0, 1.0
)
FIGURE_7_LONGITUDINAL = Axes(
    303.028, 491.598, 70.072, 249.663, 0.8, 1.2, 0.0, -1.2
)


def render_svg(pdf_path: Path, page: int, output_directory: Path) -> Path:
    output_pattern = output_directory / "page-%d.svg"
    subprocess.run(
        [
            "mutool",
            "draw",
            "-F",
            "svg",
            "-o",
            str(output_pattern),
            str(pdf_path),
            str(page),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return output_directory / f"page-{page}.svg"


def marker_centers(svg_path: Path) -> list[tuple[float, float]]:
    svg = svg_path.read_text(encoding="utf-8")
    pattern = re.compile(
        r'<path transform="matrix\(1,0,0,-1,([^,]+),([^)]+)\)"'
        rf'[^>]*stroke="{MARKER_COLOR}"'
    )
    return [
        (float(match.group(1)) + MARKER_CENTER_OFFSET_X,
         float(match.group(2)))
        for match in pattern.finditer(svg)
    ]


def paired_rows(
    transverse_points: list[tuple[float, float]],
    longitudinal_points: list[tuple[float, float]],
    transverse_axes: Axes,
    longitudinal_axes: Axes,
    acceleration_rms_g: float,
    figure: int,
) -> list[dict[str, float | int]]:
    if len(transverse_points) != len(longitudinal_points):
        raise ValueError("transverse and longitudinal marker counts differ")

    rows = []
    for transverse_point, longitudinal_point in zip(
        transverse_points, longitudinal_points, strict=True
    ):
        x_transverse, transverse = transverse_axes.transform(transverse_point)
        x_longitudinal, longitudinal = longitudinal_axes.transform(
            longitudinal_point
        )
        rows.append(
            {
                "acceleration_rms_g": acceleration_rms_g,
                "normalized_frequency": 0.5
                * (x_transverse + x_longitudinal),
                "transverse_peak": transverse,
                "longitudinal_minimum": longitudinal,
                "figure": figure,
            }
        )
    return rows


def extract(pdf_path: Path) -> list[dict[str, float | int]]:
    if shutil.which("mutool") is None:
        raise RuntimeError("mutool is required (provided by MuPDF)")

    with tempfile.TemporaryDirectory(prefix="farokhi-2022-") as directory:
        temporary_directory = Path(directory)
        page_8 = marker_centers(render_svg(pdf_path, 8, temporary_directory))
        page_11 = marker_centers(render_svg(pdf_path, 11, temporary_directory))

    # Page 8 contains 33 red centerline markers from Figure 3, followed by
    # 30 transverse and 30 longitudinal frequency-response markers.
    if len(page_8) != 93:
        raise ValueError(f"expected 93 red paths on page 8, found {len(page_8)}")
    figure_4 = paired_rows(
        page_8[33:63],
        page_8[63:93],
        FIGURE_4_TRANSVERSE,
        FIGURE_4_LONGITUDINAL,
        0.2,
        4,
    )

    # The first 68 red paths on page 11 are the 34+34 markers of Figure 7.
    if len(page_11) < 68:
        raise ValueError(
            f"expected at least 68 red paths on page 11, found {len(page_11)}"
        )
    figure_7 = paired_rows(
        page_11[0:34],
        page_11[34:68],
        FIGURE_7_TRANSVERSE,
        FIGURE_7_LONGITUDINAL,
        0.5,
        7,
    )
    return figure_4 + figure_7


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("pdf", type=Path, help="Farokhi et al. article PDF")
    parser.add_argument("output", type=Path, help="output CSV path")
    arguments = parser.parse_args()

    rows = extract(arguments.pdf.resolve())
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    with arguments.output.open("w", newline="", encoding="utf-8") as stream:
        fieldnames = [
            "acceleration_rms_g",
            "normalized_frequency",
            "transverse_peak",
            "longitudinal_minimum",
            "figure",
        ]
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} markers to {arguments.output}")


if __name__ == "__main__":
    main()
