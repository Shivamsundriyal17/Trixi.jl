#!/usr/bin/env python3
"""Plot and animate one refined extreme-cantilever cycle."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation, colors


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    example_directory = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-directory",
        type=Path,
        default=(
            example_directory
            / "results"
            / "base_excited_cantilever_extreme_cycle"
        ),
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=example_directory / "results",
    )
    parser.add_argument("--fps", type=int, default=24)
    arguments = parser.parse_args()

    geometry_rows = read_rows(
        arguments.input_directory / "centerline_cycle.csv"
    )
    tip_rows = read_rows(arguments.input_directory / "tip_cycle.csv")
    ledger = read_rows(arguments.input_directory / "cycle_ledger.csv")[0]

    frames: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in geometry_rows:
        frames[int(row["frame"])].append(row)
    frame_indices = sorted(frames)
    phases = np.array([float(row["phase"]) for row in tip_rows])
    longitudinal = np.array(
        [float(row["longitudinal"]) for row in tip_rows]
    )
    transverse = np.array(
        [float(row["transverse"]) for row in tip_rows]
    )
    rotation = np.array([float(row["rotation"]) for row in tip_rows])

    def frame_coordinates(index: int) -> tuple[np.ndarray, np.ndarray]:
        rows = frames[index]
        return (
            np.array([float(row["transverse"]) for row in rows]),
            np.array([float(row["vertical"]) for row in rows]),
        )

    all_transverse = np.concatenate(
        [frame_coordinates(index)[0] for index in frame_indices]
    )
    all_vertical = np.concatenate(
        [frame_coordinates(index)[1] for index in frame_indices]
    )
    x_padding = 0.08 * max(
        np.ptp(all_transverse), np.ptp(all_vertical), 1.0
    )
    y_padding = x_padding
    x_limits = (
        float(np.min(all_transverse) - x_padding),
        float(np.max(all_transverse) + x_padding),
    )
    y_limits = (
        float(np.min(all_vertical) - y_padding),
        float(np.max(all_vertical) + y_padding),
    )

    arguments.output_directory.mkdir(parents=True, exist_ok=True)
    static_png = (
        arguments.output_directory
        / "base_excited_cantilever_extreme_cycle.png"
    )
    static_pdf = static_png.with_suffix(".pdf")
    gif_path = static_png.with_suffix(".gif")

    figure, axes = plt.subplots(
        1, 3, figsize=(11.2, 3.8), constrained_layout=True
    )
    snapshot_indices = np.unique(
        np.round(
            np.linspace(0, len(frame_indices) - 2, 9)
        ).astype(int)
    )
    colormap = plt.get_cmap("viridis")
    normalization = colors.Normalize(vmin=0.0, vmax=1.0)
    for index in snapshot_indices:
        x_values, y_values = frame_coordinates(frame_indices[index])
        axes[0].plot(
            x_values,
            y_values,
            color=colormap(phases[index]),
            linewidth=1.6,
        )
    axes[0].plot(
        [0.0, 0.0],
        [0.0, 1.0],
        color="0.65",
        linestyle=":",
        linewidth=1.0,
        label="undeformed",
    )
    axes[0].scatter([0.0], [0.0], marker="s", color="black", s=20, zorder=5)
    axes[0].set(
        xlabel=r"$w/L$",
        ylabel=r"vertical coordinate$/L$",
        xlim=x_limits,
        ylim=y_limits,
        aspect="equal",
        title="(a) Refined centerline snapshots",
    )
    axes[0].grid(True, linewidth=0.4, alpha=0.3)
    scalar_mappable = plt.cm.ScalarMappable(
        norm=normalization, cmap=colormap
    )
    colorbar = figure.colorbar(
        scalar_mappable, ax=axes[0], fraction=0.045, pad=0.03
    )
    colorbar.set_label("cycle phase")

    axes[1].plot(
        transverse,
        1.0 + longitudinal,
        color="#2864a5",
        linewidth=1.7,
    )
    phase_markers = np.unique(
        np.round(
            np.linspace(0, len(phases) - 2, 9)
        ).astype(int)
    )
    axes[1].scatter(
        transverse[phase_markers],
        1.0 + longitudinal[phase_markers],
        c=phases[phase_markers],
        cmap=colormap,
        norm=normalization,
        s=28,
        edgecolors="white",
        linewidths=0.5,
        zorder=3,
    )
    axes[1].set(
        xlabel=r"$w_{\mathrm{tip}}/L$",
        ylabel=r"$(L+u_{\mathrm{tip}})/L$",
        title="(b) Tip orbit",
        aspect="equal",
    )
    axes[1].grid(True, linewidth=0.4, alpha=0.3)

    axes[2].plot(
        phases,
        transverse,
        color="#2864a5",
        linewidth=1.5,
        label=r"$w_{\mathrm{tip}}/L$",
    )
    axes[2].plot(
        phases,
        longitudinal,
        color="#c43c39",
        linewidth=1.5,
        label=r"$u_{\mathrm{tip}}/L$",
    )
    axes[2].plot(
        phases,
        rotation / np.pi,
        color="#2f8f5b",
        linewidth=1.5,
        label=r"$\psi_{\mathrm{tip}}/\pi$",
    )
    axes[2].set(
        xlabel="cycle phase",
        ylabel="normalized response",
        xlim=(0.0, 1.0),
        title="(c) Periodic tip response",
    )
    axes[2].grid(True, linewidth=0.4, alpha=0.3)
    axes[2].legend(frameon=False, fontsize=8, loc="lower left")

    normalized_frequency = float(ledger["normalized_frequency"])
    figure.suptitle(
        rf"Extreme base-excited cantilever: $0.5g$, "
        rf"$f/f_1={normalized_frequency:.5f}$, "
        rf"$k=3$, four cells",
        fontsize=11,
    )
    figure.savefig(static_png, dpi=240)
    figure.savefig(static_pdf)
    plt.close(figure)

    animation_figure, axis = plt.subplots(figsize=(5.2, 5.2))
    axis.plot(
        [0.0, 0.0],
        [0.0, 1.0],
        color="0.75",
        linestyle=":",
        linewidth=1.0,
    )
    axis.plot(
        transverse,
        1.0 + longitudinal,
        color="0.75",
        linestyle="--",
        linewidth=0.9,
        label="tip orbit",
    )
    beam_line, = axis.plot([], [], color="#2864a5", linewidth=3.0)
    tip_marker, = axis.plot(
        [], [], marker="o", color="#c43c39", markersize=6
    )
    phase_text = axis.text(
        0.03, 0.96, "", transform=axis.transAxes, va="top"
    )
    axis.scatter([0.0], [0.0], marker="s", color="black", s=28, zorder=5)
    axis.set(
        xlabel=r"$w/L$",
        ylabel=r"vertical coordinate$/L$",
        xlim=x_limits,
        ylim=y_limits,
        aspect="equal",
        title=rf"Refined extreme orbit: $0.5g$, $f/f_1={normalized_frequency:.5f}$",
    )
    axis.grid(True, linewidth=0.4, alpha=0.3)
    axis.legend(frameon=False, fontsize=8, loc="upper right")

    animation_frames = frame_indices[:-1]

    def update(frame_number: int):
        frame_index = animation_frames[frame_number]
        x_values, y_values = frame_coordinates(frame_index)
        beam_line.set_data(x_values, y_values)
        tip_marker.set_data([x_values[-1]], [y_values[-1]])
        phase_text.set_text(f"cycle phase = {phases[frame_number]:.3f}")
        return beam_line, tip_marker, phase_text

    movie = animation.FuncAnimation(
        animation_figure,
        update,
        frames=len(animation_frames),
        interval=1000 / arguments.fps,
        blit=True,
    )
    movie.save(
        gif_path,
        writer=animation.PillowWriter(fps=arguments.fps),
        dpi=120,
    )
    plt.close(animation_figure)

    print(f"Wrote {static_png}")
    print(f"Wrote {static_pdf}")
    print(f"Wrote {gif_path}")


if __name__ == "__main__":
    main()
