"""Figures: one chromatogram panel per read and one editing heatmap per guide."""

from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from fast_sanger_be.chromatogram import Chromatogram  # noqa: E402
from fast_sanger_be.config import BASE_COLORS, BASES, EDITABLE_BASE, Parameters  # noqa: E402

PERCENT_CMAP = LinearSegmentedColormap.from_list(
    "editing", ["#ffffff", "#c5d2e9", "#6a92c7", "#3e6bb6", "#2e3868"]
)
GUIDE_SHADE = "#d9d9d9"
MISSING_CELL = "#ececec"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "svg.fonttype": "none", "pdf.fonttype": 42})


def _text_color(value: float, vmax: float) -> str:
    return "white" if value > 0.6 * vmax else "black"


def plot_read(
    chrom: Chromatogram,
    start: int,
    guide: str,
    percentages: np.ndarray | None,
    title: str,
    subtitle: str,
    path: Path,
    params: Parameters,
    reliable: np.ndarray | None = None,
) -> None:
    """Chromatogram of the guide +/- qc_flank window (guide shaded) above a base-percentage grid.

    Guide positions that quality control could not quantify are drawn as empty grey cells.

    ``chrom`` must already be in guide orientation.
    """
    end = start + len(guide)
    lo, hi = max(0, start - params.qc_flank), min(len(chrom), end + params.qc_flank)
    locs = chrom.peak_locations
    spacing = float(np.median(np.diff(locs[lo:hi]))) if hi - lo > 1 else 12.0
    x0, x1 = locs[lo] - spacing / 2, locs[hi - 1] + spacing / 2
    s0, s1 = locs[start] - spacing / 2, locs[end - 1] + spacing / 2

    width = max(7.0, 0.19 * (hi - lo) + 1.5)
    fig, (ax_trace, ax_table) = plt.subplots(
        2, 1, figsize=(width, 4.4), gridspec_kw={"height_ratios": [2.3, 1.4], "hspace": 0.12}
    )
    ax_table.set_xlim(x0, x1)

    i0, i1 = int(max(0, np.floor(x0))), int(min(chrom.trace_length, np.ceil(x1) + 1))
    xs = np.arange(i0, i1)
    ymax = max(float(chrom.traces[b][i0:i1].max()) for b in BASES) * 1.08 or 1.0
    ax_trace.axvspan(s0, s1, color=GUIDE_SHADE, zorder=0, lw=0)
    for base in ("G", "A", "T", "C"):
        ax_trace.plot(xs, chrom.traces[base][i0:i1], color=BASE_COLORS[base], lw=0.8, label=base, zorder=2)
    ax_trace.set_ylim(0, ymax)
    ax_trace.set_xlim(x0, x1)
    ax_trace.set_yticks([])
    for side in ("top", "right", "left"):
        ax_trace.spines[side].set_visible(False)
    ax_trace.legend(loc="upper left", frameon=False, handlelength=1.2, fontsize=7, bbox_to_anchor=(1.005, 1))
    wrapped = "\n".join(textwrap.wrap(subtitle, width=max(60, int(width * 11))))
    ax_trace.set_title(f"{title}\n{wrapped}", loc="left", fontsize=8.5)
    ax_trace.text(
        (s0 + s1) / 2, ymax * 0.985, "guide", ha="center", va="top", fontsize=7, color="#555555"
    )

    # Called bases along the window, guide bases in bold.
    ax_trace.set_xticks(locs[lo:hi])
    ax_trace.set_xticklabels(list(chrom.bases[lo:hi]), fontsize=7, family="DejaVu Sans Mono")
    for idx, label in zip(range(lo, hi), ax_trace.get_xticklabels()):
        label.set_color(BASE_COLORS.get(label.get_text(), "#777777"))
        if start <= idx < end:
            label.set_fontweight("bold")
    ax_trace.tick_params(axis="x", length=0, pad=2, labelbottom=True)

    # Percentage grid aligned with the guide peaks.
    ax_table.set_ylim(4, 0)
    ax_table.set_yticks(np.arange(4) + 0.5)
    ax_table.set_yticklabels([f"%{b}" for b in BASES], fontsize=7)
    for tick, base in zip(ax_table.get_yticklabels(), BASES):
        tick.set_color(BASE_COLORS[base])
        tick.set_fontweight("bold")
    ax_table.tick_params(length=0)
    for side in ax_table.spines.values():
        side.set_visible(False)
    ax_table.set_xticks([locs[i] for i in range(start, end)])
    ax_table.set_xticklabels([f"{guide[i - start]}{i - start + 1}" for i in range(start, end)], fontsize=6.5, rotation=90)
    for i in range(start, end):
        cx = locs[i]
        left = (locs[i - 1] + cx) / 2 if i > 0 else cx - spacing / 2
        right = (cx + locs[i + 1]) / 2 if i + 1 < len(locs) else cx + spacing / 2
        usable = reliable is None or bool(reliable[i - start])
        for row, base in enumerate(BASES):
            value = np.nan if percentages is None or not usable else percentages[i - start, row]
            face = MISSING_CELL if np.isnan(value) else PERCENT_CMAP(value / 100)
            ax_table.add_patch(
                Rectangle((left, row), right - left, 1, facecolor=face, edgecolor="white", lw=0.8)
            )
            if not np.isnan(value):
                ax_table.text(
                    cx, row + 0.5, f"{value:.0f}", ha="center", va="center", fontsize=5,
                    color=_text_color(value, 100),
                )

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=params.figure_dpi, bbox_inches="tight")
    plt.close(fig)


def plot_heatmap(
    guide_name: str,
    guide: str,
    samples: list[str],
    values: dict[str, dict[str, float]],
    editor: str,
    path: Path,
    params: Parameters,
) -> None:
    """Samples x guide positions; target positions coloured by % editing, others grey."""
    labels = [f"{b}{i + 1}" for i, b in enumerate(guide)]
    matrix = np.full((len(samples), len(guide)), np.nan)
    for r, sample in enumerate(samples):
        for c, label in enumerate(labels):
            value = values.get(sample, {}).get(label)
            if value is not None and not np.isnan(value):
                matrix[r, c] = value

    vmax = params.heatmap_max
    cmap = PERCENT_CMAP.with_extremes(bad=MISSING_CELL)
    fig_w = max(6.0, 0.36 * len(guide) + 2.6)
    fig_h = max(2.2, 0.32 * len(samples) + 1.4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    image = ax.imshow(np.ma.masked_invalid(matrix), cmap=cmap, vmin=0, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(guide)), labels, rotation=90, fontsize=7, family="DejaVu Sans Mono")
    row_labels = [s if values.get(s) else f"{s} (no valid read)" for s in samples]
    ax.set_yticks(range(len(samples)), row_labels, fontsize=7.5)
    ax.set_xticks(np.arange(-0.5, len(guide)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(samples)), minor=True)
    ax.grid(which="minor", color="white", lw=1.2)
    ax.tick_params(which="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    editable = EDITABLE_BASE[editor]
    hatched = False
    for r, sample in enumerate(samples):
        for c in range(len(guide)):
            value = matrix[r, c]
            if not np.isnan(value):
                ax.text(c, r, f"{value:.1f}", ha="center", va="center", fontsize=6,
                        color=_text_color(value, vmax))
            elif guide[c] == editable and values.get(sample):
                # A target position QC could not quantify: hatched, to tell it apart from non-targets.
                ax.add_patch(Rectangle((c - 0.5, r - 0.5), 1, 1, facecolor="none", edgecolor="#9a9a9a",
                                       hatch="////", lw=0))
                hatched = True
    bar = fig.colorbar(image, ax=ax, fraction=0.03, pad=0.02)
    bar.set_label("Editing (%)", fontsize=7.5)
    bar.ax.tick_params(labelsize=7)
    conversion = "A>G" if editor == "ABE" else "C>T"
    ax.set_title(f"{guide_name} ({editor}, {conversion})  5'-{guide}-3'", loc="left", fontsize=8.5)
    key = "Grey: not a target base."
    if hatched:
        key += "  Hatched: target not quantified (overlapping peaks nearby)."
    ax.annotate(key,
                xy=(0, 0), xycoords="axes fraction", xytext=(0, -34), textcoords="offset points",
                fontsize=6.5, color="#555555", va="top")

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=params.figure_dpi, bbox_inches="tight")
    plt.close(fig)
