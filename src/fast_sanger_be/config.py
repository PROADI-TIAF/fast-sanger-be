"""Analysis parameters and their documentation."""

from __future__ import annotations

from dataclasses import dataclass, fields

EDITORS = ("ABE", "CBE")
READS = ("F", "R")

# Base targeted by each editor and the base it is converted to (guide strand).
EDITABLE_BASE = {"ABE": "A", "CBE": "C"}
PRODUCT_BASE = {"ABE": "G", "CBE": "T"}

BASES = "ACGT"
BASE_COLORS = {"A": "#1a9850", "C": "#2c61c5", "G": "#000000", "T": "#d7301f"}


@dataclass(frozen=True)
class Parameters:
    """All tunable settings of a run. Every value is written to the output workbooks."""

    qc_flank: int = 10
    background_flank: int = 50
    peak_half_width: int = 2
    max_guide_mismatches: int = 3
    min_background_phred: float = 20.0
    min_peak_height: float = 100.0
    min_background_positions: int = 20
    mixed_peak_ratio: float = 0.25
    position_flank: int = 5
    max_local_mixed_fraction: float = 0.35
    max_window_mixed_fraction: float = 0.60
    partial_reads: bool = False
    mad_cutoff: float = 3.0
    discordance_threshold: float = 10.0
    heatmap_max: float = 100.0
    figure_format: str = "png"
    figure_dpi: int = 300

    def rows(self) -> list[tuple[str, object, str]]:
        """Return (name, value, description) for every parameter."""
        return [(f.name, getattr(self, f.name), PARAMETER_HELP[f.name]) for f in fields(self)]


PARAMETER_HELP = {
    "qc_flank": "Nucleotides on each side of the guide included in the quality-control window.",
    "background_flank": "Nucleotides on each side of the guide used to estimate background noise "
    "(the guide itself is excluded).",
    "peak_half_width": "Trace samples on each side of a peak location searched for the peak maximum.",
    "max_guide_mismatches": "Maximum mismatches allowed when locating the guide in the base calls "
    "(expected edits are not counted as mismatches).",
    "min_background_phred": "Minimum median Phred quality of the region flanking the guide, which is where "
    "background noise is measured.",
    "min_peak_height": "Minimum median primary peak height (RFU) of the region flanking the guide.",
    "min_background_positions": "Minimum number of usable positions for the background estimate.",
    "mixed_peak_ratio": "Secondary/primary peak height ratio above which a position counts as an overlapping peak.",
    "position_flank": "Nucleotides on each side of a guide position whose neighbours decide whether that "
    "position can be quantified.",
    "max_local_mixed_fraction": "A target position is noisy when more than this fraction of its neighbours "
    "(positions the editor cannot convert) have overlapping peaks. One noisy target rejects the whole read "
    "unless partial_reads is set.",
    "max_window_mixed_fraction": "A read is rejected when more than this fraction of the whole QC window, "
    "excluding the positions the editor can convert, has overlapping peaks.",
    "partial_reads": "If true, a read with noisy target positions is kept and only those positions are left "
    "empty. If false (default), any noisy target position rejects the whole read.",
    "mad_cutoff": "Background values further than this many scaled MADs from the median are discarded.",
    "discordance_threshold": "Forward/reverse difference (percentage points) that is flagged in the final results.",
    "heatmap_max": "Upper limit of the heatmap colour scale (% editing).",
    "figure_format": "File format of the figures (png, pdf, svg or tiff).",
    "figure_dpi": "Resolution of raster figures.",
}
