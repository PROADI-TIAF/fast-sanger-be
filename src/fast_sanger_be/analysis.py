"""Guide location, quality control and quantification of base editing.

The quantification follows the background-subtraction principle introduced by BEAT
(Xu et al., CRISPR J 2019): noise is estimated per dye channel from positions where
that channel is not the primary peak, outliers are removed with the median absolute
deviation (MAD), and the noise is subtracted before computing base percentages. This
implementation differs in that background comes from the regions flanking the guide
(never from the guide itself), signal and background use the same peak window, the
dye order is read from the file, and the guide search tolerates the expected edits.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from Bio.Seq import reverse_complement

from fast_sanger_be.chromatogram import Chromatogram
from fast_sanger_be.config import BASES, EDITABLE_BASE, PRODUCT_BASE, Parameters

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
}


# ---------------------------------------------------------------- guide location


@dataclass(frozen=True)
class GuideHit:
    start: int  # 0-based index of the first guide base in the oriented read
    strand: str  # '+' guide on the read as sequenced, '-' guide on its reverse complement
    mismatches: int
    ambiguous: bool


def _allowed_bases(guide: str, editor: str) -> list[set[str]]:
    editable, product = EDITABLE_BASE[editor], PRODUCT_BASE[editor]
    return [{b, product} if b == editable else {b} for b in guide]


def _scan(read: str, allowed: list[set[str]], max_mismatches: int) -> list[tuple[int, int]]:
    hits = []
    width = len(allowed)
    for start in range(len(read) - width + 1):
        mismatches = 0
        for offset, ok in enumerate(allowed):
            # A mixed call (R, Y, K, ...) matches only where every base it stands for is expected,
            # e.g. R = A/G at an ABE target; elsewhere it counts as a mismatch.
            called = IUPAC.get(read[start + offset])
            if not called or not ok.issuperset(called):
                mismatches += 1
                if mismatches > max_mismatches:
                    break
        else:
            hits.append((mismatches, start))
    return hits


def locate_guide(bases: str, guide: str, editor: str, max_mismatches: int) -> GuideHit | None:
    """Find the guide on either strand of the base calls, treating expected edits as matches."""
    allowed = _allowed_bases(guide, editor)
    hits = [(mm, start, "+") for mm, start in _scan(bases, allowed, max_mismatches)]
    hits += [(mm, start, "-") for mm, start in _scan(reverse_complement(bases), allowed, max_mismatches)]
    if not hits:
        return None
    hits.sort()
    mismatches, start, strand = hits[0]
    # Another site is a real alternative only if it matches about as well as the best one.
    ambiguous = len(hits) > 1 and hits[1][0] - mismatches <= 1
    return GuideHit(start, strand, mismatches, ambiguous)


def orient(chromatogram: Chromatogram, hit: GuideHit) -> Chromatogram:
    """Return the chromatogram in guide orientation (guide read 5'->3' left to right)."""
    return chromatogram if hit.strand == "+" else chromatogram.reverse_complement()


# ---------------------------------------------------------------- quality control


@dataclass
class QCResult:
    passed: bool
    reasons: list[str] = field(default_factory=list)
    background_median_phred: float | None = None
    background_median_peak_height: float | None = None
    flank_mixed_fraction: float | None = None
    window_mixed_fraction: float | None = None
    reliable: np.ndarray | None = None  # one flag per guide position
    local_mixed: np.ndarray | None = None  # fraction of mixed neighbours per guide position

    def unreliable_labels(self, guide: str, editor: str) -> list[str]:
        if self.reliable is None:
            return []
        return [position_label(guide, i) for i in target_positions(guide, editor) if not self.reliable[i]]


def secondary_ratio(heights: np.ndarray) -> np.ndarray:
    """Secondary/primary peak height ratio for each row of a (n, 4) height matrix."""
    ordered = np.sort(heights, axis=1)
    primary, secondary = ordered[:, -1], ordered[:, -2]
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(primary > 0, secondary / primary, 1.0)


def position_reliability(
    chrom: Chromatogram, start: int, guide: str, editor: str, params: Parameters
) -> tuple[np.ndarray, np.ndarray]:
    """Decide, position by position, whether the chromatogram can be quantified there.

    A position is judged by its neighbours rather than by itself, because a genuine edit
    produces two peaks at the edited position: only positions the editor cannot convert
    carry information about noise. Returns (reliable, fraction of mixed neighbours).
    """
    mixed = secondary_ratio(chrom.peak_heights(params.peak_half_width)) > params.mixed_peak_ratio
    editable = np.zeros(len(chrom), dtype=bool)
    editable[start : start + len(guide)] = [b == EDITABLE_BASE[editor] for b in guide]

    reliable = np.ones(len(guide), dtype=bool)
    local = np.zeros(len(guide))
    for i in range(len(guide)):
        position = start + i
        lo = max(0, position - params.position_flank)
        hi = min(len(chrom), position + params.position_flank + 1)
        neighbours = np.arange(lo, hi)
        neighbours = neighbours[(neighbours != position) & ~editable[neighbours]]
        if neighbours.size == 0:
            continue
        local[i] = float(mixed[neighbours].mean())
        reliable[i] = local[i] <= params.max_local_mixed_fraction
    return reliable, local


def assess_quality(
    chrom: Chromatogram, start: int, guide: str, editor: str, params: Parameters
) -> QCResult:
    """Quality of the guide +/- qc_flank window of an oriented chromatogram.

    Overlapping peaks around any target position reject the whole read. With
    ``params.partial_reads`` the read is kept instead, and only the affected target
    positions are left unquantified.
    """
    end = start + len(guide)
    lo, hi = start - params.qc_flank, end + params.qc_flank
    if lo < 0 or hi > len(chrom):
        return QCResult(False, [f"read does not cover the guide +/- {params.qc_flank} nt"])

    all_heights = chrom.peak_heights(params.peak_half_width)
    heights = all_heights[lo:hi]
    mixed = secondary_ratio(heights) > params.mixed_peak_ratio
    flank = np.ones(hi - lo, dtype=bool)
    flank[start - lo : end - lo] = False
    editable = np.zeros(hi - lo, dtype=bool)
    editable[start - lo : end - lo] = [b == EDITABLE_BASE[editor] for b in guide]

    # Read-level signal and quality are judged on the region where background is measured,
    # using the median, so that one dead flank does not condemn a read whose guide is clean.
    background_region = background_positions(chrom, start, len(guide), params)

    qc = QCResult(True)
    qc.background_median_peak_height = float(np.median(all_heights[background_region].max(axis=1)))
    qc.flank_mixed_fraction = float(mixed[flank].mean())
    qc.window_mixed_fraction = float(mixed[~editable].mean())
    if chrom.phred is not None:
        qc.background_median_phred = float(np.median(chrom.phred[background_region]))
    qc.reliable, qc.local_mixed = position_reliability(chrom, start, guide, editor, params)

    if background_region.size < params.min_background_positions:
        qc.reasons.append(f"only {background_region.size} positions available to measure background")
    if qc.background_median_phred is not None and qc.background_median_phred < params.min_background_phred:
        qc.reasons.append(f"low base quality around the guide (median Phred {qc.background_median_phred:.0f})")
    if qc.background_median_peak_height < params.min_peak_height:
        qc.reasons.append(f"weak signal around the guide (median peak {qc.background_median_peak_height:.0f} RFU)")
    if qc.window_mixed_fraction > params.max_window_mixed_fraction:
        qc.reasons.append(
            f"overlapping peaks across the window ({qc.window_mixed_fraction:.0%} of non-editable positions)"
        )
    targets = target_positions(guide, editor)
    noisy = qc.unreliable_labels(guide, editor)
    if noisy and not params.partial_reads:
        qc.reasons.append(f"overlapping peaks near {', '.join(noisy)}")
    elif targets and not qc.reliable[targets].any():
        qc.reasons.append("no target position has a clean enough neighbourhood to be quantified")
    qc.passed = not qc.reasons
    return qc


# ---------------------------------------------------------------- quantification


@dataclass
class Quantification:
    percentages: np.ndarray  # (guide length, 4) percent of A, C, G, T at each guide position
    background: dict[str, float]
    background_positions: int


def mad_filter(values: np.ndarray, cutoff: float) -> np.ndarray:
    """Drop values further than ``cutoff`` scaled MADs from the median."""
    if values.size < 3:
        return values
    median = np.median(values)
    mad = 1.4826 * np.median(np.abs(values - median))
    if mad == 0:
        return values
    return values[np.abs(values - median) <= cutoff * mad]


def estimate_background(
    heights: np.ndarray, phred: np.ndarray | None, positions: np.ndarray, params: Parameters
) -> tuple[dict[str, float], int]:
    """Mean non-primary signal per channel over ``positions`` after MAD outlier removal."""
    if phred is not None:
        good = positions[phred[positions] >= params.min_background_phred]
        if good.size >= 10:
            positions = good
    sub = heights[positions]
    primary = sub.argmax(axis=1)
    background = {}
    for column, base in enumerate(BASES):
        values = mad_filter(sub[primary != column, column], params.mad_cutoff)
        background[base] = float(values.mean()) if values.size else 0.0
    return background, int(positions.size)


def background_positions(chrom: Chromatogram, start: int, guide_length: int, params: Parameters) -> np.ndarray:
    """Positions flanking the guide that background is measured on. The guide is never included."""
    end = start + guide_length
    return np.r_[
        max(0, start - params.background_flank) : start,
        end : min(len(chrom), end + params.background_flank),
    ]


def quantify(chrom: Chromatogram, start: int, guide_length: int, params: Parameters) -> Quantification:
    """Background-corrected base percentages along the guide of an oriented chromatogram."""
    heights = chrom.peak_heights(params.peak_half_width)
    end = start + guide_length
    flank_positions = background_positions(chrom, start, guide_length, params)
    background, n_positions = estimate_background(heights, chrom.phred, flank_positions, params)

    signal = heights[start:end] - np.array([background[b] for b in BASES])
    signal = np.clip(signal, 0, None)
    totals = signal.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        percentages = np.where(totals > 0, 100 * signal / totals, np.nan)
    return Quantification(percentages, background, n_positions)


def target_positions(guide: str, editor: str) -> list[int]:
    """0-based guide positions carrying the base the editor converts."""
    return [i for i, b in enumerate(guide) if b == EDITABLE_BASE[editor]]


def position_label(guide: str, index: int) -> str:
    return f"{guide[index]}{index + 1}"


def editing_values(
    quant: Quantification, guide: str, editor: str, reliable: np.ndarray | None = None
) -> dict[str, float]:
    """Percent conversion (A->G for ABE, C->T for CBE) at each target position.

    Positions whose neighbourhood is too noisy to quantify are returned as NaN.
    """
    product = BASES.index(PRODUCT_BASE[editor])
    values = {}
    for i in target_positions(guide, editor):
        usable = reliable is None or bool(reliable[i])
        values[position_label(guide, i)] = float(quant.percentages[i, product]) if usable else float("nan")
    return values
