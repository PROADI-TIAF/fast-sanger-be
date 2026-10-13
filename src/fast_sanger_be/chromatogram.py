"""Reading .ab1 files into a strand-aware chromatogram."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Seq import reverse_complement

from fast_sanger_be.config import BASES

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


class Ab1Error(Exception):
    """The file cannot be used as a Sanger chromatogram."""


@dataclass(frozen=True)
class Chromatogram:
    """Base calls, peak locations, per-base quality and the four dye traces of one read."""

    name: str
    bases: str
    peak_locations: np.ndarray  # trace index of each base call
    traces: dict[str, np.ndarray]  # A, C, G, T -> processed trace
    phred: np.ndarray | None  # per-base quality, None if absent

    def __len__(self) -> int:
        return len(self.bases)

    @property
    def trace_length(self) -> int:
        return len(self.traces["A"])

    def peak_heights(self, half_width: int = 1) -> np.ndarray:
        """Maximum of each channel around every peak location. Shape (n_bases, 4), columns ACGT."""
        stacked = np.vstack([self.traces[b] for b in BASES])
        heights = np.empty((len(self.bases), 4))
        last = self.trace_length - 1
        for i, loc in enumerate(self.peak_locations):
            lo, hi = max(0, loc - half_width), min(last, loc + half_width)
            heights[i] = stacked[:, lo : hi + 1].max(axis=1)
        return heights

    def reverse_complement(self) -> Chromatogram:
        last = self.trace_length - 1
        return Chromatogram(
            name=self.name,
            bases=reverse_complement(self.bases),
            peak_locations=(last - self.peak_locations)[::-1].copy(),
            traces={b: self.traces[COMPLEMENT[b]][::-1].copy() for b in BASES},
            phred=None if self.phred is None else self.phred[::-1].copy(),
        )


def _as_text(value) -> str:
    return value.decode("ascii", errors="replace") if isinstance(value, bytes) else str(value)


def read_ab1(path: Path) -> Chromatogram:
    """Read an ABIF file. Dye order is taken from the FWO_ tag instead of being assumed."""
    try:
        record = SeqIO.read(str(path), "abi")
    except Exception as exc:  # Biopython raises several exception types for corrupt files
        raise Ab1Error(f"cannot be read as an ABIF file ({exc})") from exc
    raw = record.annotations.get("abif_raw", {})

    dye_order = _as_text(raw.get("FWO_1", b"GATC")).upper()
    if sorted(dye_order) != sorted(BASES):
        raise Ab1Error(f"unexpected dye order '{dye_order}'")
    traces = {}
    for channel, base in zip(("DATA9", "DATA10", "DATA11", "DATA12"), dye_order):
        if channel not in raw:
            raise Ab1Error(f"trace {channel} is missing")
        traces[base] = np.asarray(raw[channel], dtype=float)

    bases = raw.get("PBAS2") or raw.get("PBAS1")
    locations = raw.get("PLOC2") or raw.get("PLOC1")
    if not bases or not locations:
        raise Ab1Error("base calls or peak locations are missing")
    bases = _as_text(bases).upper()
    locations = np.atleast_1d(np.asarray(locations, dtype=int))

    quality = raw.get("PCON2") or raw.get("PCON1")
    phred = None
    if quality:
        phred = np.frombuffer(quality, dtype=np.uint8).astype(float) if isinstance(quality, bytes) \
            else np.asarray(quality, dtype=float)

    n = min(len(bases), len(locations), len(phred) if phred is not None else len(bases))
    if n < 50:
        raise Ab1Error(f"only {n} base calls")
    locations = np.clip(locations[:n], 0, len(traces["A"]) - 1)
    return Chromatogram(
        name=path.name,
        bases=bases[:n],
        peak_locations=locations,
        traces=traces,
        phred=None if phred is None else phred[:n],
    )
