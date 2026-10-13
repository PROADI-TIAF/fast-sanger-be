"""Synthetic Sanger chromatograms with known editing, for testing and benchmarking.

Usage::

    python -m fast_sanger_be.simulate OUTPUT_FOLDER

writes an example dataset (.ab1 files, guides.txt and truth.tsv with the simulated
editing at every target position).
"""

from __future__ import annotations

import argparse
import csv
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from fast_sanger_be.chromatogram import COMPLEMENT

DYE_ORDER = "GATC"
PEAK_SPACING = 12
PEAK_SIGMA = 2.4


# ---------------------------------------------------------------- ABIF writer


def _entry(name: str, number: int, elem_type: int, elem_size: int, data: bytes, offset: int) -> tuple[bytes, bytes]:
    count = len(data) // elem_size
    if len(data) <= 4:
        inline = data.ljust(4, b"\0")
        return struct.pack(">4sIHHII", name.encode(), number, elem_type, elem_size, count, len(data)) + inline + b"\0\0\0\0", b""
    return struct.pack(">4sIHHIIII", name.encode(), number, elem_type, elem_size, count, len(data), offset, 0), data


def write_abif(path: Path, traces: dict[str, np.ndarray], bases: str, locations: list[int],
               quality: list[int], sample_name: str) -> None:
    """Write a minimal ABIF file readable by Biopython and most chromatogram viewers."""
    tags: list[tuple[str, int, int, int, bytes]] = []
    for number, base in zip((9, 10, 11, 12), DYE_ORDER):
        values = np.clip(np.round(traces[base]), 0, 32767).astype(">i2")
        tags.append(("DATA", number, 4, 2, values.tobytes()))
    tags.append(("FWO_", 1, 2, 1, DYE_ORDER.encode()))
    for number in (1, 2):
        tags.append(("PBAS", number, 2, 1, bases.encode()))
        tags.append(("PLOC", number, 4, 2, np.asarray(locations, dtype=">i2").tobytes()))
    tags.append(("PCON", 2, 2, 1, bytes(quality)))
    name = sample_name.encode()[:255]
    tags.append(("SMPL", 1, 18, 1, bytes([len(name)]) + name))

    header_size = 128
    blob, entries = b"", b""
    for name_, number, elem_type, elem_size, data in tags:
        entry, payload = _entry(name_, number, elem_type, elem_size, data, header_size + len(blob))
        entries += entry
        blob += payload
    dir_offset = header_size + len(blob)
    header = b"ABIF" + struct.pack(">H4sIHHIIII", 101, b"tdir", 1, 1023, 28, len(tags), len(entries), dir_offset, 0)
    path.write_bytes(header.ljust(header_size, b"\0") + blob + entries)


# ---------------------------------------------------------------- chromatogram simulation


@dataclass
class SimulatedRead:
    traces: dict[str, np.ndarray]
    bases: str
    locations: list[int]
    quality: list[int]


def simulate_read(composition: list[dict[str, float]], rng: np.random.Generator, noise: float = 12.0,
                  crosstalk: float = 0.03, contaminant: list[dict[str, float]] | None = None,
                  contaminant_fraction: float = 0.0) -> SimulatedRead:
    """Simulate a read whose position i carries base fractions ``composition[i]``."""
    n = len(composition)
    length = (n + 2) * PEAK_SPACING
    x = np.arange(length)
    traces = {b: np.abs(rng.normal(0, noise, length)) for b in "ACGT"}
    locations = []
    heights = rng.lognormal(np.log(1000), 0.25, n)
    for i in range(n):
        centre = (i + 1) * PEAK_SPACING + int(rng.integers(-1, 2))
        locations.append(centre)
        shape = np.exp(-0.5 * ((x - centre) / PEAK_SIGMA) ** 2)
        mixture = dict.fromkeys("ACGT", 0.0)
        for base, frac in composition[i].items():
            mixture[base] += frac * (1 - contaminant_fraction)
        if contaminant:
            for base, frac in contaminant[i].items():
                mixture[base] += frac * contaminant_fraction
        for base in "ACGT":
            amplitude = heights[i] * (mixture[base] + crosstalk * rng.uniform(0.5, 1.5))
            traces[base] += amplitude * shape
    bases, quality = [], []
    for i, loc in enumerate(locations):
        values = {b: traces[b][loc] for b in "ACGT"}
        ordered = sorted(values.values())
        bases.append(max(values, key=values.get))
        ratio = ordered[-2] / ordered[-1]
        quality.append(int(np.clip(50 - 60 * ratio, 5, 50)))
    return SimulatedRead(traces, "".join(bases), locations, quality)


def reverse_composition(composition: list[dict[str, float]]) -> list[dict[str, float]]:
    return [{COMPLEMENT[b]: f for b, f in pos.items()} for pos in reversed(composition)]


def edited_composition(reference: str, guide_start: int, guide: str, editor: str,
                       efficiencies: dict[int, float], guide_strand: str = "+") -> list[dict[str, float]]:
    """Composition of the amplicon (reference orientation) after editing.

    ``efficiencies`` maps 0-based guide positions to edited fractions.
    """
    source, product = ("A", "G") if editor == "ABE" else ("C", "T")
    composition = [{b: 1.0} for b in reference]
    for pos, eff in efficiencies.items():
        if guide_strand == "+":
            idx, ref_base, src, prod = guide_start + pos, reference[guide_start + pos], source, product
        else:
            idx = guide_start + len(guide) - 1 - pos
            ref_base, src, prod = reference[idx], COMPLEMENT[source], COMPLEMENT[product]
        assert ref_base == src, "edited position must carry the editable base"
        composition[idx] = {src: 1 - eff, prod: eff}
    return composition


def make_example(out_dir: Path, seed: int = 7) -> None:
    """Example dataset: two guides, ABE and CBE, multi-guide sample, a missing read and a noisy read."""
    rng = np.random.default_rng(seed)
    out_dir.mkdir(parents=True, exist_ok=True)

    def random_seq(n):
        return "".join(rng.choice(list("ACGT"), n))

    guides = {"sgHBB": "CTTGACACAACTGTGTTCAC", "sgEMX1": "GAGTCCGAGCAGAAGAAGAA"}
    amplicons = {}
    for name, guide in guides.items():
        amplicons[name] = random_seq(150) + guide + "TGG" + random_seq(180)
    # A third amplicon carrying both guides, sgEMX1 on the minus strand.
    rc = "".join(COMPLEMENT[b] for b in reversed(guides["sgEMX1"]))
    both = random_seq(140) + guides["sgHBB"] + "AGG" + random_seq(25) + "CCA" + rc + random_seq(160)

    truth_rows = []

    def emit(sample, guide_names, editor, reference, edits, reads=("F", "R"), noisy_read=None):
        composition = [{b: 1.0} for b in reference]
        for guide_name, (start, strand, eff) in edits.items():
            comp = edited_composition(reference, start, guides[guide_name], editor, eff, strand)
            for i, c in enumerate(comp):
                if c != {reference[i]: 1.0}:
                    composition[i] = c
            source = "A" if editor == "ABE" else "C"
            for pos, base in enumerate(guides[guide_name]):
                if base == source:
                    truth_rows.append([sample, guide_name, editor, f"{base}{pos + 1}", 100 * eff.get(pos, 0.0)])
        for read in reads:
            comp = composition if read == "F" else reverse_composition(composition)
            contaminant = None
            fraction = 0.0
            if read == noisy_read:
                shifted = comp[7:] + comp[:7]
                contaminant, fraction = shifted, 0.45
            sim = simulate_read(comp, rng, contaminant=contaminant, contaminant_fraction=fraction)
            name = f"{sample}_{'_'.join(guide_names)}_{editor}_{read}.ab1"
            write_abif(out_dir / name, sim.traces, sim.bases, sim.locations, sim.quality, sample)

    for i, level in enumerate((0.0, 0.15, 0.45, 0.8), start=1):
        hbb = guides["sgHBB"]
        eff = {p: level * w for p, w in zip([k for k, b in enumerate(hbb) if b == "A"], (0.3, 1.0, 0.8, 0.5, 0.2, 0.1))}
        emit(f"ABE-{i}", ["sgHBB"], "ABE", amplicons["sgHBB"], {"sgHBB": (150, "+", eff)},
             reads=("F", "R") if i != 3 else ("F",))
    for i, level in enumerate((0.1, 0.6), start=1):
        emx = guides["sgEMX1"]
        eff = {p: level * w for p, w in zip([k for k, b in enumerate(emx) if b == "C"], (1.0, 0.6))}
        emit(f"CBE-{i}", ["sgEMX1"], "CBE", amplicons["sgEMX1"], {"sgEMX1": (150, "+", eff)},
             noisy_read="R" if i == 2 else None)
    hbb_a = [k for k, b in enumerate(guides["sgHBB"]) if b == "A"]
    emx_a = [k for k, b in enumerate(guides["sgEMX1"]) if b == "A"]
    emit("Dual-1", ["sgHBB", "sgEMX1"], "ABE", both, {
        "sgHBB": (140, "+", {p: 0.4 for p in hbb_a[:3]}),
        "sgEMX1": (140 + 20 + 3 + 25 + 3, "-", {p: 0.25 for p in emx_a[:4]}),
    })

    (out_dir / "guides.txt").write_text("".join(f">{n}\n{s}\n" for n, s in guides.items()), encoding="utf-8")
    with (out_dir / "truth.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sample", "guide", "editor", "position", "simulated_editing_percent"])
        writer.writerows(truth_rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Write a synthetic example dataset for Fast-Sanger-BE.")
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args(argv)
    make_example(args.output_dir, args.seed)
    print(f"Example dataset written to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
