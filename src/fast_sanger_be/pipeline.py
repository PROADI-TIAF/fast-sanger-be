"""Run guide location, QC and quantification over all files, then combine forward and reverse reads."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from fast_sanger_be.analysis import (
    GuideHit,
    QCResult,
    Quantification,
    assess_quality,
    editing_values,
    locate_guide,
    orient,
    position_label,
    quantify,
    target_positions,
)
from fast_sanger_be.chromatogram import Ab1Error, read_ab1
from fast_sanger_be.config import READS, Parameters
from fast_sanger_be.naming import Guide, SampleFile
from fast_sanger_be.plotting import plot_read

log = logging.getLogger("fast_sanger_be")

OK = "OK"
LOW_QUALITY = "Low quality"
NOT_FOUND = "Guide not found"
UNREADABLE = "Unreadable file"
FAILED = "Analysis error"
MISSING = "Missing"


@dataclass
class ReadResult:
    file: SampleFile
    guide: Guide
    status: str
    reasons: list[str] = field(default_factory=list)
    hit: GuideHit | None = None
    qc: QCResult | None = None
    quant: Quantification | None = None
    values: dict[str, float] = field(default_factory=dict)
    figure: Path | None = None

    @property
    def used(self) -> bool:
        return self.status == OK

    @property
    def status_text(self) -> str:
        return f"{self.status}: {'; '.join(self.reasons)}" if self.reasons else self.status

    @property
    def partial(self) -> bool:
        return self.used and bool(self.values) and any(np.isnan(v) for v in self.values.values())


@dataclass
class SampleSummary:
    sample: str
    reads: dict[str, ReadResult | None]
    values: dict[str, float]
    reads_used: str
    max_difference: float | None
    notes: list[str]


def natural_key(text: str) -> list:
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", text)]


def raw_dir(out_root: Path, editor: str) -> Path:
    return out_root / editor / "Raw_Data"


def _quantified_mask(guide: str, editor: str, qc: QCResult) -> np.ndarray | None:
    """Per guide position: False only for target positions that QC could not quantify.

    Non-target positions stay visible in the figure: their mixed peaks are the evidence
    for why a neighbouring target was not reported.
    """
    if qc.reliable is None:
        return None
    mask = np.ones(len(guide), dtype=bool)
    for i in target_positions(guide, editor):
        mask[i] = bool(qc.reliable[i])
    return mask


def analyse_read(
    sample_file: SampleFile, chrom, guide: Guide, params: Parameters, out_root: Path
) -> ReadResult:
    result = ReadResult(sample_file, guide, OK)
    hit = locate_guide(chrom.bases, guide.sequence, sample_file.editor, params.max_guide_mismatches)
    if hit is None:
        result.status = NOT_FOUND
        result.reasons = [f"no match with <= {params.max_guide_mismatches} mismatches on either strand"]
        return result
    result.hit = hit
    oriented = orient(chrom, hit)

    qc = assess_quality(oriented, hit.start, guide.sequence, sample_file.editor, params)
    if hit.ambiguous:
        qc.passed = False
        qc.reasons.insert(0, "guide matches more than one location in the read")
    result.qc = qc
    if not qc.passed:
        result.status = LOW_QUALITY
        result.reasons = list(qc.reasons)

    covered = hit.start + len(guide.sequence) <= len(oriented)
    if covered:
        result.quant = quantify(oriented, hit.start, len(guide.sequence), params)
        if result.used:
            result.values = editing_values(result.quant, guide.sequence, sample_file.editor, qc.reliable)
            skipped = qc.unreliable_labels(guide.sequence, sample_file.editor)
            if skipped:
                result.reasons.append(f"not quantified at {', '.join(skipped)} (overlapping peaks nearby)")

    folder = "Figures" if result.used else "Low_Quality"
    path = raw_dir(out_root, sample_file.editor) / folder / guide.name / f"{sample_file.path.stem}.{params.figure_format}"
    strand = "guide on forward strand of read" if hit.strand == "+" else "guide on reverse complement of read"
    subtitle = f"{guide.name} 5'-{guide.sequence}-3' | {strand} | {result.status_text}"
    if covered:
        plot_read(
            oriented, hit.start, guide.sequence,
            result.quant.percentages if result.quant else None,
            sample_file.path.name, subtitle, path, params,
            reliable=_quantified_mask(guide.sequence, sample_file.editor, qc),
        )
        result.figure = path
    return result


def analyse_files(
    files: list[SampleFile], guides: dict[str, Guide], params: Parameters, out_root: Path
) -> list[ReadResult]:
    results: list[ReadResult] = []
    for index, sample_file in enumerate(files, start=1):
        log.info("[%d/%d] %s", index, len(files), sample_file.path.name)
        try:
            chrom = read_ab1(sample_file.path)
        except Ab1Error as exc:
            for name in sample_file.guides:
                results.append(ReadResult(sample_file, guides[name], UNREADABLE, [str(exc)]))
            log.warning("    unreadable: %s", exc)
            continue
        for name in sample_file.guides:
            try:
                result = analyse_read(sample_file, chrom, guides[name], params, out_root)
            except Exception as exc:  # one bad read must not stop the batch
                log.exception("    %s: unexpected error", name)
                result = ReadResult(sample_file, guides[name], FAILED, [f"{type(exc).__name__}: {exc}"])
            results.append(result)
            if result.used:
                note = "; ".join(result.reasons)
                log.info("    %s: OK (strand %s)%s", name, result.hit.strand, f" - {note}" if note else "")
            else:
                log.warning("    %s: %s", name, result.status_text)
    return results


def summarise_guide(results: list[ReadResult], guide: Guide, editor: str, params: Parameters) -> list[SampleSummary]:
    """Combine forward and reverse reads of each sample for one guide (mean of the reads that passed QC)."""
    labels = [position_label(guide.sequence, i) for i in target_positions(guide.sequence, editor)]
    by_sample: dict[str, dict[str, ReadResult]] = {}
    for r in results:
        if r.guide.name == guide.name and r.file.editor == editor:
            by_sample.setdefault(r.file.sample, {})[r.file.read] = r

    summaries = []
    for sample in sorted(by_sample, key=natural_key):
        reads = {read: by_sample[sample].get(read) for read in READS}
        used = [r for r in reads.values() if r is not None and r.used]
        notes = []
        for read, r in reads.items():
            if r is None:
                notes.append(f"{read} read missing")
            elif not r.used:
                notes.append(f"{read} read excluded ({r.status_text})")
            elif r.partial:
                skipped = [label for label in labels if np.isnan(r.values.get(label, np.nan))]
                notes.append(f"{read} read not quantified at {', '.join(skipped)} (overlapping peaks nearby)")

        values: dict[str, float] = {}
        max_difference = None
        if used:
            for label in labels:
                read_values = [r.values[label] for r in used if not np.isnan(r.values.get(label, np.nan))]
                values[label] = float(np.mean(read_values)) if read_values else float("nan")
        if len(used) == 2:
            differences = {
                label: abs(used[0].values[label] - used[1].values[label])
                for label in labels
                if not np.isnan(used[0].values.get(label, np.nan))
                and not np.isnan(used[1].values.get(label, np.nan))
            }
            if differences:
                max_difference = max(differences.values())
                discordant = [label for label, d in differences.items() if d > params.discordance_threshold]
                if discordant:
                    notes.append(
                        f"F/R differ by more than {params.discordance_threshold:g} pp at {', '.join(discordant)}"
                    )
        reads_used = "+".join(r.file.read for r in used) if used else "None"
        summaries.append(SampleSummary(sample, reads, values, reads_used, max_difference, notes))
    return summaries
