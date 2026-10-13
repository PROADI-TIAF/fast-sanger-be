"""Entry point for the browser version, which runs this package with Pyodide.

The page writes the chromatograms (already under their convention names) and guides.txt
into a folder of the in-memory file system, calls :func:`run`, shows the returned summary
and offers the zipped results folder for download.
"""

from __future__ import annotations

import csv
import io
import json
import logging
import math
import shutil
import zipfile
from dataclasses import fields
from pathlib import Path

from fast_sanger_be import __version__
from fast_sanger_be.analysis import position_label, target_positions
from fast_sanger_be.batch import run_analysis
from fast_sanger_be.config import PRODUCT_BASE, Parameters
from fast_sanger_be.naming import find_conflicts, list_ab1_files, parse_sample_files, read_guide_file

log = logging.getLogger("fast_sanger_be")


class WebInputError(Exception):
    """The inputs sent by the page cannot be analysed."""


def parameters_from_dict(values: dict) -> Parameters:
    """Build Parameters from the page's settings, ignoring unknown keys and casting to each field's type."""
    defaults = Parameters()
    kwargs = {}
    for f in fields(Parameters):
        if f.name in values and values[f.name] not in (None, ""):
            kind = type(getattr(defaults, f.name))
            kwargs[f.name] = values[f.name] if kind is bool else kind(values[f.name])
    return Parameters(**kwargs)


def _number(value):
    return None if value is None or (isinstance(value, float) and math.isnan(value)) else round(float(value), 4)


def _relative(path: Path | None, root: Path) -> str | None:
    return None if path is None else path.relative_to(root).as_posix()


class _CallbackHandler(logging.Handler):
    def __init__(self, callback):
        super().__init__()
        self.callback = callback

    def emit(self, record: logging.LogRecord) -> None:
        self.callback(record.levelname, self.format(record))


def run(input_dir: str, output_dir: str, settings_json: str = "{}", names_json: str = "{}", on_log=None) -> str:
    """Analyse every .ab1 in ``input_dir`` and write the results folder to ``output_dir``.

    ``names_json`` maps each analysed file name to the name the user's file had, so the
    results record which original file each read came from. ``on_log(level, message)`` is
    called for every progress message. Returns a JSON summary.
    """
    input_path, out_root = Path(input_dir), Path(output_dir)
    params = parameters_from_dict(json.loads(settings_json))
    original_names: dict[str, str] = json.loads(names_json)

    files, invalid = parse_sample_files(list_ab1_files(input_path))
    if invalid:
        raise WebInputError("; ".join(f"{p.name}: {reason}" for p, reason in invalid.items()))
    if not files:
        raise WebInputError("no .ab1 files to analyse")
    conflicts = find_conflicts(files)
    if conflicts:
        raise WebInputError("several files provide the same read: " + "; ".join(conflicts))
    guides, errors = read_guide_file(input_path / "guides.txt")
    missing = sorted({g for f in files for g in f.guides} - set(guides))
    if errors or missing:
        raise WebInputError("; ".join(errors + [f"no sequence for guide '{g}'" for g in missing]))
    used = {g for f in files for g in f.guides}
    guides = {name: g for name, g in guides.items() if name in used}

    if out_root.exists():
        shutil.rmtree(out_root)
    out_root.mkdir(parents=True)
    log_buffer = io.StringIO()
    handler = logging.StreamHandler(log_buffer)
    handler.setFormatter(logging.Formatter("%(message)s"))
    log.setLevel(logging.INFO)
    log.handlers.clear()
    log.addHandler(handler)
    if on_log is not None:
        live = _CallbackHandler(on_log)
        live.setFormatter(logging.Formatter("%(message)s"))
        log.addHandler(live)
    try:
        results, summaries = run_analysis(files, guides, params, out_root, "browser upload")
    finally:
        log.handlers.clear()
    (out_root / "run_log.txt").write_text(log_buffer.getvalue(), encoding="utf-8")
    with (out_root / "input_files.tsv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["analysed_as", "original_file"])
        for f in sorted(files, key=lambda f: f.path.name):
            writer.writerow([f.path.name, original_names.get(f.path.name, f.path.name)])

    editors = {}
    for editor, by_guide in summaries.items():
        guide_blocks = []
        for name, rows in by_guide.items():
            guide = guides[name]
            labels = [position_label(guide.sequence, i) for i in target_positions(guide.sequence, editor)]
            reads = [r for r in results if r.file.editor == editor and r.guide.name == name]
            heatmap = out_root / editor / "Heatmaps" / f"{name}_heatmap.{params.figure_format}"
            guide_blocks.append({
                "name": name,
                "sequence": guide.sequence,
                "product": PRODUCT_BASE[editor],
                "labels": labels,
                "samples": [
                    {
                        "sample": s.sample,
                        "reads_used": s.reads_used,
                        "values": {label: _number(s.values.get(label)) for label in labels},
                        "max_difference": _number(s.max_difference),
                        "notes": s.notes,
                    }
                    for s in rows
                ],
                "reads": [
                    {
                        "file": r.file.path.name,
                        "original_file": original_names.get(r.file.path.name, r.file.path.name),
                        "sample": r.file.sample,
                        "read": r.file.read,
                        "status": r.status,
                        "used": r.used,
                        "reasons": r.reasons,
                        "strand": r.hit.strand if r.hit else None,
                        "figure": _relative(r.figure, out_root),
                    }
                    for r in reads
                ],
                "heatmap": _relative(heatmap, out_root) if heatmap.exists() else None,
            })
        editors[editor] = guide_blocks

    counts: dict[str, int] = {}
    for r in results:
        counts[r.status] = counts.get(r.status, 0) + 1
    return json.dumps({
        "version": __version__,
        "counts": counts,
        "editors": editors,
        "log": log_buffer.getvalue(),
        "figure_format": params.figure_format,
    })


def zip_folder(folder: str, zip_path: str) -> str:
    """Zip ``folder`` (keeping its name as the top-level directory) and return the zip path."""
    root = Path(folder)
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(root.rglob("*")):
            if path.is_file():
                zf.write(path, Path(root.name) / path.relative_to(root))
    return zip_path
