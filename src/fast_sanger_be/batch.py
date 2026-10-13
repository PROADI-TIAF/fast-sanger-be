"""Quality control, quantification and output writing for a set of named files (shared by the CLI and the web version)."""

from __future__ import annotations

import logging
from collections import Counter
from pathlib import Path

from fast_sanger_be import __version__
from fast_sanger_be.config import EDITORS, Parameters
from fast_sanger_be.naming import Guide, SampleFile
from fast_sanger_be.pipeline import OK, ReadResult, SampleSummary, analyse_files, natural_key
from fast_sanger_be.reporting import run_info, write_editor_outputs

log = logging.getLogger("fast_sanger_be")


def run_analysis(
    files: list[SampleFile], guides: dict[str, Guide], params: Parameters, out_root: Path, input_label: Path | str
) -> tuple[list[ReadResult], dict[str, dict[str, list[SampleSummary]]]]:
    """Analyse every file against its guides and write all outputs under ``out_root``.

    Returns the per-read results and, for each editor, the sample summaries of each guide.
    """
    log.info("\n== Step 3/4  Quality control and quantification")
    log.info("Fast-Sanger-BE %s, %d file(s), %d guide(s)", __version__, len(files), len(guides))
    files = sorted(files, key=lambda f: (f.editor, natural_key(f.sample), f.read))
    results = analyse_files(files, guides, params, out_root)

    log.info("\n== Step 4/4  Writing results")
    info = run_info(input_label, len(files))
    summaries: dict[str, dict[str, list[SampleSummary]]] = {}
    for editor in EDITORS:
        editor_results = [r for r in results if r.file.editor == editor]
        if not editor_results:
            continue
        used = {r.guide.name for r in editor_results}
        editor_guides = [g for name, g in guides.items() if name in used]
        summaries[editor] = write_editor_outputs(out_root, editor, editor_guides, editor_results, params, info)
        log.info("  %s: %s", editor, Path(out_root.name) / editor / f"{editor}_final_results.xlsx")

    counts = Counter(r.status for r in results)
    log.info("\nDone. %d read/guide analyses: %s", len(results),
             ", ".join(f"{n} {status}" for status, n in counts.most_common()))
    excluded = [r for r in results if r.status != OK]
    if excluded:
        log.warning("Excluded from the results (see Excluded_reads in the final workbooks):")
        for r in excluded:
            log.warning("  - %s [%s]: %s", r.file.path.name, r.guide.name, r.status_text)
    return results, summaries
