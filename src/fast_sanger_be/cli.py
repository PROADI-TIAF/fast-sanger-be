"""Command-line entry point: input checks, interactive editors, analysis and reports."""

from __future__ import annotations

import argparse
import logging
import sys
from datetime import datetime
from pathlib import Path

from fast_sanger_be import __version__, banner
from fast_sanger_be.batch import run_analysis
from fast_sanger_be.config import Parameters
from fast_sanger_be.editors import UserAbort, run_editor
from fast_sanger_be.naming import (
    GUIDE_FILE_NAMES,
    RENAME_COLUMNS,
    RENAME_LOG_NAME,
    RENAME_TABLE_NAME,
    Guide,
    SampleFile,
    apply_renames,
    find_conflicts,
    find_guide_file,
    format_guides,
    guess_fields,
    list_ab1_files,
    parse_rename_table,
    parse_sample_files,
    raw_guide_entries,
    read_guide_file,
)
from fast_sanger_be.pipeline import natural_key

log = logging.getLogger("fast_sanger_be")

EXIT_OK, EXIT_ERROR, EXIT_NEEDS_INPUT = 0, 1, 2


class InputError(Exception):
    """The inputs cannot be analysed as they are."""


class NeedsInput(Exception):
    """Non-interactive run: a template was written for the user to complete."""


# ---------------------------------------------------------------- arguments


class _Parser(argparse.ArgumentParser):
    """Argument parser that shows the banner above the help text."""

    def print_help(self, file=None) -> None:
        if file is None or file is sys.stdout:
            banner.show()
        super().print_help(file)


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    defaults = Parameters()
    parser = _Parser(
        prog="fast-sanger-be",
        description="Quantify adenine (ABE) and cytosine (CBE) base editing from Sanger .ab1 files.",
        epilog="File names: <Sample>_<Guide>[_<Guide2>...]_<ABE|CBE>_<F|R>.ab1 . "
        "Guides: FASTA-like text file (>name / sequence), by default guides.txt in the input folder.",
    )
    parser.add_argument("input_dir", nargs="?", default=".", type=Path,
                        help="folder with the .ab1 files (default: current folder)")
    parser.add_argument("-g", "--guides", type=Path, help="guide file (default: guides.txt in the input folder)")
    parser.add_argument("-o", "--output-name", default="Final_Results",
                        help="name of the results folder created inside the input folder (default: Final_Results)")
    parser.add_argument("--no-browser", action="store_true", help="print the editor address instead of opening a browser")
    parser.add_argument("--non-interactive", action="store_true",
                        help="never prompt: write templates for missing information and exit with code 2")
    parser.add_argument("-V", "--version", action="version", version=f"%(prog)s {__version__}")

    qc = parser.add_argument_group("quality control")
    qc.add_argument("--qc-flank", type=int, default=defaults.qc_flank, help="nt on each side of the guide checked by QC (default: %(default)s)")
    qc.add_argument("--min-phred", type=float, default=defaults.min_background_phred, help="minimum median Phred around the guide (default: %(default)s)")
    qc.add_argument("--min-peak-height", type=float, default=defaults.min_peak_height, help="minimum median peak height around the guide, RFU (default: %(default)s)")
    qc.add_argument("--mixed-peak-ratio", type=float, default=defaults.mixed_peak_ratio, help="secondary/primary ratio that defines an overlapping peak (default: %(default)s)")
    qc.add_argument("--position-flank", type=int, default=defaults.position_flank, help="nt on each side of a target position whose neighbours decide if it is noisy (default: %(default)s)")
    qc.add_argument("--max-local-mixed", type=float, default=defaults.max_local_mixed_fraction, help="fraction of a target's neighbours with overlapping peaks above which it counts as noisy (default: %(default)s)")
    qc.add_argument("--partial-reads", action="store_true", help="keep reads with noisy target positions and leave only those positions empty (default: reject the whole read)")
    qc.add_argument("--max-window-mixed", type=float, default=defaults.max_window_mixed_fraction, help="fraction of overlapping peaks in the window above which the whole read is rejected (default: %(default)s)")
    qc.add_argument("--max-mismatches", type=int, default=defaults.max_guide_mismatches, help="mismatches allowed when locating the guide (default: %(default)s)")

    quant = parser.add_argument_group("quantification and output")
    quant.add_argument("--background-flank", type=int, default=defaults.background_flank, help="nt on each side of the guide used for background (default: %(default)s)")
    quant.add_argument("--discordance", type=float, default=defaults.discordance_threshold, help="F/R difference in percentage points that is flagged (default: %(default)s)")
    quant.add_argument("--heatmap-max", type=float, default=defaults.heatmap_max, help="upper limit of the heatmap colour scale (default: %(default)s)")
    quant.add_argument("--figure-format", choices=("png", "pdf", "svg", "tiff"), default=defaults.figure_format, help="figure format (default: %(default)s)")
    quant.add_argument("--dpi", type=int, default=defaults.figure_dpi, help="resolution of raster figures (default: %(default)s)")
    return parser.parse_args(argv)


def parameters_from(args: argparse.Namespace) -> Parameters:
    return Parameters(
        qc_flank=args.qc_flank,
        background_flank=args.background_flank,
        max_guide_mismatches=args.max_mismatches,
        min_background_phred=args.min_phred,
        min_peak_height=args.min_peak_height,
        mixed_peak_ratio=args.mixed_peak_ratio,
        position_flank=args.position_flank,
        max_local_mixed_fraction=args.max_local_mixed,
        partial_reads=args.partial_reads,
        max_window_mixed_fraction=args.max_window_mixed,
        discordance_threshold=args.discordance,
        heatmap_max=args.heatmap_max,
        figure_format=args.figure_format,
        figure_dpi=args.dpi,
    )


# ---------------------------------------------------------------- input preparation


def _step(title: str) -> None:
    print(f"\n== {title}")


def resolve_file_names(input_dir: Path, known_guides: set[str], args) -> list[SampleFile]:
    """Parse .ab1 names; let the user rename files that do not follow the convention."""
    paths = list_ab1_files(input_dir)
    if not paths:
        raise InputError(f"no .ab1 files found in {input_dir}")
    valid, invalid = parse_sample_files(paths)
    print(f"  {len(paths)} .ab1 file(s): {len(valid)} correctly named, {len(invalid)} to rename.")
    if not invalid:
        return valid

    for path, reason in invalid.items():
        print(f"    - {path.name}: {reason}")
    table_path = input_dir / RENAME_TABLE_NAME
    existing = {p.name for p in paths if p not in invalid}

    def validate(path: Path):
        return parse_rename_table(path.read_text(encoding="utf-8-sig"), input_dir, list(invalid), existing)

    if args.non_interactive:
        if table_path.exists():
            plan, errors = validate(table_path)
            if errors:
                raise InputError(f"{RENAME_TABLE_NAME} has problems:\n    - " + "\n    - ".join(errors))
        else:
            lines = ["\t".join(RENAME_COLUMNS)]
            for path in invalid:
                g = guess_fields(path.name, known_guides)
                lines.append("\t".join([path.name, "yes", g["sample"], g["guides"], g["editor"], g["read"]]))
            table_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            raise NeedsInput(f"complete {table_path} (one row per file) and run again")
    else:
        print("  Opening the rename editor.")
        rows = [{"current_file": p.name, "include": True, **guess_fields(p.name, known_guides)} for p in invalid]
        config = {"rows": rows, "known_guides": sorted(known_guides, key=natural_key), "existing_names": sorted(existing)}
        plan = run_editor("rename_editor.html", config, table_path, validate, not args.no_browser,
                          "  Rename the files and continue?")

    stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    apply_renames(plan, input_dir / RENAME_LOG_NAME, stamp)
    for old, new in plan.renames:
        print(f"    renamed {old.name} -> {new.name}")
    for path in plan.excluded:
        print(f"    skipped {path.name}")
    renamed, still_invalid = parse_sample_files([new for _, new in plan.renames])
    if still_invalid:  # cannot happen when the table validated, kept as a safeguard
        raise InputError("renamed files still do not follow the naming convention")
    return valid + renamed


def resolve_guides(input_dir: Path, files: list[SampleFile], args) -> dict[str, Guide]:
    """Load the guide file; open the guide editor when it is missing, invalid or incomplete."""
    guide_path = args.guides.resolve() if args.guides else find_guide_file(input_dir)
    needed = sorted({g for f in files for g in f.guides}, key=natural_key)
    guides, errors, entries = {}, [], []
    if guide_path and guide_path.exists():
        guides, errors = read_guide_file(guide_path)
        entries = raw_guide_entries(guide_path.read_text(encoding="utf-8-sig", errors="replace"))
        print(f"  Guide file: {guide_path.name} ({len(guides)} valid guide(s))")
    else:
        print(f"  No guide file found (looked for {', '.join(GUIDE_FILE_NAMES)}).")
        guide_path = guide_path or input_dir / "guides.txt"
    missing = [g for g in needed if g not in guides]
    if not errors and not missing:
        return guides

    for error in errors:
        print(f"    - {error}")
    if missing:
        print(f"    - no valid sequence for: {', '.join(missing)}")

    def validate(path: Path):
        parsed, problems = read_guide_file(path)
        problems += [f"no sequence for guide '{g}'" for g in needed if g not in parsed]
        return parsed, problems

    if args.non_interactive:
        if not guide_path.exists():
            guide_path.write_text("".join(f">{g}\n\n" for g in needed), encoding="utf-8")
            raise NeedsInput(f"add the sequences to {guide_path} and run again")
        raise InputError(f"fix {guide_path.name} and run again")

    print("  Opening the guide editor.")
    listed = {name for name, _ in entries}
    rows = [[name, seq] for name, seq in entries] + [[g, ""] for g in missing if g not in listed]
    return run_editor("guide_editor.html", {"rows": rows, "required": needed}, guide_path, validate,
                      not args.no_browser, "  Guide sequences saved. Continue the analysis?")


# ---------------------------------------------------------------- output


def prepare_output(input_dir: Path, name: str) -> Path:
    out_root = input_dir / name
    if out_root.exists():
        backup = input_dir / f"{name}_{datetime.now().strftime('%Y%m%d-%H%M%S')}"
        out_root.rename(backup)
        print(f"  Previous results moved to {backup.name}")
    out_root.mkdir()
    return out_root


def setup_logging(log_file: Path) -> None:
    log.setLevel(logging.INFO)
    log.handlers.clear()
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter("%(message)s"))
    file_handler = logging.FileHandler(log_file, encoding="utf-8")
    file_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)-7s %(message)s"))
    log.addHandler(console)
    log.addHandler(file_handler)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    params = parameters_from(args)
    input_dir = args.input_dir.resolve()
    banner.show()
    print(f"Fast-Sanger-BE {__version__}  |  input: {input_dir}")
    if not input_dir.is_dir():
        print(f"ERROR: {input_dir} is not a folder")
        return EXIT_ERROR

    try:
        _step("Step 1/4  Checking file names")
        guide_path = args.guides.resolve() if args.guides else find_guide_file(input_dir)
        known = set(read_guide_file(guide_path)[0]) if guide_path and guide_path.exists() else set()
        files = resolve_file_names(input_dir, known, args)
        conflicts = find_conflicts(files)
        if conflicts:
            raise InputError("several files provide the same read:\n    - " + "\n    - ".join(conflicts))

        _step("Step 2/4  Checking guide sequences")
        guides = resolve_guides(input_dir, files, args)
    except UserAbort as exc:
        print(f"\nStopped: {exc}.")
        return EXIT_NEEDS_INPUT
    except NeedsInput as exc:
        print(f"\nAction needed: {exc}.")
        return EXIT_NEEDS_INPUT
    except InputError as exc:
        print(f"\nERROR: {exc}")
        return EXIT_ERROR

    unused = [g for g in guides if g not in {name for f in files for name in f.guides}]
    if unused:
        print(f"  Guides not used by any file (ignored): {', '.join(unused)}")

    out_root = prepare_output(input_dir, args.output_name)
    setup_logging(out_root / "run_log.txt")
    run_analysis(files, guides, params, out_root, input_dir)
    log.info("Results folder: %s", out_root)
    return EXIT_OK


if __name__ == "__main__":
    raise SystemExit(main())
