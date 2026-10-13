"""Sample file naming convention, guide files and rename tables.

File names follow::

    <Sample>_<Guide>[_<Guide2>...]_<ABE|CBE>_<F|R>.ab1

Sample and guide names may contain letters, digits and hyphens, never underscores.
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from fast_sanger_be.config import EDITORS, READS

NAME_RE = re.compile(r"^[A-Za-z0-9-]+$")
SEQUENCE_RE = re.compile(r"^[ACGT]+$")
RESERVED_NAMES = {"ABE", "CBE", "F", "R"}
GUIDE_FILE_NAMES = ("guides.txt", "guides.fasta", "guides.fa")
RENAME_TABLE_NAME = "rename_table.tsv"
RENAME_LOG_NAME = "rename_log.tsv"
RENAME_COLUMNS = ("current_file", "include", "sample", "guides", "editor", "read")
MIN_GUIDE_LENGTH = 15
MAX_GUIDE_LENGTH = 30


@dataclass(frozen=True)
class SampleFile:
    path: Path
    sample: str
    guides: tuple[str, ...]
    editor: str
    read: str


@dataclass(frozen=True)
class Guide:
    name: str
    sequence: str


def list_ab1_files(directory: Path) -> list[Path]:
    return sorted(p for p in directory.iterdir() if p.is_file() and p.suffix.lower() == ".ab1")


def build_filename(sample: str, guides: list[str] | tuple[str, ...], editor: str, read: str) -> str:
    return f"{sample}_{'_'.join(guides)}_{editor}_{read}.ab1"


def check_name(value: str, what: str) -> str | None:
    """Return an error message if ``value`` is not a valid sample or guide name."""
    if not value:
        return f"{what} is empty"
    if not NAME_RE.match(value):
        return f"{what} '{value}' may only contain letters, digits and hyphens"
    if value.upper() in RESERVED_NAMES:
        return f"{what} '{value}' is a reserved word (ABE, CBE, F, R)"
    return None


def parse_filename(path: Path) -> tuple[SampleFile | None, str | None]:
    """Parse an .ab1 file name. Returns (SampleFile, None) or (None, reason)."""
    parts = path.stem.split("_")
    if len(parts) < 4:
        return None, "expected <Sample>_<Guide>_<ABE|CBE>_<F|R>.ab1"
    sample, guides, editor, read = parts[0], parts[1:-2], parts[-2].upper(), parts[-1].upper()
    if editor not in EDITORS:
        return None, f"'{parts[-2]}' should be ABE or CBE"
    if read not in READS:
        return None, f"'{parts[-1]}' should be F or R"
    error = check_name(sample, "sample name")
    if error:
        return None, error
    for guide in guides:
        error = check_name(guide, "guide name")
        if error:
            return None, error
    if len(set(guides)) != len(guides):
        return None, "the same guide is listed twice"
    return SampleFile(path, sample, tuple(guides), editor, read), None


def parse_sample_files(paths: list[Path]) -> tuple[list[SampleFile], dict[Path, str]]:
    valid, invalid = [], {}
    for path in paths:
        parsed, reason = parse_filename(path)
        if parsed:
            valid.append(parsed)
        else:
            invalid[path] = reason
    return valid, invalid


def find_conflicts(files: list[SampleFile]) -> list[str]:
    """Report files that would provide the same read of the same sample for the same guide."""
    seen: dict[tuple, list[str]] = defaultdict(list)
    for f in files:
        for guide in f.guides:
            seen[(f.editor, guide, f.sample, f.read)].append(f.path.name)
    return [
        f"{editor} / guide {guide} / sample {sample} / read {read}: {', '.join(names)}"
        for (editor, guide, sample, read), names in seen.items()
        if len(names) > 1
    ]


# ---------------------------------------------------------------- guide files


def find_guide_file(directory: Path) -> Path | None:
    """Find the guide file: a standard name first, otherwise the only FASTA-like text file."""
    by_lower = {p.name.lower(): p for p in directory.iterdir() if p.is_file()}
    for name in GUIDE_FILE_NAMES:
        if name in by_lower:
            return by_lower[name]
    candidates = [
        p
        for p in by_lower.values()
        if p.suffix.lower() in (".txt", ".fa", ".fasta") and _looks_like_fasta(p)
    ]
    return candidates[0] if len(candidates) == 1 else None


def _looks_like_fasta(path: Path) -> bool:
    try:
        with path.open(encoding="utf-8-sig") as fh:
            for line in fh:
                if line.strip():
                    return line.lstrip().startswith((">", "["))
    except (OSError, UnicodeDecodeError):
        return False
    return False


def raw_guide_entries(text: str) -> list[tuple[str, str]]:
    """Lenient (name, sequence) pairs, used to pre-fill the editor even when the file has errors."""
    entries: list[tuple[str, list[str], int]] = []
    _collect_entries(text, entries, [])
    return [(name, "".join("".join(lines).split()).upper()) for name, lines, _ in entries]


def _collect_entries(text: str, entries: list, errors: list[str]) -> None:
    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            entries.append((line[1:].strip(), [], line_no))
        elif line.startswith("[") and line.endswith("]"):
            entries.append((line[1:-1].strip(), [], line_no))
        elif entries:
            entries[-1][1].append(line)
        else:
            errors.append(f"line {line_no}: sequence found before any guide name")


def parse_guides(text: str) -> tuple[dict[str, Guide], list[str]]:
    """Parse a FASTA-like guide list. Headers may be '>name' or '[name]'."""
    entries: list[tuple[str, list[str], int]] = []
    errors: list[str] = []
    _collect_entries(text, entries, errors)

    guides: dict[str, Guide] = {}
    seen_lower: dict[str, str] = {}
    for name, seq_lines, line_no in entries:
        sequence = "".join("".join(seq_lines).split()).upper()
        error = check_name(name, "guide name")
        if error:
            errors.append(f"line {line_no}: {error}")
            continue
        if name.lower() in seen_lower:
            errors.append(f"line {line_no}: guide '{name}' is defined twice")
            continue
        if not sequence:
            errors.append(f"guide '{name}': sequence is missing")
            continue
        if not SEQUENCE_RE.match(sequence):
            errors.append(f"guide '{name}': sequence may only contain A, C, G and T")
            continue
        if not MIN_GUIDE_LENGTH <= len(sequence) <= MAX_GUIDE_LENGTH:
            errors.append(
                f"guide '{name}': length {len(sequence)} nt is outside "
                f"{MIN_GUIDE_LENGTH}-{MAX_GUIDE_LENGTH} nt (enter the protospacer without PAM)"
            )
            continue
        seen_lower[name.lower()] = name
        guides[name] = Guide(name, sequence)
    return guides, errors


def read_guide_file(path: Path) -> tuple[dict[str, Guide], list[str]]:
    try:
        text = path.read_text(encoding="utf-8-sig")
    except (OSError, UnicodeDecodeError) as exc:
        return {}, [f"cannot read {path.name}: {exc}"]
    return parse_guides(text)


def format_guides(guides: list[Guide]) -> str:
    return "".join(f">{g.name}\n{g.sequence}\n" for g in guides)


# ---------------------------------------------------------------- rename tables


@dataclass(frozen=True)
class RenamePlan:
    renames: list[tuple[Path, Path]]
    excluded: list[Path]


def parse_rename_table(
    text: str, directory: Path, invalid_files: list[Path], existing_names: set[str]
) -> tuple[RenamePlan | None, list[str]]:
    """Validate a rename table written by the rename editor (tab-separated, with header)."""
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    errors: list[str] = []
    if not rows or any(col not in rows[0] for col in RENAME_COLUMNS):
        return None, [f"the table must have the columns: {', '.join(RENAME_COLUMNS)}"]

    by_name = {p.name: p for p in invalid_files}
    renames: list[tuple[Path, Path]] = []
    excluded: list[Path] = []
    seen_files: set[str] = set()
    new_names: dict[str, str] = {}

    for row in rows:
        current = (row["current_file"] or "").strip()
        label = current or "(empty file name)"
        if current not in by_name:
            errors.append(f"{label}: not one of the files that need renaming")
            continue
        seen_files.add(current)
        if (row["include"] or "").strip().lower() in ("no", "n", "false", "0"):
            excluded.append(by_name[current])
            continue
        sample = (row["sample"] or "").strip()
        guides = [g.strip() for g in re.split(r"[,;+\s]+", row["guides"] or "") if g.strip()]
        editor = (row["editor"] or "").strip().upper()
        read = (row["read"] or "").strip().upper()

        row_errors = []
        if err := check_name(sample, "sample name"):
            row_errors.append(err)
        if not guides:
            row_errors.append("at least one guide name is required")
        for guide in guides:
            if err := check_name(guide, "guide name"):
                row_errors.append(err)
        if len(set(guides)) != len(guides):
            row_errors.append("the same guide is listed twice")
        if editor not in EDITORS:
            row_errors.append("editor must be ABE or CBE")
        if read not in READS:
            row_errors.append("read must be F or R")
        if row_errors:
            errors.extend(f"{current}: {e}" for e in row_errors)
            continue

        new_name = build_filename(sample, guides, editor, read)
        if new_name in new_names:
            errors.append(f"{current}: new name {new_name} is also used for {new_names[new_name]}")
            continue
        if new_name != current and new_name in existing_names:
            errors.append(f"{current}: a file named {new_name} already exists")
            continue
        new_names[new_name] = current
        renames.append((by_name[current], directory / new_name))

    for name in by_name:
        if name not in seen_files:
            errors.append(f"{name}: missing from the table")
    if errors:
        return None, errors
    return RenamePlan(renames, excluded), []


def apply_renames(plan: RenamePlan, log_path: Path, timestamp: str) -> None:
    """Rename files and append every change to a tab-separated log so it can be undone."""
    new_log = not log_path.exists()
    with log_path.open("a", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        if new_log:
            writer.writerow(["timestamp", "old_name", "new_name"])
        for old, new in plan.renames:
            if old != new:
                old.rename(new)
                writer.writerow([timestamp, old.name, new.name])


READ_WORDS = {"F": "F", "FW": "F", "FWD": "F", "FOR": "F", "FORWARD": "F",
              "R": "R", "RV": "R", "REV": "R", "REVERSE": "R"}


def guess_fields(filename: str, known_guides: set[str]) -> dict[str, str]:
    """Best-effort pre-fill of the rename editor from a non-conforming file name.

    Tokens are split on underscores, spaces and dots. Hyphenated tokens are only
    split to peel recognised words (editor, read, known guide) off their end, so
    sample names such as ``ABE-1`` stay intact.
    """
    known_lower = {g.lower(): g for g in known_guides}
    found = {"guides": [], "editor": "", "read": ""}
    rest: list[str] = []

    def classify(word: str) -> bool:
        upper = word.upper()
        if upper in EDITORS:
            found["editor"] = upper
        elif upper in READ_WORDS:
            found["read"] = READ_WORDS[upper]
        elif word.lower() in known_lower:
            found["guides"].append(known_lower[word.lower()])
        else:
            return False
        return True

    for token in (t for t in re.split(r"[_\s.]+", Path(filename).stem) if t):
        if classify(token):
            continue
        parts = token.split("-")
        while len(parts) > 1 and classify(parts[-1]):
            parts.pop()
        rest.append("-".join(parts))
    return {
        "sample": rest[0] if rest else "",
        "guides": ", ".join(dict.fromkeys(found["guides"])),
        "editor": found["editor"],
        "read": found["read"],
    }
