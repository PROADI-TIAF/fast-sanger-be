"""Prepare the web version for serving: package wheels and example data.

    python web/build.py

Writes web/wheels/ (Fast-Sanger-BE plus openpyxl and its dependency, which Pyodide does
not ship) with a manifest, and web/examples/ with a small synthetic dataset.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

WEB = Path(__file__).resolve().parent
ROOT = WEB.parent
EXTRA_WHEELS = ["openpyxl==3.1.5", "et-xmlfile==2.0.0"]


def build_wheels() -> list[str]:
    target = WEB / "wheels"
    shutil.rmtree(target, ignore_errors=True)
    target.mkdir()
    pip = [sys.executable, "-m", "pip", "--disable-pip-version-check", "-q"]
    subprocess.run(pip + ["wheel", "--no-deps", "-w", str(target), str(ROOT)], check=True)
    shutil.rmtree(ROOT / "build", ignore_errors=True)  # setuptools leftovers
    for leftover in (ROOT / "src").glob("*.egg-info"):
        shutil.rmtree(leftover, ignore_errors=True)
    subprocess.run(pip + ["download", "--no-deps", "--only-binary", ":all:", "-d", str(target), *EXTRA_WHEELS],
                   check=True)
    wheels = sorted(p.name for p in target.glob("*.whl"))
    version = next(w.split("-")[1] for w in wheels if w.startswith("fast_sanger_be-"))
    (target / "manifest.json").write_text(json.dumps({"version": version, "wheels": wheels}, indent=1) + "\n")
    return wheels


def build_examples() -> list[str]:
    sys.path.insert(0, str(ROOT / "src"))
    from fast_sanger_be.simulate import make_example

    target = WEB / "examples"
    shutil.rmtree(target, ignore_errors=True)
    target.mkdir()
    with tempfile.TemporaryDirectory() as tmp:
        make_example(Path(tmp))
        for path in Path(tmp).iterdir():
            if path.suffix == ".ab1" or path.name == "guides.txt":
                shutil.copy(path, target / path.name)
    names = sorted(p.name for p in target.iterdir())
    (target / "manifest.json").write_text(json.dumps({"files": names}, indent=1) + "\n")
    return names


if __name__ == "__main__":
    print("wheels:", ", ".join(build_wheels()))
    print("examples:", len(build_examples()), "files")
