"""Prepare the web version for serving: package wheels, example data and Pyodide.

    python web/build.py [--no-pyodide]

Writes web/wheels/ (Fast-Sanger-BE plus openpyxl and its dependency, which Pyodide does
not ship) with a manifest, web/examples/ with a small synthetic dataset, and web/pyodide/
with the Pyodide runtime and the packages the analysis needs. Serving Pyodide from the
site itself means the page works on networks that block the jsDelivr CDN; when
web/pyodide/ is missing, the page falls back to the CDN.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

WEB = Path(__file__).resolve().parent
ROOT = WEB.parent
EXTRA_WHEELS = ["openpyxl==3.1.5", "et-xmlfile==2.0.0"]
PYODIDE_CORE = ["pyodide.mjs", "pyodide.asm.mjs", "pyodide.asm.wasm", "python_stdlib.zip", "pyodide-lock.json"]
PYODIDE_PACKAGES = ["numpy", "matplotlib", "biopython", "micropip"]  # must match loadPackage in worker.js


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


def pyodide_version() -> str:
    match = re.search(r'PYODIDE_VERSION = "([^"]+)"', (WEB / "worker.js").read_text())
    if not match:
        raise RuntimeError("PYODIDE_VERSION not found in web/worker.js")
    return match.group(1)


def _download(url: str, path: Path, sha256: str | None = None) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "fast-sanger-be-build"})
    with urllib.request.urlopen(request, timeout=120) as response:
        data = response.read()
    if sha256 and hashlib.sha256(data).hexdigest() != sha256:
        raise RuntimeError(f"checksum mismatch for {url}")
    path.write_bytes(data)


def _with_dependencies(packages: dict, names: list[str]) -> list[str]:
    found: list[str] = []
    pending = [n.lower() for n in names]
    while pending:
        name = pending.pop()
        if name not in found:
            found.append(name)
            pending.extend(d.lower() for d in packages[name].get("depends", []))
    return sorted(found)


def build_pyodide() -> list[str]:
    """Copy the Pyodide runtime and the needed packages from the CDN into web/pyodide/."""
    version = pyodide_version()
    base = f"https://cdn.jsdelivr.net/pyodide/v{version}/full/"
    target = WEB / "pyodide"
    shutil.rmtree(target, ignore_errors=True)
    target.mkdir()
    try:
        for name in PYODIDE_CORE:
            _download(base + name, target / name)
        packages = json.loads((target / "pyodide-lock.json").read_text())["packages"]
        names = _with_dependencies(packages, PYODIDE_PACKAGES)
        for name in names:
            entry = packages[name]
            _download(base + entry["file_name"], target / entry["file_name"], entry.get("sha256"))
    except Exception:
        shutil.rmtree(target, ignore_errors=True)
        raise
    return names


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare the web version for serving.")
    parser.add_argument("--no-pyodide", action="store_true",
                        help="do not copy Pyodide into web/pyodide/; the page then loads it from the CDN")
    args = parser.parse_args()
    print("wheels:", ", ".join(build_wheels()))
    print("examples:", len(build_examples()), "files")
    if not args.no_pyodide:
        try:
            print("pyodide packages:", ", ".join(build_pyodide()))
        except Exception as err:  # offline or CDN blocked: the page falls back to the CDN at run time
            print(f"warning: could not copy Pyodide ({err}); the page will load it from the CDN", file=sys.stderr)
