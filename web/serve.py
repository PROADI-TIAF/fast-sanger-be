"""Open the web version on this computer.

    python web/serve.py [--port 8000] [--no-browser]

Serves the web/ folder at http://localhost:8000 (building it first if needed) and opens it
in the default browser. The analysis runs in the browser; an internet connection is needed
the first time to download Python (Pyodide).
"""

from __future__ import annotations

import argparse
import functools
import subprocess
import sys
import webbrowser
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

WEB = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--no-browser", action="store_true")
    args = parser.parse_args()

    if not (WEB / "wheels" / "manifest.json").exists() or not (WEB / "examples" / "manifest.json").exists():
        subprocess.run([sys.executable, str(WEB / "build.py")], check=True)

    handler = functools.partial(SimpleHTTPRequestHandler, directory=str(WEB))
    server = ThreadingHTTPServer(("127.0.0.1", args.port), handler)
    url = f"http://localhost:{args.port}/"
    print(f"Fast-Sanger-BE web version at {url}  (Ctrl+C to stop)")
    if not args.no_browser:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
