"""Browser-based editors served from a local web server.

The page saves directly into the analysis folder (no trip through the Downloads
folder); the terminal waits for the user to confirm with y/n.
"""

from __future__ import annotations

import json
import secrets
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from importlib import resources
from pathlib import Path

MAX_UPLOAD_BYTES = 2_000_000


class UserAbort(Exception):
    """The user chose not to continue."""


def load_template(name: str, config: dict) -> str:
    html = resources.files("fast_sanger_be").joinpath("templates", name).read_text(encoding="utf-8")
    payload = json.dumps(config).replace("</", "<\\/")
    return html.replace("/*__CONFIG__*/null", payload)


class _Handler(BaseHTTPRequestHandler):
    server: "_EditorServer"

    def _send(self, status: int, body: bytes, content_type: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _json(self, status: int, payload: dict) -> None:
        self._send(status, json.dumps(payload).encode(), "application/json")

    def do_GET(self) -> None:  # noqa: N802
        if self.path.split("?")[0] in ("/", f"/{self.server.token}"):
            self._send(200, self.server.html.encode("utf-8"), "text/html; charset=utf-8")
        else:
            self._send(404, b"Not found", "text/plain")

    def do_POST(self) -> None:  # noqa: N802
        length = int(self.headers.get("Content-Length") or 0)
        if self.path != "/save" or length > MAX_UPLOAD_BYTES:
            self._json(400, {"error": "invalid request"})
            return
        try:
            data = json.loads(self.rfile.read(length))
        except json.JSONDecodeError:
            self._json(400, {"error": "invalid JSON"})
            return
        if data.get("token") != self.server.token or not isinstance(data.get("content"), str):
            self._json(403, {"error": "invalid token"})
            return
        self.server.target.write_text(data["content"], encoding="utf-8")
        self._json(200, {"path": str(self.server.target)})

    def log_message(self, *args) -> None:  # keep the terminal clean
        pass


class _EditorServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, html: str, target: Path, token: str):
        super().__init__(("127.0.0.1", 0), _Handler)
        self.html, self.target, self.token = html, target, token


class EditorSession:
    """Serve an editor page on localhost for as long as the ``with`` block runs."""

    def __init__(self, template: str, config: dict, target: Path, open_browser: bool = True):
        self.token = secrets.token_urlsafe(16)
        config = {**config, "token": self.token, "save_url": "/save", "target_name": target.name}
        self.server = _EditorServer(load_template(template, config), target, self.token)
        self.url = f"http://127.0.0.1:{self.server.server_port}/{self.token}"
        self.open_browser = open_browser

    def __enter__(self) -> "EditorSession":
        threading.Thread(target=self.server.serve_forever, daemon=True).start()
        opened = self.open_browser and webbrowser.open(self.url)
        if opened:
            print(f"  The editor opened in your browser. If it did not, open: {self.url}")
        else:
            print(f"  Open this address in a browser: {self.url}")
        return self

    def __exit__(self, *exc) -> None:
        self.server.shutdown()
        self.server.server_close()


def ask_yes_no(question: str) -> bool:
    while True:
        try:
            answer = input(f"{question} [y/n]: ").strip().lower()
        except EOFError:
            raise UserAbort("no answer (input closed)") from None
        if answer in ("y", "yes"):
            return True
        if answer in ("n", "no"):
            return False
        print("  Please answer y (yes) or n (no).")


def run_editor(template: str, config: dict, target: Path, validate, open_browser: bool, prompt: str):
    """Open an editor and loop until the saved file validates. Returns the validated result.

    ``validate(path)`` must return (result, errors).
    """
    with EditorSession(template, config, target, open_browser):
        while True:
            if not ask_yes_no(prompt):
                raise UserAbort("aborted by the user")
            if not target.exists():
                print(f"  {target.name} was not found next to the .ab1 files. Save it from the editor first.")
                continue
            result, errors = validate(target)
            if not errors:
                return result
            print(f"  {target.name} still has problems:")
            for error in errors:
                print(f"    - {error}")
            print("  Fix them in the editor, save again, then answer y.")
