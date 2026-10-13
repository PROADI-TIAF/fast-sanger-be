"""Start-up banner shown when Fast-Sanger-BE runs in an interactive terminal."""

from __future__ import annotations

import shutil
import sys

from fast_sanger_be import __version__

TITLE = "~~~ Fast Sanger BE ~~~"
TAGLINE = "A script for quick and easy analysis focused on base editors"

WIDE = r"""
                 _    _                                                                                   _    _
                (_\__/(,_                     _________________________________                          (_\__/(,_
                | \  _////-._                |                                 |                         | \  _////-._
 _    _         L_/__  => __/ \              |{name}|          _    _         L_/__  => __/ \
(_\__/(,_       |=====;__/___./              |   Sanger Base Editing Pipeline  |         (_\__/(,_       |=====;__/___./
| \  _////-._   '-'-'-''''''''               |_________________________________|         | \  _////-._   '-'-'-''''''''
J_/___'=> __/ \                                                                          J_/___'=> __/ \
|=====;__/___./                                                                          |=====;__/___./
'-'-'-''''''''                                                                           '-'-'-''''''''
"""

NARROW = r"""
 _    _                _________________________________
(_\__/(,_             |                                 |
| \  _////-._         |{name}|
L_/__  => __/ \       |   Sanger Base Editing Pipeline  |
|=====;__/___./       |_________________________________|
'-'-'-''''''''
"""


def render(columns: int) -> str:
    """The banner for a terminal ``columns`` wide: two pairs of sneakers if it fits, else one."""
    name = f"     FAST SANGER BE  v{__version__}".ljust(33)
    wide, narrow = (a.strip("\n").format(name=name) for a in (WIDE, NARROW))
    art = wide if columns > _width(wide) else narrow
    width = _width(art)
    return "\n".join(["", TITLE.center(width).rstrip(), TAGLINE.center(width).rstrip(), "", art, ""])


def show() -> None:
    """Print the banner, but only to a terminal so that logs and pipelines stay clean."""
    if sys.stdout.isatty():
        print(render(shutil.get_terminal_size().columns))


def _width(art: str) -> int:
    return max(len(line) for line in art.splitlines())
