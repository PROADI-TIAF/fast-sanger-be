from fast_sanger_be import __version__
from fast_sanger_be.banner import render


def test_banner_fits_the_terminal_and_shows_the_version():
    for columns in (80, 200):
        text = render(columns)
        assert f"v{__version__}" in text
        assert max(len(line) for line in text.splitlines()) < columns
    assert render(200).count("(_\\__/(,_") == 4
    assert render(80).count("(_\\__/(,_") == 1


def test_help_shows_the_banner_in_a_terminal(monkeypatch, capsys):
    import sys

    import pytest

    from fast_sanger_be.cli import parse_args

    monkeypatch.setattr(sys.stdout, "isatty", lambda: True, raising=False)
    with pytest.raises(SystemExit):
        parse_args(["-h"])
    out = capsys.readouterr().out
    assert "Sanger Base Editing Pipeline" in out and "usage: fast-sanger-be" in out
