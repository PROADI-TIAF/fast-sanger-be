import json
import shutil
import zipfile
from pathlib import Path

import pytest

from fast_sanger_be.simulate import make_example
from fast_sanger_be.web import WebInputError, parameters_from_dict, run, zip_folder


def test_parameters_from_dict_casts_and_ignores_unknown_keys():
    params = parameters_from_dict({"qc_flank": "12", "partial_reads": True, "mixed_peak_ratio": "0.3", "nope": 1})
    assert params.qc_flank == 12 and params.partial_reads is True and params.mixed_peak_ratio == 0.3


def test_run_returns_summary_and_writes_results(tmp_path):
    data = tmp_path / "input"
    make_example(data)
    for extra in data.iterdir():
        if extra.suffix not in (".ab1", ".txt") or extra.name == "truth.tsv":
            extra.unlink() if extra.is_file() else shutil.rmtree(extra)
    seen = []
    summary = json.loads(run(str(data), str(tmp_path / "Final_Results"), json.dumps({"figure_dpi": 60}),
                             json.dumps({}), on_log=lambda level, msg: seen.append(msg)))
    assert summary["counts"] and summary["editors"]
    block = next(iter(summary["editors"].values()))[0]
    assert block["labels"] and block["samples"] and block["reads"]
    for read in block["reads"]:
        if read["figure"]:
            assert (tmp_path / "Final_Results" / read["figure"]).exists()
    assert any("Step 3/4" in m for m in seen)
    assert (tmp_path / "Final_Results" / "run_log.txt").read_text()
    archive = zip_folder(str(tmp_path / "Final_Results"), str(tmp_path / "out.zip"))
    names = zipfile.ZipFile(archive).namelist()
    assert "Final_Results/run_log.txt" in names and "Final_Results/input_files.tsv" in names


def test_run_rejects_missing_guide(tmp_path):
    data = tmp_path / "input"
    make_example(data)
    (data / "guides.txt").write_text(">other\nACGTACGTACGTACGTACGT\n")
    with pytest.raises(WebInputError, match="no sequence for guide"):
        run(str(data), str(tmp_path / "out"))
