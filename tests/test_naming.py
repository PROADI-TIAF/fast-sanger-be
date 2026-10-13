from pathlib import Path

from fast_sanger_be.naming import (
    find_conflicts,
    guess_fields,
    parse_filename,
    parse_guides,
    parse_rename_table,
    parse_sample_files,
)


def test_parse_single_and_multiple_guides():
    f, err = parse_filename(Path("Clone-3_sgHBB_ABE_F.ab1"))
    assert err is None and (f.sample, f.guides, f.editor, f.read) == ("Clone-3", ("sgHBB",), "ABE", "F")
    f, err = parse_filename(Path("S1_g1_g2-b_cbe_r.ab1"))
    assert err is None and f.guides == ("g1", "g2-b") and f.editor == "CBE" and f.read == "R"


def test_parse_rejects_bad_names():
    for name, fragment in [
        ("S1_g1_F.ab1", "expected"),
        ("S1_g1_XBE_F.ab1", "ABE or CBE"),
        ("S1_g1_ABE_X.ab1", "F or R"),
        ("S 1_g1_ABE_F.ab1", "letters, digits"),
        ("S1_g1_g1_ABE_F.ab1", "twice"),
    ]:
        f, err = parse_filename(Path(name))
        assert f is None and fragment in err, name


def test_conflicts_detected_per_guide():
    files, _ = parse_sample_files([Path("S1_g1_ABE_F.ab1"), Path("S1_g1_g2_ABE_F.ab1"), Path("S1_g2_CBE_F.ab1")])
    conflicts = find_conflicts(files)
    assert len(conflicts) == 1 and "guide g1" in conflicts[0]


def test_parse_guides_fasta_and_brackets():
    guides, errors = parse_guides(">g1\nACGTACGTACGTACGTACGT\n[g2]\nacgtacgtac gtacgtacgt\n")
    assert not errors and guides["g2"].sequence == "ACGTACGTACGTACGTACGT"


def test_parse_guides_errors():
    _, errors = parse_guides(">g1\nACGTNACGT\n>g1\nACGTACGTACGTACGTACGT\n>bad_name\nACGTACGTACGTACGTACGT\n>g3\nACGT\n")
    text = " ".join(errors)
    assert "A, C, G and T" in text and "defined twice" not in text  # first g1 is invalid, second is accepted
    assert "letters, digits" in text and "15-30" in text


def test_rename_table_validation(tmp_path):
    invalid = [tmp_path / "a.ab1", tmp_path / "b.ab1"]
    header = "current_file\tinclude\tsample\tguides\teditor\tread\n"
    ok = header + "a.ab1\tyes\tS1\tg1, g2\tABE\tF\nb.ab1\tno\t\t\t\t\n"
    plan, errors = parse_rename_table(ok, tmp_path, invalid, set())
    assert not errors and plan.renames[0][1].name == "S1_g1_g2_ABE_F.ab1" and plan.excluded == [invalid[1]]

    clash = header + "a.ab1\tyes\tS1\tg1\tABE\tF\nb.ab1\tyes\tS1\tg1\tABE\tF\n"
    plan, errors = parse_rename_table(clash, tmp_path, invalid, set())
    assert plan is None and any("also used" in e for e in errors)

    plan, errors = parse_rename_table(header + "a.ab1\tyes\tS1\tg1\tABE\tF\n", tmp_path, invalid, {"S1_g1_ABE_F.ab1"})
    assert any("already exists" in e for e in errors) and any("missing from the table" in e for e in errors)


def test_guess_fields():
    g = guess_fields("Clone7-sgHBB-ABE-fwd.ab1", {"sgHBB"})
    assert g == {"sample": "Clone7", "guides": "sgHBB", "editor": "ABE", "read": "F"}
