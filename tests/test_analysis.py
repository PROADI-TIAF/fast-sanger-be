import numpy as np
import pytest
from Bio.Seq import reverse_complement

from fast_sanger_be.analysis import assess_quality, editing_values, locate_guide, orient, quantify
from fast_sanger_be.chromatogram import Chromatogram
from fast_sanger_be.config import Parameters
from fast_sanger_be.simulate import edited_composition, reverse_composition, simulate_read

GUIDE = "CTTGACACAACTGTGTTCAC"
A_POSITIONS = [i for i, b in enumerate(GUIDE) if b == "A"]


def _reference(rng):
    flank = lambda n: "".join(rng.choice(list("ACGT"), n))  # noqa: E731
    return flank(120) + GUIDE + "TGG" + flank(120)


def _chrom(sim, name="x"):
    return Chromatogram(name, sim.bases, np.array(sim.locations), sim.traces, np.array(sim.quality, float))


def test_locate_guide_tolerates_edits_and_strand():
    read = "TTTT" + GUIDE.replace("A", "G") + "GGGG"
    hit = locate_guide(read, GUIDE, "ABE", 3)
    assert hit.strand == "+" and hit.start == 4 and hit.mismatches == 0
    hit = locate_guide(reverse_complement(read), GUIDE, "ABE", 3)
    assert hit.strand == "-" and hit.start == 4
    assert locate_guide(read, GUIDE, "CBE", 3) is None  # A>G changes are mismatches for a CBE


@pytest.mark.parametrize("level", [0.0, 0.1, 0.35, 0.7, 0.95])
@pytest.mark.parametrize("read", ["F", "R"])
def test_quantification_recovers_simulated_editing(level, read):
    rng = np.random.default_rng(int(level * 100) + (read == "R"))
    reference = _reference(rng)
    efficiencies = {p: level for p in A_POSITIONS}
    composition = edited_composition(reference, 120, GUIDE, "ABE", efficiencies)
    if read == "R":
        composition = reverse_composition(composition)
    chrom = _chrom(simulate_read(composition, rng))
    params = Parameters()
    hit = locate_guide(chrom.bases, GUIDE, "ABE", params.max_guide_mismatches)
    assert hit is not None and hit.strand == ("+" if read == "F" else "-")
    oriented = orient(chrom, hit)
    qc = assess_quality(oriented, hit.start, GUIDE, "ABE", params)
    assert qc.passed, qc.reasons
    values = editing_values(quantify(oriented, hit.start, len(GUIDE), params), GUIDE, "ABE")
    errors = [abs(v - 100 * level) for v in values.values()]
    # Simulated cross-talk varies 1.5-4.5% per peak, so single positions can be off by a few points.
    assert np.mean(errors) < 2.5 and max(errors) < 5.0, values


def test_qc_rejects_mixed_template():
    rng = np.random.default_rng(3)
    reference = _reference(rng)
    composition = [{b: 1.0} for b in reference]
    contaminant = composition[9:] + composition[:9]
    chrom = _chrom(simulate_read(composition, rng, contaminant=contaminant, contaminant_fraction=0.45))
    params = Parameters()
    hit = locate_guide(chrom.bases, GUIDE, "ABE", 6)
    if hit is None:
        return  # guide not even found: also excluded
    qc = assess_quality(orient(chrom, hit), hit.start, GUIDE, "ABE", params)
    assert not qc.passed and any("overlapping peaks" in r for r in qc.reasons)


def test_qc_rejects_guide_near_read_end():
    rng = np.random.default_rng(5)
    reference = "".join(rng.choice(list("ACGT"), 5)) + GUIDE + "".join(rng.choice(list("ACGT"), 80))
    chrom = _chrom(simulate_read([{b: 1.0} for b in reference], rng))
    hit = locate_guide(chrom.bases, GUIDE, "ABE", 3)
    qc = assess_quality(orient(chrom, hit), hit.start, GUIDE, "ABE", Parameters())
    assert not qc.passed and "does not cover" in qc.reasons[0]


def test_mixed_base_calls_do_not_create_false_sites():
    # Degraded read ends from mixed-base calling are full of IUPAC codes, which must not
    # match the guide; only R (A/G) at ABE targets is accepted.
    guide = "AGATATTTGCATTGAGATAG"
    real_site = reverse_complement("AGATRTTTGCRTTGAGATAG")
    read = "ACGTTGCA" * 10 + real_site + "TGCA" * 20 + "RRAATTTTTSKKKKKKKTSS" + "KSWMRY" * 10
    hit = locate_guide(read, guide, "ABE", 3)
    assert hit.strand == "-" and hit.mismatches == 0 and not hit.ambiguous
    assert locate_guide("TTTT" + "KKKKSSSSRRRRYYYYWWWW" + "TTTT", "GGGGCCCCAAAATTTTAAAA", "ABE", 3) is None


def _blob_read():
    """40% edited read with overlapping peaks over the 10 nt before the guide and its first 8 nt."""
    rng = np.random.default_rng(11)
    reference = _reference(rng)
    composition = edited_composition(reference, 120, GUIDE, "ABE", {p: 0.4 for p in A_POSITIONS})
    shifted = composition[6:] + composition[:6]
    mixed = [dict(c) for c in composition]
    for i in range(110, 128):
        mixed[i] = {b: 0.5 * f for b, f in composition[i].items()}
        for b, f in shifted[i].items():
            mixed[i][b] = mixed[i].get(b, 0.0) + 0.5 * f
    chrom = _chrom(simulate_read(mixed, rng))
    hit = locate_guide(chrom.bases, GUIDE, "ABE", 3)
    assert hit is not None
    return orient(chrom, hit), hit.start


def test_local_noise_rejects_the_whole_read_by_default():
    oriented, start = _blob_read()
    qc = assess_quality(oriented, start, GUIDE, "ABE", Parameters())
    assert not qc.passed
    skipped = qc.unreliable_labels(GUIDE, "ABE")
    assert skipped and any(", ".join(skipped) in reason for reason in qc.reasons), qc.reasons


def test_partial_reads_keep_the_clean_positions():
    # With partial_reads, a dye blob over the 5' half of the protospacer leaves the 3' half usable.
    oriented, start = _blob_read()
    params = Parameters(partial_reads=True)
    qc = assess_quality(oriented, start, GUIDE, "ABE", params)
    assert qc.passed, qc.reasons

    skipped = qc.unreliable_labels(GUIDE, "ABE")
    assert skipped, "positions under the blob should not be quantified"
    values = editing_values(quantify(oriented, start, len(GUIDE), params), GUIDE, "ABE", qc.reliable)
    assert all(np.isnan(values[label]) for label in skipped)
    kept = {label: v for label, v in values.items() if label not in skipped}
    assert kept, "the clean half must still be quantified"
    assert max(abs(v - 40.0) for v in kept.values()) < 6.0, kept


def test_partial_read_rejected_when_nothing_is_quantifiable():
    rng = np.random.default_rng(13)
    reference = _reference(rng)
    composition = [{b: 1.0} for b in reference]
    chrom = _chrom(simulate_read(composition, rng, contaminant=composition[9:] + composition[:9],
                                 contaminant_fraction=0.5))
    params = Parameters(partial_reads=True)
    hit = locate_guide(chrom.bases, GUIDE, "ABE", 6)
    if hit is None:
        return
    qc = assess_quality(orient(chrom, hit), hit.start, GUIDE, "ABE", params)
    assert not qc.passed and qc.unreliable_labels(GUIDE, "ABE")
