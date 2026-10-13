"""Tests for pseudoautosomal region coordinates.

``NCBI`` below is transcribed from each build's ``assembly_regions.txt`` exactly
as NCBI publishes it (1-based inclusive), so that a change to the shipped data has
to be justified against the published source. The shipped data is interbase
(0-based, half-open); the conversion is ``[start - 1, end]``.
"""

import pytest

from bioutils.assemblies import get_assembly
from bioutils.par import (
    get_par_label,
    get_par_map,
    get_par_maps,
    get_par_meta,
    get_par_names,
    get_par_regions,
    in_par,
)

# The PAR rows of NCBI's assembly_regions.txt, 1-based inclusive, as published
# ("PAR#1" / "PAR#2"). The T2T values match the consortium's chm13v2.0_PAR.bed.
NCBI = {
    "GRCh37": {
        "X": {"PAR1": [60001, 2699520], "PAR2": [154931044, 155260560]},
        "Y": {"PAR1": [10001, 2649520], "PAR2": [59034050, 59363566]},
    },
    "GRCh38": {
        "X": {"PAR1": [10001, 2781479], "PAR2": [155701383, 156030895]},
        "Y": {"PAR1": [10001, 2781479], "PAR2": [56887903, 57217415]},
    },
    "T2T-CHM13v2.0": {
        "X": {"PAR1": [1, 2394410], "PAR2": [153925835, 154259566]},
        "Y": {"PAR1": [1, 2458320], "PAR2": [62122810, 62460029]},
    },
}

# RefSeq accessions of each build's X and Y, from the assembly reports.
ACCESSIONS = {
    "GRCh37": {"X": "NC_000023.10", "Y": "NC_000024.9"},
    "GRCh38": {"X": "NC_000023.11", "Y": "NC_000024.10"},
    "T2T-CHM13v2.0": {"X": "NC_060947.1", "Y": "NC_060948.1"},
}
GRCH37_X, GRCH37_Y = ACCESSIONS["GRCh37"]["X"], ACCESSIONS["GRCh37"]["Y"]
GRCH38_X = ACCESSIONS["GRCh38"]["X"]
T2T_X = ACCESSIONS["T2T-CHM13v2.0"]["X"]

# The same regions in interbase coordinates, as shipped.
EXPECTED = {
    name: {
        chrom: {par: [start - 1, end] for par, (start, end) in pars.items()}
        for chrom, pars in chroms.items()
    }
    for name, chroms in NCBI.items()
}


def test_get_par_names():
    assert sorted(get_par_names()) == ["GRCh37", "GRCh38", "T2T-CHM13v2.0"]


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_coordinates_match_converted_source(name):
    assert get_par_map(name) == EXPECTED[name]


@pytest.mark.parametrize("name", sorted(NCBI))
def test_ncbi_values_recoverable(name):
    # [start_i + 1, end_i] gives back exactly what NCBI publishes.
    recovered = {
        chrom: {par: [s + 1, e] for par, (s, e) in pars.items()}
        for chrom, pars in get_par_map(name).items()
    }
    assert recovered == NCBI[name]


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_meta_records_provenance(name):
    meta = get_par_meta(name)
    assembly = get_assembly(meta["assembly"])
    assert set(meta) == {
        "assembly",
        "refseq_ac",
        "accessions",
        "coordinates",
        "source",
        "source_url",
    }
    assert "regions" not in meta
    assert meta["assembly"].startswith(name)
    assert meta["refseq_ac"] == assembly["refseq_ac"]
    assert set(meta["accessions"]) == {"X", "Y"}
    for chrom, ac in meta["accessions"].items():
        assert any(s["name"] == chrom and s["refseq_ac"] == ac for s in assembly["sequences"])
    assert meta["coordinates"] == "interbase (0-based, left-closed, right-open)"
    assert meta["source_url"].endswith("_assembly_regions.txt")
    assert meta["assembly"] in meta["source_url"]
    assert meta["refseq_ac"] in meta["source_url"]


@pytest.mark.parametrize(
    ("assembly_name", "build"),
    [
        ("GRCh37.p13", "GRCh37"),
        ("GRCh37.p2", "GRCh37"),
        ("GRCh38.p14", "GRCh38"),
        ("GRCh38.p1", "GRCh38"),
    ],
)
def test_patch_release_resolves_to_build(assembly_name, build):
    # PARs are stable across patches, so any patch gets its build's data.
    assert get_par_map(assembly_name) == EXPECTED[build]
    assert get_par_meta(assembly_name) == get_par_meta(build)
    assert get_par_meta(assembly_name)["accessions"] == ACCESSIONS[build]


@pytest.mark.parametrize("assembly_name", ["hg38", "GRCh39", "GRCh38.p14x", "NCBI36"])
def test_unknown_assembly_is_rejected(assembly_name):
    with pytest.raises(ValueError, match="no PAR data for assembly"):
        get_par_map(assembly_name)


def test_get_par_maps_returns_all():
    assert get_par_maps() == EXPECTED


def test_get_par_maps_subset():
    assert get_par_maps(["GRCh38"]) == {"GRCh38": EXPECTED["GRCh38"]}


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_meta_accessions_match_assembly(name):
    assert get_par_meta(name)["accessions"] == ACCESSIONS[name]


@pytest.mark.parametrize("name", sorted(EXPECTED))
@pytest.mark.parametrize("chrom", ["X", "Y"])
def test_get_par_regions_by_accession(name, chrom):
    assert get_par_regions(ACCESSIONS[name][chrom]) == EXPECTED[name][chrom]


@pytest.mark.parametrize("name", sorted(NCBI))
@pytest.mark.parametrize("chrom", ["X", "Y"])
@pytest.mark.parametrize("par", ["PAR1", "PAR2"])
def test_boundaries_half_open(name, chrom, par):
    ac = ACCESSIONS[name][chrom]
    start_1, end_1 = NCBI[name][chrom][par]
    # First and last NCBI bases, as interbase positions, are inside.
    assert get_par_label(ac, start_1 - 1) == par
    assert get_par_label(ac, end_1 - 1) == par
    # One before the first base and the position at ``end`` are outside.
    assert get_par_label(ac, start_1 - 2) is None
    assert get_par_label(ac, end_1) is None


@pytest.mark.parametrize("name", sorted(NCBI))
@pytest.mark.parametrize("chrom", ["X", "Y"])
def test_in_par_agrees_with_get_par_label(name, chrom):
    ac = ACCESSIONS[name][chrom]
    for pos_i in (0, 10000, 1000000, 100000000, 155000000, 156100000):
        assert in_par(ac, pos_i) is (get_par_label(ac, pos_i) is not None)


def test_non_par_x_is_excluded():
    # The hemizygous bulk of X in males, between PAR1 and PAR2.
    assert in_par(GRCH38_X, 100000000) is False


@pytest.mark.parametrize("ac", ["NC_000001.11", "NC_000022.10", "NC_012920.1", "NC_060925.1"])
def test_autosomes_are_never_par(ac):
    # GRCh38 1, GRCh37 22, MT and T2T 1.
    assert get_par_regions(ac) == {}
    assert in_par(ac, 1000000) is False
    assert get_par_label(ac, 1000000) is None


@pytest.mark.parametrize("ac", ["X", "chrX", "chrx", "NC_000023", "nc_000023.11", "NC_000023.11 "])
def test_non_accession_is_rejected(ac):
    # A chromosome name or unversioned accession does not identify a sequence,
    # so it is a caller error rather than something to silently treat as non-PAR.
    with pytest.raises(ValueError, match="versioned sequence accession"):
        in_par(ac, 1000000)


@pytest.mark.parametrize("ac", ["NC_000023.9", "NC_000024.7"])
def test_sex_chromosome_without_par_data_is_rejected(ac):
    # NCBI36 X and Y: known sex chromosomes, but no PAR data shipped, so
    # answering False would be wrong rather than merely unknown.
    with pytest.raises(ValueError, match="no PAR data for sequence"):
        in_par(ac, 1000000)


def test_vcf_position_conversion():
    # A VCF POS is 1-based; the first PAR1 base in GRCh38 is POS 10001.
    assert in_par(GRCH38_X, 10001 - 1) is True
    assert in_par(GRCH38_X, 10000 - 1) is False


def test_par1_differs_between_x_and_y_in_grch37():
    # GRCh37 PAR1 does not start at the same coordinate on X and Y, so callers
    # must not assume the two are interchangeable.
    par = get_par_map("GRCh37")
    assert par["X"]["PAR1"] != par["Y"]["PAR1"]
    assert in_par(GRCH37_Y, 10000) is True
    assert in_par(GRCH37_X, 10000) is False


def test_par2_moved_between_builds():
    # GRCh38 relocated PAR2 on both chromosomes; the two builds' PAR2
    # intervals do not overlap at all.
    x37 = get_par_map("GRCh37")["X"]["PAR2"]
    x38 = get_par_map("GRCh38")["X"]["PAR2"]
    assert x37[1] <= x38[0]
    y37 = get_par_map("GRCh37")["Y"]["PAR2"]
    y38 = get_par_map("GRCh38")["Y"]["PAR2"]
    assert y38[1] <= y37[0]


def test_t2t_par1_starts_at_first_base():
    # T2T-CHM13v2.0 has no telomeric gap, so PAR1 begins at the first base of
    # both X and Y, unlike GRCh37/38 where it starts after an N run.
    par = get_par_map("T2T-CHM13v2.0")
    assert par["X"]["PAR1"][0] == 0
    assert par["Y"]["PAR1"][0] == 0
    assert in_par(T2T_X, 0) is True
    assert in_par(GRCH38_X, 0) is False


def test_t2t_par1_lengths_differ_between_x_and_y():
    # The T2T Y comes from a different individual (HG002) than the X (CHM13),
    # so PAR1 is not the same length on the two chromosomes.
    par = get_par_map("T2T-CHM13v2.0")
    assert par["X"]["PAR1"][1] != par["Y"]["PAR1"][1]


@pytest.mark.parametrize("chrom", ["X", "Y"])
def test_t2t_par2_runs_to_chromosome_end(chrom):
    # With no telomeric gap, PAR2 ends at the last base of the chromosome.
    length = next(
        s["length"] for s in get_assembly("T2T-CHM13v2.0")["sequences"] if s["name"] == chrom
    )
    assert get_par_map("T2T-CHM13v2.0")[chrom]["PAR2"][1] == length
