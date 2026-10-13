"""Pseudoautosomal region (PAR) coordinates for human reference assemblies.

The pseudoautosomal regions are the segments at either end of X and Y that
recombine and are diploid in both sexes. Code that treats X or Y as hemizygous
in males, such as sex inference, copy-number estimation or variant calling
with per-sample ploidy, has to mask them first.

Coordinates are **interbase** (0-based, left-closed, right-open), the same
convention as :mod:`bioutils.cytobands` and the ``_i`` system described in
:mod:`bioutils.coordinates`. NCBI publishes the regions 1-based inclusive; the
generator converts ``[start, end]`` to ``[start - 1, end]``, so a
1-based position ``pos`` (for example a VCF ``POS``) is inside a region
``[s, e]`` when ``s <= pos - 1 < e``. The original NCBI values are recoverable
as ``[start + 1, end]`` and are recorded in ``tests/test_par.py``.

Positions are looked up by RefSeq sequence accession (``"NC_000023.11"``), as
elsewhere in biocommons. A versioned accession identifies one sequence of one
build, so no assembly name is needed and there is no ambiguity over ``chr``
prefixes. The whole-assembly views (:func:`get_par_map`, :func:`get_par_maps`)
are keyed by the assembly's own chromosome names (``"X"``, ``"Y"``), and
:func:`get_par_meta` gives the accession of each.

Region labels are ``"PAR1"`` (the p-terminal region, shared by the short arms)
and ``"PAR2"`` (the q-terminal region): NCBI's own "PAR#1" and "PAR#2" names
without the ``#``, which gives the plain forms used in the literature and by
other resources.

Data is generated from each assembly's ``assembly_regions.txt`` report,
published by NCBI beside the assembly report, whose ``PAR`` rows give the
regions::

  ./sbin/ncbi-par-to-json GRCh38.p14 GCF_000001405.40_GRCh38.p14_assembly_regions.txt \\
      | gzip -c >src/bioutils/_data/par/GRCh38.json.gz

PAR coordinates are stable across patch releases within a build, so files are
named for the build (``"GRCh38"``). Functions take an assembly name as used by
:func:`bioutils.assemblies.get_assembly`, either a build or a patch release
(``"GRCh38.p14"``), which resolves to its build. Each file is a flat JSON document in the style of
``_data/assemblies``: provenance fields (``assembly``, ``refseq_ac``,
``accessions``, ``coordinates``, ``source``, ``source_url``) beside the
``regions`` data (see :func:`get_par_meta`).
"""

import functools
import gzip
import json
import re
from importlib import resources
from pathlib import Path
from typing import TypedDict, cast

from bioutils.assemblies import get_assemblies

_data_dir = Path(str(resources.files("bioutils") / "_data" / "par"))


class ParMeta(TypedDict):
    """Provenance fields of a PAR dataset; see :func:`get_par_meta`."""

    assembly: str
    refseq_ac: str
    accessions: dict[str, str]
    coordinates: str
    source: str
    source_url: str


def _build_name(assembly_name: str) -> str:
    """Maps an assembly name to the build its PAR data is stored under."""
    build = re.sub(r"\.p\d+$", "", assembly_name)
    if build not in get_par_names():
        msg = f"no PAR data for assembly {assembly_name!r}; available builds: {sorted(get_par_names())}"
        raise ValueError(msg)
    return build


@functools.cache
def _load(assembly_name: str) -> dict:
    fn = _data_dir / f"{_build_name(assembly_name)}.json.gz"
    with gzip.open(fn, mode="rt", encoding="utf-8") as fh:
        return json.load(fh)


def get_par_names() -> list[str]:
    """Retrieves available PAR datasets from the ``_data/par`` directory.

    Returns:
        list of str: The builds with PAR data. Any patch release of these
        (e.g. ``"GRCh38.p14"``) is also accepted by the other functions.

    Examples:
        >>> sorted(get_par_names())
        ['GRCh37', 'GRCh38', 'T2T-CHM13v2.0']
    """

    return [n.name.replace(".json.gz", "") for n in _data_dir.glob("*.json.gz")]


def get_par_map(assembly_name: str) -> dict[str, dict[str, list[int]]]:
    """Retrieves PAR coordinates for an assembly.

    Args:
        assembly_name (str): A build (``"GRCh38"``) or patch release
            (``"GRCh38.p14"``); patches share their build's coordinates.

    Returns:
        dict: A dictionary of the form ``{chromosome: {label: [start_i, end_i]}}``,
        with interbase (0-based, half-open) coordinates.

    Examples:
        >>> par = get_par_map("GRCh38")
        >>> par["X"]["PAR1"]
        [10000, 2781479]
        >>> par["Y"]["PAR2"]
        [56887902, 57217415]

        A patch release gives its build's coordinates:

        >>> get_par_map("GRCh38.p14") == par
        True

        PAR1 is not at the same X and Y coordinates in every build:

        >>> get_par_map("GRCh37")["X"]["PAR1"]
        [60000, 2699520]
        >>> get_par_map("GRCh37")["Y"]["PAR1"]
        [10000, 2649520]
    """

    return _load(assembly_name)["regions"]


def get_par_meta(assembly_name: str) -> ParMeta:
    """Retrieves provenance metadata for a PAR dataset.

    Args:
        assembly_name (str): A build or patch release, as for :func:`get_par_map`.

    Returns:
        dict: Every top-level field of the dataset except ``regions``: the
        assembly the data was generated from and its RefSeq accession, the X and Y sequence accessions,
        the coordinate convention, and the upstream source and URL.

    Examples:
        >>> meta = get_par_meta("GRCh38")
        >>> meta["assembly"]
        'GRCh38.p14'
        >>> meta["refseq_ac"]
        'GCF_000001405.40'
        >>> meta["accessions"]["X"]
        'NC_000023.11'
        >>> meta["coordinates"]
        'interbase (0-based, left-closed, right-open)'
    """

    return cast("ParMeta", {k: v for k, v in _load(assembly_name).items() if k != "regions"})


def get_par_maps(
    assembly_names: list[str] | None = None,
) -> dict[str, dict[str, dict[str, list[int]]]]:
    """Retrieves PAR coordinates for multiple datasets.

    If datasets are not specified, retrieves data for all available ones.

    Args:
        assembly_names (list of str, optional): Builds or patch releases, as
            for :func:`get_par_map`.

    Returns:
        dict: A dictionary of the form ``{assembly_name: par_data}``.

    Examples:
        >>> maps = get_par_maps()
        >>> maps["GRCh38"]["X"]["PAR1"]
        [10000, 2781479]
    """

    if assembly_names is None:
        assembly_names = get_par_names()
    return {name: get_par_map(name) for name in assembly_names}


@functools.cache
def _ac_regions() -> dict[str, dict[str, list[int]]]:
    """Maps each X and Y accession with PAR data to its regions."""
    ac_regions = {}
    for name in get_par_names():
        doc = _load(name)
        for chrom, ac in doc["accessions"].items():
            ac_regions[ac] = doc["regions"][chrom]
    return ac_regions


@functools.cache
def _sex_chromosome_acs() -> frozenset[str]:
    """X and Y accessions of every assembly in bioutils, with or without PAR data."""
    return frozenset(
        s["refseq_ac"]
        for assembly in get_assemblies().values()
        for s in assembly["sequences"]
        if s["name"] in ("X", "Y")
    )


def get_par_regions(ac: str) -> dict[str, list[int]]:
    """Retrieves the pseudoautosomal regions of a sequence.

    Args:
        ac (str): Versioned RefSeq accession of the sequence, e.g.
            ``"NC_000023.11"`` (GRCh38 X).

    Returns:
        dict: A dictionary of the form ``{label: [start_i, end_i]}``, with
        interbase coordinates. Empty for a sequence with no pseudoautosomal
        regions, such as an autosome.

    Raises:
        ValueError: If ``ac`` is not a versioned accession (e.g. ``"X"``,
            ``"chrX"`` or ``"NC_000023"``), or is the X or Y of an assembly
            with no PAR data, such as NCBI36.

    Examples:
        >>> get_par_regions("NC_000023.11")
        {'PAR1': [10000, 2781479], 'PAR2': [155701382, 156030895]}
        >>> get_par_regions("NC_000001.11")
        {}
    """

    if not re.fullmatch(r"[A-Z]{2}_\d+\.\d+", ac):
        msg = f"expected a versioned sequence accession such as 'NC_000023.11', got {ac!r}"
        raise ValueError(msg)
    if ac in _ac_regions():
        return _ac_regions()[ac]
    if ac in _sex_chromosome_acs():
        msg = f"no PAR data for sequence {ac!r}; available builds: {sorted(get_par_names())}"
        raise ValueError(msg)
    return {}


def get_par_label(ac: str, pos_i: int) -> str | None:
    """Names the pseudoautosomal region containing a position, if any.

    Args:
        ac (str): Versioned RefSeq accession of the sequence, as for
            :func:`get_par_regions`.
        pos_i (int): Interbase (0-based) position of a nucleotide, i.e. the
            interval ``[pos_i, pos_i + 1)``. For a 1-based position such as a
            VCF ``POS``, pass ``POS - 1``.

    Returns:
        str or None: ``"PAR1"`` or ``"PAR2"`` when the position falls inside a
        region, otherwise None. None for any sequence that has no
        pseudoautosomal regions.

    Raises:
        ValueError: As for :func:`get_par_regions`.

    Examples:
        >>> get_par_label("NC_000023.11", 1000000)
        'PAR1'
        >>> get_par_label("NC_000023.11", 155800000)
        'PAR2'

        The non-PAR bulk of X, where males are hemizygous:

        >>> get_par_label("NC_000023.11", 100000000) is None
        True

        Autosomes never are:

        >>> get_par_label("NC_000001.11", 1000000) is None
        True
    """

    for label, (start_i, end_i) in get_par_regions(ac).items():
        if start_i <= pos_i < end_i:
            return label
    return None


def in_par(ac: str, pos_i: int) -> bool:
    """Tests whether a position falls in a pseudoautosomal region.

    Args:
        ac (str): Versioned RefSeq accession of the sequence, as for
            :func:`get_par_regions`.
        pos_i (int): Interbase (0-based) position of a nucleotide. For a
            1-based position such as a VCF ``POS``, pass ``POS - 1``.

    Returns:
        bool: True if the position is within a PAR on ``ac``. Always False
        for sequences that have no pseudoautosomal regions.

    Raises:
        ValueError: As for :func:`get_par_regions`.

    Examples:
        The first and last bases of GRCh38 PAR1 are NCBI 10,001 and 2,781,479
        (1-based), so in interbase terms:

        >>> in_par("NC_000023.11", 10000)
        True
        >>> in_par("NC_000023.11", 9999)
        False
        >>> in_par("NC_000023.11", 2781478)
        True
        >>> in_par("NC_000023.11", 2781479)
        False

        From a VCF record:

        >>> vcf_pos = 10001
        >>> in_par("NC_000023.11", vcf_pos - 1)
        True
    """

    return get_par_label(ac, pos_i) is not None
