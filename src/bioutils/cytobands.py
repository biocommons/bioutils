"""
``./sbin/ucsc-cytoband-to-json cytoband-hg38.txt.gz | gzip -c >bioutils/_data/cytobands/ucsc-hg38.json.gz``
"""

import gzip
import json
from importlib import resources
from pathlib import Path

_data_dir = Path(str(resources.files(__package__) / "_data" / "cytobands"))


def get_cytoband_names():
    """Retrieves available cytobands from the ``_data/cytobands`` directory.

    Returns:
        list of str: The names of the available cytobands.

    Examples:
        >>> sorted(get_cytoband_names())
        ['ucsc-hg19', 'ucsc-hg38']
    """

    return [n.name.replace(".json.gz", "") for n in _data_dir.glob("*.json.gz")]


def get_cytoband_map(name):
    """Retrives a cytoband by name.

    Args:
        name (str): The name of the cytoband to retrieve.

    Returns:
        dict: A dictionary of the cytoband data.


    Examples:
        >>> map = get_cytoband_map("ucsc-hg38")
        >>> map["1"]["p32.2"]
        [55600000, 58500000, 'gpos50']
    """

    fn = _data_dir / f"{name}.json.gz"
    return json.load(gzip.open(fn, mode="rt", encoding="utf-8"))


def get_cytoband_maps(names=[]):
    """Retrieves data from multiple cytobands.

    If cytobands are not specified, retrieves data from all available ones.

    Args:
        names (list of str, optional): The names of cytobands to retrieve data for.

    Returns:
        dict: A dictionary of the form ``{cytoband_name, cytoband_data}``.

    Examples:
        >>> maps = get_cytoband_maps()
        >>> maps["ucsc-hg38"]["1"]["p32.2"]
        [55600000, 58500000, 'gpos50']
        >>> maps["ucsc-hg19"]["1"]["p32.2"]
        [56100000, 59000000, 'gpos50']
    """

    if names == []:
        names = get_cytoband_names()
    return {name: get_cytoband_map(name) for name in names}
