"""
``./sbin/ucsc-cytoband-to-json cytoband-hg38.txt.gz | gzip -c >bioutils/_data/cytobands/ucsc-hg38.json.gz``
"""

import gzip
import json

import pkg_resources

_data_dir = "_data/cytobands"
_data_path_fmt = _data_dir + "/" + "{name}.json.gz"


def get_cytoband_names():
    """Retrieves available cytobands from the ``_data/cytobands`` directory.

    Returns:
        list of str: The names of the available cytobands.

    Examples:
        >>> sorted(get_cytoband_names())
        ['ucsc-hg19', 'ucsc-hg38']
    """

    return [
        n.replace(".json.gz", "")
        for n in pkg_resources.resource_listdir(__name__, _data_dir)
        if n.endswith(".json.gz")
    ]


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

    fn = pkg_resources.resource_filename(__name__, _data_path_fmt.format(name=name))
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
