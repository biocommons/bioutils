"""Creates dictionaries of genome assembly data as provided by  

ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/*.assembly.txt

Assemblies are stored in json files with the package in
``_data/assemblies/``. Those files are built with sbin/assembly-to-json,
also in this package.

Definitions:

* accession ``ac``: symbol used to refer to a sequence (e.g., NC_000001.10)
* name: human-label (e.g., '1', 'MT', 'HSCHR6_MHC_APD_CTG1') that
  refers to a sequence, unique within some domain (e.g., GRCh37.p10)
* chromosome (chr): subset of names that refer to chromosomes 1..22, X, Y, MT
* aliases: list of other names; uniqueness unknown

.. note:: Some users prefer using a 'chr' prefix for chromosomes and
    some don't.  Some prefer upper case and others prefer lower.  This
    rift is unfortunate and creates unnecessary friction in sharing data.
    You say TO-my-to and I say TO-mah-to doesn't apply here.  This code
    favors using the authoritative names exactly as defined in the
    assembly records.  Users are encouraged to use sequence names
    verbatim, without prefixes or case changes.
"""

import gzip
import json

import pkg_resources

_assy_dir = "_data/assemblies"
_assy_path_fmt = _assy_dir + "/" + "{name}.json.gz"


def get_assembly_names():
    """Retrieves available assemblies from the ``_data/assemblies`` directory.

    Returns:
        list of str: The names of the available assemblies.

    Examples:
        >>> assy_names = get_assembly_names()

        >>> 'GRCh37.p13' in assy_names
        True
    """

    return [
        n.replace(".json.gz", "")
        for n in pkg_resources.resource_listdir(__name__, _assy_dir)
        if n.endswith(".json.gz")
    ]


def get_assembly(name):
    """Retreives the assembly data for a given assembly.

    Args:
        name (str): The name of the assembly to retrieve data for.

    Returns:
        dict: A dictionary of the assembly data. See examples for details.


    Examples:
        >>> assy = get_assembly('GRCh37.p13')

        >>> assy['name']
        'GRCh37.p13'

        >>> assy['description']
        'Genome Reference Consortium Human Build 37 patch release 13 (GRCh37.p13)'

        >>> assy['refseq_ac']
        'GCF_000001405.25'

        >>> assy['genbank_ac']
        'GCA_000001405.14'

        >>> len(assy['sequences'])
        297

        >>> import pprint
        >>> pprint.pprint(assy['sequences'][0])
        {'aliases': ['chr1'],
        'assembly_unit': 'Primary Assembly',
        'genbank_ac': 'CM000663.1',
        'length': 249250621,
        'name': '1',
        'refseq_ac': 'NC_000001.10',
        'relationship': '=',
        'sequence_role': 'assembled-molecule'}
    """

    fn = pkg_resources.resource_filename(__name__, _assy_path_fmt.format(name=name))
    return json.load(gzip.open(fn, mode="rt", encoding="utf-8"))


def get_assemblies(names=[]):
    """Retrieves data from multiple assemblies.

    If assemblies are not specified, retrieves data from all available ones.

    Args:
        names (list of str, optional): The names of the assemblies to retrieve data for.

    Returns:
        dict: A dictionary of the form ``{assembly_name, : assembly_data}``, where the values
            are the dictionaries of assembly data as described in ``get_assembly()``.

    Examples:
        >>> assemblies = get_assemblies(names=['GRCh37.p13'])
        >>> assy = assemblies['GRCh37.p13']

        >>> assemblies = get_assemblies()
        >>> 'GRCh38.p2' in assemblies
        True
    """

    if names == []:
        names = get_assembly_names()
    return {a["name"]: a for a in (get_assembly(n) for n in names)}


def make_name_ac_map(assy_name, primary_only=False):
    """Creates a map from sequence names to accessions for a given assembly.

    Args:
        assy_name (str): The name of the assembly to make a map for.
        primary_only (bool, optional): Whether to include only primary sequences.
            Defaults to False.

    Returns:
        dict: A dictionary of the form ``{sequence_name : accession}`` for sequences in the given assembly,
            Where sequence_name and accession are both strings.

    Examples:
        >>> grch38p5_name_ac_map = make_name_ac_map('GRCh38.p5')
        >>> grch38p5_name_ac_map['1']
        'NC_000001.11'
    """

    return {
        s["name"]: s["refseq_ac"]
        for s in get_assembly(assy_name)["sequences"]
        if (not primary_only or _is_primary(s))
    }


def make_ac_name_map(assy_name, primary_only=False):
    """Creates a map from accessions to sequence names for a given assembly.

    Args:
        assy_name (str): The name of the assembly to make a map for.
        primary_only (bool, optional): Whether to include only primary sequences.
            Defaults to False.

    Returns:
        dict: A dictionary of the form ``{accesssion : sequence_name}`` for accessions in the given assembly,
            where accession and sequence_name are strings.


    Examples:
        >>> grch38p5_ac_name_map = make_ac_name_map('GRCh38.p5')
        >>> grch38p5_ac_name_map['NC_000001.11']
        '1'
    """

    return {
        s["refseq_ac"]: s["name"]
        for s in get_assembly(assy_name)["sequences"]
        if (not primary_only or _is_primary(s))
    }


############################################################################
# Internal functions


def _is_primary(s):
    """Indicates whether a sequence is a part of the primary assembly.

    Args:
        s (dict): A dictionary of sequence data, e.g. those in assembly['sequences'].

    Returns:
        bool: True if the sequence is part of the primary assembly, False otherwise.


    Examples:
        >>> _is_primary({'assembly_unit': 'Primary Assembly'})
        True

        >>> _is_primary({'assembly_unit': 'Something else entirely'})
        False
    """

    return s["assembly_unit"] == "Primary Assembly"
