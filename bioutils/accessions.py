# -*- coding: utf-8 -*-
"""simple routines to deal with accessions, identifiers, etc.

biocommons terminology: an identifier is composed of a *namespace* and
an *accession*. The namespace is a string, composed of any character
other than colon (:). The accession is a string without character set
restriction.  An accession is expected to be unique within the
namespace; there is no expectation of uniqueness of accessions across
namespaces.

Identifier := <Namespace, Accession>
Namespace := [^:]+
Accession := \w+

Some sample serializations of Identifiers:

json: {"namespace": "RefSeq", "accession": "NM_000551.3"}
xml: <Identifier namespace="RefSeq" accession="NM_000551.3"/>
string: "RefSeq:NM_000551.3"

The string form may be used as a CURIE, in which case the document in
which the CURIE is used must contain a map of {namespace => uri}.

"""

from __future__ import absolute_import, division, print_function, unicode_literals

import re

from six import iteritems

from .exceptions import BioutilsError


_ensembl_species_prefixes = "|".join("""ENS ENSACA ENSAME ENSAMX
ENSANA ENSAPL ENSBTA ENSCAF ENSCAN ENSCAP ENSCAT ENSCCA ENSCEL ENSCGR
ENSCGR ENSCHI ENSCHO ENSCIN ENSCJA ENSCLA ENSCPO ENSCSA ENSCSAV ENSDAR
ENSDNO ENSDOR ENSEBU ENSECA ENSEEU ENSETE ENSFAL ENSFCA ENSFDA ENSGAC
ENSGAL ENSGGO ENSGMO ENSHGLF ENSHGLM ENSJJA ENSLAC ENSLAF ENSLOC
ENSMAU ENSMEU ENSMFA ENSMGA ENSMIC ENSMLE ENSMLU ENSMMU ENSMNE ENSMOC
ENSMOD ENSMPU ENSMUS ENSNGA ENSNLE ENSOAN ENSOAR ENSOCU ENSODE ENSOGA
ENSONI ENSOPR ENSORL ENSPAN ENSPCA ENSPCO ENSPEM ENSPFO ENSPMA ENSPPA
ENSPPR ENSPPY ENSPSI ENSPTI ENSPTR ENSPVA ENSRBI ENSRNO ENSRRO ENSSAR
ENSSBO ENSSCE ENSSHA ENSSSC ENSSTO ENSTBE ENSTGU ENSTNI ENSTRU ENSTSY
ENSTTR ENSVPA ENSXET ENSXMA FB MGP_129S1SvImJ_ MGP_AJ_ MGP_AKRJ_
MGP_BALBcJ_ MGP_C3HHeJ_ MGP_C57BL6NJ_ MGP_CAROLIEiJ_ MGP_CASTEiJ_
MGP_CBAJ_ MGP_DBA2J_ MGP_FVBNJ_ MGP_LPJ_ MGP_NODShiLtJ_ MGP_NZOHlLtJ_
MGP_PWKPhJ_ MGP_PahariEiJ_ MGP_SPRETEiJ_ MGP_WSBEiJ_""".split())
_ensembl_feature_types_re = r"E|FM|G|GT|P|R|T"
_ensembl_re = r"^(?:{})(?:{}){}$".format(
    _ensembl_species_prefixes, _ensembl_feature_types_re, r"\d{11}(?:\.\d+)?")

# map of regexp => namespace
# TODO: make this namespace => [regexps] for clarity
# namespaces follow convention of identifiers.org
ac_namespace_regexps = {
    # https://uswest.ensembl.org/info/genome/stable_ids/prefixes.html
    # [species prefix][feature type prefix][a unique eleven digit number]
    # N.B. The regexp at http://identifiers.org/ensembl appears broken:
    # 1) Human only; 2) escaped backslashes (\\d rather than \d).
    _ensembl_re: "ensembl",

    # http://identifiers.org/insdc/
    # P12345, a UniProtKB accession matches the miriam regexp but shouldn't (I think)
    r"^([A-Z]\d{5}|[A-Z]{2}\d{6}|[A-Z]{4}\d{8}|[A-J][A-Z]{2}\d{5})(\.\d+)?$":
    "insdc",

    # http://identifiers.org/refseq/
    # https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/
    r"^((AC|AP|NC|NG|NM|NP|NR|NT|NW|XM|XP|XR|YP|ZP)_\d+|(NZ\_[A-Z]{4}\d+))(\.\d+)?$":
    "refseq",

    # Uniprot
    # http://identifiers.org/uniprot/
    # https://www.uniprot.org/help/accession_numbers
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$":
    "uniprot",
}

ac_namespace_regexps = {re.compile(k): v for k, v in iteritems(ac_namespace_regexps)}


def chr22XY(c):
    """force to name from 1..22, 23, 24, X, Y, M 
    to in chr1..chr22, chrX, chrY, chrM
    str or ints accepted

    >>> chr22XY('1')
    'chr1'
    >>> chr22XY(1)
    'chr1'
    >>> chr22XY('chr1')
    'chr1'
    >>> chr22XY(23)
    'chrX'
    >>> chr22XY(24)
    'chrY'
    >>> chr22XY("X")
    'chrX'
    >>> chr22XY("23")
    'chrX'
    >>> chr22XY("M")
    'chrM'

    """
    c = str(c)
    if c[0:3] == 'chr':
        c = c[3:]
    if c == '23':
        c = 'X'
    if c == '24':
        c = 'Y'
    return 'chr' + c


def infer_namespace(ac):
    """Infer the single namespace of the given accession

    This function is convenience wrapper around infer_namespaces().
    Returns:
      * None if no namespaces are inferred
      * The (single) namespace if only one namespace is inferred
      * Raises an exception if more than one namespace is inferred

    >>> infer_namespace("ENST00000530893.6")
    'ensembl'

    >>> infer_namespace("NM_01234.5")
    'refseq'

    >>> infer_namespace("A2BC19")
    'uniprot'

    N.B. The following test is disabled because Python 2 and Python 3
    handle doctest exceptions differently. :-(
    X>>> infer_namespace("P12345")
    Traceback (most recent call last):
    ...
    bioutils.exceptions.BioutilsError: Multiple namespaces possible for P12345

    >>> infer_namespace("BOGUS99") is None
    True

    """

    namespaces = infer_namespaces(ac)
    if not namespaces:
        return None
    if len(namespaces) > 1:
        raise BioutilsError("Multiple namespaces possible for {}".format(ac))
    return namespaces[0]


def infer_namespaces(ac):
    """infer possible namespaces of given accession based on syntax
    Always returns a list, possibly empty

    >>> infer_namespaces("ENST00000530893.6")
    ['ensembl']
    >>> infer_namespaces("ENST00000530893")
    ['ensembl']
    >>> infer_namespaces("ENSQ00000530893")
    []
    >>> infer_namespaces("NM_01234")
    ['refseq']
    >>> infer_namespaces("NM_01234.5")
    ['refseq']
    >>> infer_namespaces("NQ_01234.5")
    []
    >>> infer_namespaces("A2BC19")
    ['uniprot']
    >>> sorted(infer_namespaces("P12345"))
    ['insdc', 'uniprot']
    >>> infer_namespaces("A0A022YWF9")
    ['uniprot']


    """
    return [v for k, v in iteritems(ac_namespace_regexps) if k.match(ac)]


def prepend_chr(chr):
    """prefix chr with 'chr' if not present

    Users are strongly encouraged to NOT use this function. Added a
    'chr' prefix means that you're using a name that is not consistent
    with authoritative assembly records.

    >>> prepend_chr('22')
    'chr22'

    >>> prepend_chr('chr22')
    'chr22'

    """
    return chr if chr[0:3] == 'chr' else 'chr' + chr


def strip_chr(chr):
    """remove 'chr' prefix if it exists

    >>> strip_chr('22')
    '22'

    >>> strip_chr('chr22')
    '22'

    """
    return chr[3:] if chr[0:3] == 'chr' else chr


## <LICENSE>
## Copyright 2014 Bioutils Contributors (https://bitbucket.org/biocommons/bioutils)
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>
