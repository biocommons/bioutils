# -*- coding: utf-8 -*-
"""interface to identifiers.org data


N.B. This module is on a probation of sorts.  Although it works, the
regular expressions in identifiers.org have significant problems that
impede their utility, and therefore the utility of this module.  As of
Oct 2018, seven regular expressions aren't even valid.  Many others
are way too permissive, which leads to many non-specific hits for
extremely niche databases. Even removing these (see namespace_sets
below), there are important collections with overly permissive
regexps. Examples:

  In [3]: infer_namespaces("1")
  Out[3]: ['clinvar', 'ensembl', 'hgnc']
  
  In [4]: infer_namespaces("1111")
  Out[4]: ['clinvar', 'ensembl', 'hgnc', 'medgen']
  
  In [5]: infer_namespaces("A")
  Out[5]: ['ensembl']



The following loosey-goosey grammar is imagined:

  FQIdentifier := Namespace ':' Accession
  UQIdentifier := Accession
  Identifier := FQIdentifier | UQIdentifier
  Namespace := [^:]+
  Accession := \w+

where

  * Namespace is a string that refers to a database authority
  * Accession is a unique key within a Namespace 
  * FQIdentifier (Fully Qualified Identifier) is a string composed of a Namespace and Accession that is globally unique
  * UQIdentifier (Unqualified Identifier) is a string that is NOT guaranteed globally unique
  * Identifier is a FQIdentifier or UQIdentifier

Ideally, the community will migrated to using FQIdentifiers
exclusively. Not coincidentally, these look like CURIEs.

"""

from __future__ import absolute_import, division, print_function, unicode_literals

import datetime
import gzip
import json
import logging
import re

import pkg_resources

_logger = logging.getLogger()

_miriam_registry_path = '_data/identifiers/registry.json.gz'
_identfier_data = None

_namespace_sets = {
    # the common set is intended to contain namespaces that are
    # frequently encountered and for which the regexp is relatively
    # specific. E.g., pubmed is common, but the regexp is extremely
    # non-specific.
    "common": set(['clinvar', 'dbsnp', 'ensembl', 'go',
                   'hgnc', 'insdc', 'lrg', 'medgen',
                   'pharmgkb.disease', 'pharmgkb.drug',
                   'pharmgkb.gene', 'pharmgkb.pathways',
                   'refseq', 'so'])
}
default_namespace_set = "common"


class _IdentifierData:

    """**Internal Use Only**

    This class is a container for runtime data relevant for dealing
    with identifiers and accessions.  The module creates an instance
    of this object on first use.

    """

    def __init__(self, include_obsolete=False, namespace_filter=None):
        self._load_miriam(include_obsolete=include_obsolete, namespace_filter=namespace_filter)
        self.re_to_namespace = {re.compile(c["pattern"]):
                                c["namespace"] for c in self.miriam["collections"].values()}


    def _load_miriam(self, include_obsolete, namespace_filter):
        stream = gzip.GzipFile(fileobj=pkg_resources.resource_stream(__name__, _miriam_registry_path))
        self.miriam = json.load(stream)

        collections = self.miriam["collections"]
        _logger.info("Loaded miriam data (modified {}); {} namespaces".format(
            datetime.datetime.fromtimestamp(self.miriam["last_modified"]).isoformat(),
            len(collections)))
        if not include_obsolete:
            delkeys = [k for k, c in collections.items() if c["obsolete"]]
            for k in delkeys:
                del collections[k]
            _logger.debug("removed {} obsolete namespaces".format(len(delkeys)))
        if namespace_filter:
            delkeys = [k for k, c in collections.items() if c["namespace"] not in namespace_filter]
            for k in delkeys:
                del collections[k]
            _logger.debug("removed {} namespaces not in namespace filter".format(len(delkeys)))
        _logger.debug("{} namespaces remain".format(len(collections.keys())))
        self.miriam["collections"] = collections


    def infer_namespaces(self, acc):
        """infer possible source namespaces based on syntax of accession, acc

        """
        return [v for k, v in self.re_to_namespace.items() if k.match(acc)]

    def infer_namespace(self, acc):
        """return single namespace of given accession, None if none could be
        found, or raise exception if multiple found

        """
        namespaces = self.infer_namespaces(acc)
        if len(namespaces) > 1:
            raise ValueError("{} matches multiple namespace patterns".format(acc))
        if namespaces:
            return namespaces[0]
        return None


def _get_identifierdata():
    global _identfier_data
    if not _identfier_data:
        _identfier_data = _IdentifierData(namespace_filter=_namespace_sets[default_namespace_set])
    return _identfier_data


def infer_namespaces(acc):
    """infer possible source namespaces based on syntax of accession,
    acc

    >>> infer_namespaces("NM_01234.5")
    ['refseq']

    """

    return _get_identifierdata().infer_namespaces(acc)


def infer_namespace(acc):
    """return single namespace of given accession, None if none could be
    found, or raise exception if multiple found

    >>> infer_namespace("NM_01234.5")
    'refseq'

    >>> infer_namespace("!")

    >>> infer_namespace("01234")
    Traceback (most recent call last):
    ...
    ValueError: 01234 matches multiple namespace patterns

    """
    return _get_identifierdata().infer_namespace(acc)

