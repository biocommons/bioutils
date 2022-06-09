# -*- coding: utf-8 -*-
"""Provides sequence fetching from NCBI and Ensembl.

"""

from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import os
import random
import re
import time

import requests

_logger = logging.getLogger(__name__)

# Reece requested registration on 2017-09-03
ncbi_tool = "bioutils"
ncbi_email = "biocommons-dev@googlegroups.com"
retry_limit = 3


def fetch_seq(ac, start_i=None, end_i=None):
    """Fetches sequences and subsequences from NCBI eutils and Ensembl REST interfaces.

    Args:
        ac (str): The accession of the sequence to fetch.
        start_i (int, optional): The start index (interbase coordinates) of the subsequence to fetch. Defaults to ``None``.
            It is recommended to retrieve a subsequence by providing an index here, rather than by
            Python slicing the whole sequence.
        end_i (int, optional): The end index (interbase coordinates) of the subsequence to fetch. Defaults to ``None``.
            It is recommended to retrieve a subsequence by providing an index here, rather than by
            Python slicing the whole sequence.

    Returns:
        str: The requested sequence.

    Raises:
        RuntimeError: If the syntax doesn't match that of any of the databases.
        RuntimeError: If the request to the database fails.

    Examples:
        >>> len(fetch_seq('NP_056374.2'))
        1596

        >>> fetch_seq('NP_056374.2',0,10)   # This!
        'MESRETLSSS'

        >>> fetch_seq('NP_056374.2')[0:10]  # Not this!
        'MESRETLSSS'

        # Providing intervals is especially important for large sequences:

        >>> fetch_seq('NC_000001.10',2000000,2000030)
        'ATCACACGTGCAGGAACCCTTTTCCAAAGG'

        # This call will pull back 30 bases plus overhead; without the
        # interval, one would receive 250MB of chr1 plus overhead!

        # Essentially any RefSeq, Genbank, BIC, or Ensembl sequence may be
        # fetched.

        >>> fetch_seq('NM_9.9')
        Traceback (most recent call last):
        ...
        RuntimeError: No sequence available for NM_9.9

        >>> fetch_seq('QQ01234')
        Traceback (most recent call last):
        ...
        RuntimeError: No sequence fetcher for QQ01234
    """

    ac_dispatch = [
        {
            "re": re.compile(r"^(?:AC|N[CGMPRTW])_|^[A-L]\w\d|^U\d"),
            "fetcher": _fetch_seq_ncbi,
        },
        {"re": re.compile(r"^ENS[TP]\d+"), "fetcher": _fetch_seq_ensembl},
    ]

    eligible_fetchers = [dr["fetcher"] for dr in ac_dispatch if dr["re"].match(ac)]

    if len(eligible_fetchers) == 0:
        raise RuntimeError("No sequence fetcher for {ac}".format(ac=ac))

    if len(eligible_fetchers) >= 1:  # pragma: nocover (no way to test)
        _logger.debug(
            "Multiple sequence fetchers found for " "{ac}; using first".format(ac=ac)
        )

    fetcher = eligible_fetchers[0]
    _logger.debug("fetching {ac} with {f}".format(ac=ac, f=fetcher))

    try:
        return fetcher(ac, start_i, end_i)
    except requests.RequestException as ex:
        raise RuntimeError("Failed to fetch {ac} ({ex})".format(ac=ac, ex=ex))


# ###########################################################################
# Internal functions


def _fetch_seq_ensembl(ac, start_i=None, end_i=None):
    """Fetch sequence slice from Ensembl public REST interface.

    Args:
        ac (str): The accession of the sequence to fetch.
        start_i (int, optional): The start index (interbase coordinates) of the subsequence to fetch.
            Defaults to None.
        end_i (int, optional): The end index (interbase coordinates) of the subsequence to fetch.
            Defaults to None.

    Returns:
        str: The requested (sub)sequence

    Raises:
        RequestException: if request is unsuccessful.
        KeyError: if Ensembl API returns a different version than requested

    Note:
        The Ensembl REST interface does not currently accept intervals, so this method
        slices the sequence locally.

    Examples:
        >> len(_fetch_seq_ensembl('ENSP00000288602'))
        766

        >> _fetch_seq_ensembl('ENSP00000288602',0,10)
        u'MAALSGGGGG'

        >> _fetch_seq_ensembl('ENSP00000288602')[0:10]
        u'MAALSGGGGG'

        >> ac = 'ENSP00000288602'
        >> _fetch_seq_ensembl(ac ,0, 10) == _fetch_seq_ensembl(ac)[0:10]
        True
    """

    # Ensembl API only takes transcript IDs (without version) and returns the latest one
    # So we need to strip the transcript version, then check if what was returned was the one we requested
    version = None
    m = re.match(r"^(ENST\d+)\.(\d+)$", ac)
    if m:
        ac, version = m.groups()
        version = int(version)

    url = f"http://rest.ensembl.org/sequence/id/{ac}"
    r = requests.get(url, headers={"Content-Type": "application/json"})
    r.raise_for_status()
    data = r.json()
    if version is not None:
        latest_version = data["version"]
        if version != latest_version:
            msg = f"Ensembl API only provides {ac} version ({latest_version}), requested: {version}"
            raise KeyError(msg)

    seq = data["seq"]
    return seq if (start_i is None or end_i is None) else seq[start_i:end_i]


def _fetch_seq_ncbi(ac, start_i=None, end_i=None):
    """Fetches sequences from NCBI using the eutils interface.

    Args:
        ac (str): The accession of the sequence to fetch.
        start_i (int, optional): The start index (interbase coordinates) of the subsequence to fetch.
            Defaults to None.
        end_i (int, optional): The end index (interbase coordinates) of the subsequence to fetch.
            Defaults to None.

    Returns:
        str: The requested (sub)sequence

    Raises:
        RequestException: if request is unsuccessful.

    Notes:
        An interbase interval may be optionally provided with start_i and
        end_i. NCBI eutils will return just the requested subsequence,
        which might greatly reduce payload sizes (especially with
        chromosome-scale sequences).

        The request includes `tool` and `email` arguments to identify the
        caller as the bioutils package.  According to
        https://www.ncbi.nlm.nih.gov/books/NBK25497/, these values should
        correspond to the library, not the library client.  Using the
        defaults is recommended.  Nonetheless, callers may set
        ``bioutils.seqfetcher.ncbi_tool`` and
        ``bioutils.seqfetcher.ncbi_email`` to custom values if that is
        desired.

    Examples:
        >> len(_fetch_seq_ncbi('NP_056374.2'))
        1596

        Pass the desired interval rather than using Python's [] slice
        operator.

        >> _fetch_seq_ncbi('NP_056374.2',0,10)
        'MESRETLSSS'

        >> _fetch_seq_ncbi('NP_056374.2')[0:10]
        'MESRETLSSS'

        >> _fetch_seq_ncbi('NP_056374.2',0,10) == _fetch_seq_ncbi('NP_056374.2')[0:10]
        True
    """

    db = "protein" if ac[1] == "P" else "nucleotide"
    url_fmt = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        "db={db}&id={ac}&rettype=fasta"
    )

    if start_i is None or end_i is None:
        url = url_fmt.format(db=db, ac=ac)
    else:
        url_fmt += "&seq_start={start}&seq_stop={stop}"
        url = url_fmt.format(db=db, ac=ac, start=start_i + 1, stop=end_i)

    url += "&tool={tool}&email={email}".format(tool=ncbi_tool, email=ncbi_email)

    url = _add_eutils_api_key(url)

    n_retries = 0
    while True:
        resp = requests.get(url)
        if resp.ok:
            seq = "".join(resp.text.splitlines()[1:])
            return seq
        elif resp.status_code == 400:
            # Invalid sequence or start/stop position for that sequence
            raise RuntimeError(
                "Fetching sequence {} with start index {} and end index {} failed, invalid sequence "
                "or start or end position".format(ac, start_i, end_i)
            )
        if n_retries >= retry_limit:
            break
        if n_retries == 0:
            _logger.warning("Failed to fetch {}".format(url))
        sleeptime = random.randint(n_retries, 3) ** n_retries
        _logger.warning(
            "Failure {}/{}; retry in {} seconds".format(
                n_retries, retry_limit, sleeptime
            )
        )
        time.sleep(sleeptime)
        n_retries += 1
    # Falls through only on failure
    resp.raise_for_status()


def _add_eutils_api_key(url):
    """Adds an eutils api key to the query.

    Args:
        url (str): The query url without the api key.

    Returns:
        str: The query url with the api key, if one is stored in the environment variable
            ``NCBI_API_KEY``, otherwise it is unaltered.
    """

    apikey = os.environ.get("NCBI_API_KEY")
    if apikey:
        url += "&api_key={apikey}".format(apikey=apikey)
    return url


# So that I don't forget why I didn't use ebi too:
# $ curl 'http://www.ebi.ac.uk/ena/data/view/AM407889.1&display=fasta'
# >ENA|AM407889|AM407889.2 Medicago sativa partial mRNA ...
# AACGTATCACACTTCTTCTCCATTTCTTTTTCTTACATCTTCTCTCTACAAATTCATTTC
# Note that we requested .1, got .2.  Implicit behavior bites again.

if __name__ == "__main__":  # pragma: nocover
    import doctest

    doctest.testmod()
