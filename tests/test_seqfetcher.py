import contextlib
import multiprocessing
import os

import pytest
import vcr

from bioutils.seqfetcher import (
    _add_eutils_api_key,
    _fetch_seq_ncbi,
    fetch_seq,
)


@pytest.fixture(autouse=True)
def clear_env():
    """Clear environment.

    Some tests in this module assume that the default utils access configs are
    active. If you execute tests in an environment with an existing `NCBI_API_KEY` env
    var, those tests will fail unless we first remove that variable.
    """
    if "NCBI_API_KEY" in os.environ:
        del os.environ["NCBI_API_KEY"]


@vcr.use_cassette
def test_fetch_seq():
    assert len(fetch_seq("NP_056374.2")) == 1596

    assert fetch_seq("NP_056374.2", 0, 10) == "MESRETLSSS"
    assert fetch_seq("NP_056374.2")[0:10] == "MESRETLSSS"  # NOT RECOMMENDED

    assert fetch_seq("NC_000001.10", 2000000, 2000025) == "ATCACACGTGCAGGAACCCTTTTCC"
    assert fetch_seq("NG_032072.1", 0, 25) == "AAAATTAAATTAAAATAAATAAAAA"
    assert fetch_seq("NW_003571030.1", 0, 25) == "TTGTGTGTTAGGGTGCTCTAAGCAA"
    assert fetch_seq("NT_113901.1", 0, 25) == "GAATTCCTCGTTCACACAGTTTCTT"
    assert fetch_seq("NC_000001.10", 0, 25) == "NNNNNNNNNNNNNNNNNNNNNNNNN"
    assert fetch_seq("NP_056374.2", 0, 25) == "MESRETLSSSRQRGGESDFLPVSSA"
    assert fetch_seq("GL000191.1", 0, 25) == "GATCCACCTGCCTCAGCCTCCCAGA"
    assert fetch_seq("KB663603.1", 0, 25) == "TTTATTTATTTTAGATACTTATCTC"
    assert fetch_seq("ENST00000288602.11", 0, 25) == "CCGCTCGGGCCCCGGCTCTCGGTTA"
    assert fetch_seq("ENSP00000288602", 0, 25) == "MAALSGGGGGGAEPGQALFNGDMEP"


@vcr.use_cassette
def test_fetch_seq_ncbi_invalid_positions():
    with pytest.raises(RuntimeError) as excinfo:
        _fetch_seq_ncbi("NP_001230161.1", 3190, 3190)
    assert "invalid sequence or start or end position" in str(excinfo.value)


@vcr.use_cassette
def test_add_eutils_api_key():
    try:
        url = "http://test.com?boo=bar"
        assert _add_eutils_api_key(url) == url
        os.environ["NCBI_API_KEY"] = "test-api-key"
        assert _add_eutils_api_key(url) == url + "&api_key=test-api-key"
    finally:
        with contextlib.suppress(KeyError):
            os.environ.pop("NCBI_API_KEY")


@vcr.use_cassette
@pytest.mark.network
def test_fetch_seq_errors():
    # Traceback (most recent call last):
    #    ...
    # RuntimeError: No sequence available for NM_9.9
    with pytest.raises(RuntimeError):
        fetch_seq("NM_9.9")

    # Traceback (most recent call last):
    #    ...
    # RuntimeError: No sequence fetcher for QQ01234
    with pytest.raises(RuntimeError):
        fetch_seq("QQ01234")


def _check1(_x):
    # small, fast query
    assert fetch_seq("NP_056374.2", 0, 10) == "MESRETLSSS"


# no vcr!
@pytest.mark.network
def test_rate_limit():
    num_requests = num_threads = 5
    p = multiprocessing.Pool(num_threads)
    p.map(_check1, range(num_requests))
