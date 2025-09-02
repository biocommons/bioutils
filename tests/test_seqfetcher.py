import multiprocessing
import os

import pytest
import vcr

from bioutils.seqfetcher import _add_eutils_api_key, _fetch_seq_ncbi, fetch_seq, enst_default_seq_type


@pytest.fixture(autouse=True)
def clear_env():
    """Some tests in this module assume that the default utils access configs are
    active. If you execute tests in an environment with an existing `NCBI_API_KEY` env
    var, those tests will fail unless we first remove that variable.
    """
    if "NCBI_API_KEY" in os.environ:
        del os.environ["NCBI_API_KEY"]


@vcr.use_cassette
def test_fetch_seq():
    assert 1596 == len(fetch_seq("NP_056374.2"))

    assert "MESRETLSSS" == fetch_seq("NP_056374.2", 0, 10)
    assert "MESRETLSSS" == fetch_seq("NP_056374.2")[0:10]  # NOT RECOMMENDED

    assert "ATCACACGTGCAGGAACCCTTTTCC" == fetch_seq("NC_000001.10", 2000000, 2000025)
    assert "AAAATTAAATTAAAATAAATAAAAA" == fetch_seq("NG_032072.1", 0, 25)
    assert "TTGTGTGTTAGGGTGCTCTAAGCAA" == fetch_seq("NW_003571030.1", 0, 25)
    assert "GAATTCCTCGTTCACACAGTTTCTT" == fetch_seq("NT_113901.1", 0, 25)
    assert "NNNNNNNNNNNNNNNNNNNNNNNNN" == fetch_seq("NC_000001.10", 0, 25)
    assert "MESRETLSSSRQRGGESDFLPVSSA" == fetch_seq("NP_056374.2", 0, 25)
    assert "GATCCACCTGCCTCAGCCTCCCAGA" == fetch_seq("GL000191.1", 0, 25)
    assert "TTTATTTATTTTAGATACTTATCTC" == fetch_seq("KB663603.1", 0, 25)
    assert "CCGCTCGGGCCCCGGCTCTCGGTTA" == fetch_seq("ENST00000288602.11", 0, 25)
    assert "MAALSGGGGGGAEPGQALFNGDMEP" == fetch_seq("ENSP00000288602", 0, 25)


ENST00000617537_470_480 = {
    # In [16]: s_gen[470:480], s_cdna[470:480], s_cds[470:480]
    # Out[16]: ("TAGGTATGCA", "TAGGGTGTGT", "TGACATTTGT")
    "genomic": "TAGGTATGCA",
    "cdna": "TAGGGTGTGT",
    "cds": "TGACATTTGT",
}


@vcr.use_cassette
def test_fetch_ENST00000617537_noenv(caplog, monkeypatch):
    """ensure expected lengths for ENST00000617537 with ENST_DEFAULT_SEQ_TYPE unset"""
    monkeypatch.delenv("ENST_DEFAULT_SEQ_TYPE", raising=False)
    ac = "ENST00000617537"
    assert ENST00000617537_470_480[enst_default_seq_type] == fetch_seq(ac, start_i=470, end_i=480)
    assert "Transcript type not specified or set in ENST_DEFAULT_SEQ_TYPE" in caplog.text
    assert ENST00000617537_470_480["genomic"] == fetch_seq(ac, start_i=470, end_i=480, seq_type="genomic")
    assert ENST00000617537_470_480["cdna"] == fetch_seq(ac, start_i=470, end_i=480, seq_type="cdna")
    assert ENST00000617537_470_480["cds"] == fetch_seq(ac, start_i=470, end_i=480, seq_type="cds")


@vcr.use_cassette
def test_fetch_ENST00000617537_env(caplog, monkeypatch):
    """ensure expected lengths for ENST00000617537 with ENST_DEFAULT_SEQ_TYPE set"""
    user_enst_default_type = "cds"  # intentionally != enst_default_seq_type to ensure use
    monkeypatch.setenv("ENST_DEFAULT_SEQ_TYPE", user_enst_default_type)
    ac = "ENST00000617537"
    assert ENST00000617537_470_480[user_enst_default_type] == fetch_seq(ac, start_i=470, end_i=480)
    assert "Transcript type not specified or set in ENST_DEFAULT_SEQ_TYPE" not in caplog.text
    assert ENST00000617537_470_480["genomic"] == fetch_seq(ac, start_i=470, end_i=480, seq_type="genomic")
    assert ENST00000617537_470_480["cdna"] == fetch_seq(ac, start_i=470, end_i=480, seq_type="cdna")
    assert ENST00000617537_470_480["cds"] == fetch_seq(ac, start_i=470, end_i=480, seq_type="cds")


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
        try:
            os.environ.pop("NCBI_API_KEY")
        except KeyError:
            pass


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
    assert "MESRETLSSS" == fetch_seq("NP_056374.2", 0, 10)


# no vcr!
@pytest.mark.network
def test_rate_limit():
    num_requests = num_threads = 5
    p = multiprocessing.Pool(num_threads)
    p.map(_check1, range(num_requests))
    p.close()
    p.join()
