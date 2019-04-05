import multiprocessing
import os

import pytest
import vcr

from bioutils.seqfetcher import fetch_seq, _fetch_seq_ensembl, _fetch_seq_ncbi, _add_eutils_api_key


@vcr.use_cassette
def test_fetch_seq():
    assert 1596 == len(fetch_seq('NP_056374.2'))

    assert 'MESRETLSSS' == fetch_seq('NP_056374.2',0,10)
    assert 'MESRETLSSS' == fetch_seq('NP_056374.2')[0:10]  # NOT RECOMMENDED

    assert 'ATCACACGTGCAGGAACCCTTTTCC' == fetch_seq('NC_000001.10', 2000000, 2000025)
    assert 'AAAATTAAATTAAAATAAATAAAAA' == fetch_seq('NG_032072.1', 0, 25)
    assert 'TTGTGTGTTAGGGTGCTCTAAGCAA' == fetch_seq('NW_003571030.1', 0, 25)
    assert 'GAATTCCTCGTTCACACAGTTTCTT' == fetch_seq('NT_113901.1', 0, 25)
    assert 'NNNNNNNNNNNNNNNNNNNNNNNNN' == fetch_seq('NC_000001.10', 0, 25)
    assert 'MESRETLSSSRQRGGESDFLPVSSA' == fetch_seq('NP_056374.2', 0, 25)
    assert 'GATCCACCTGCCTCAGCCTCCCAGA' == fetch_seq('GL000191.1', 0, 25)
    assert 'TTTATTTATTTTAGATACTTATCTC' == fetch_seq('KB663603.1', 0, 25)
    assert 'CGCCTCCCTTCCCCCTCCCCGCCCG' == fetch_seq('ENST00000288602', 0, 25)
    assert 'MAALSGGGGGGAEPGQALFNGDMEP' == fetch_seq('ENSP00000288602', 0, 25)


def test_add_eutils_api_key():
    try:
        url = 'http://test.com?boo=bar'
        assert _add_eutils_api_key(url) == url
        os.environ['NCBI_API_KEY'] = 'test-api-key'
        assert _add_eutils_api_key(url) == url + '&api_key=test-api-key'
    finally:
        try:
            os.environ.pop('NCBI_API_KEY')
        except KeyError:
            pass


def test_fetch_seq_errors():
    # Traceback (most recent call last):
    #    ...
    # RuntimeError: No sequence available for NM_9.9
    with pytest.raises(RuntimeError):
        fetch_seq('NM_9.9')

    # Traceback (most recent call last):
    #    ...
    # RuntimeError: No sequence fetcher for QQ01234
    with pytest.raises(RuntimeError):
        fetch_seq('QQ01234')


def _check1(_x):
    # small, fast query
    assert 'MESRETLSSS' == fetch_seq('NP_056374.2',0,10)

# no vcr!
@pytest.mark.network
def test_rate_limit():
    num_requests = num_threads = 5
    p = multiprocessing.Pool(num_threads)
    p.map(_check1, range(num_requests))
