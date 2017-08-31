import pytest
import vcr

from bioutils.seqfetcher import fetch_seq, _fetch_seq_ensembl, _fetch_seq_ncbi


@vcr.use_cassette
def test_fetch_seq():
    assert 1596 == len(fetch_seq('NP_056374.2'))

    assert 'MESRETLSSS' == fetch_seq('NP_056374.2',0,10)
    assert 'MESRETLSSS' == fetch_seq('NP_056374.2')[0:10]  # NOT RECOMMENDED

    email = "tests@biocommons.org"
    tool = "bioutils"

    assert 'ATCACACGTGCAGGAACCCTTTTCC' == fetch_seq('NC_000001.10', 2000000, 2000025, email=email, tool=tool)
    assert 'AAAATTAAATTAAAATAAATAAAAA' == fetch_seq('NG_032072.1', 0, 25, email=email, tool=tool)
    assert 'TTGTGTGTTAGGGTGCTCTAAGCAA' == fetch_seq('NW_003571030.1', 0, 25, email=email, tool=tool)
    assert 'GAATTCCTCGTTCACACAGTTTCTT' == fetch_seq('NT_113901.1', 0, 25, email=email, tool=tool)
    assert 'NNNNNNNNNNNNNNNNNNNNNNNNN' == fetch_seq('NC_000001.10', 0, 25, email=email, tool=tool)
    assert 'MESRETLSSSRQRGGESDFLPVSSA' == fetch_seq('NP_056374.2', 0, 25, email=email, tool=tool)
    assert 'GATCCACCTGCCTCAGCCTCCCAGA' == fetch_seq('GL000191.1', 0, 25, email=email, tool=tool)
    assert 'TTTATTTATTTTAGATACTTATCTC' == fetch_seq('KB663603.1', 0, 25, email=email, tool=tool)
    assert 'CGCCTCCCTTCCCCCTCCCCGCCCG' == fetch_seq('ENST00000288602', 0, 25, email=email, tool=tool)
    assert 'MAALSGGGGGGAEPGQALFNGDMEP' == fetch_seq('ENSP00000288602', 0, 25, email=email, tool=tool)
    

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
