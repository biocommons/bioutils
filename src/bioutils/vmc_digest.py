# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import hashlib

from .digest import Digest


ENC = "UTF-8"
DEFAULT_DIGEST_SIZE = 24


def vmc_digest(data, digest_size=DEFAULT_DIGEST_SIZE):
    """Returns the VMC Digest as a Digest object, which has both bytes and str (
    URL-safe, Base 64) representations.

    >>> d = vmc_digest("")

    # I can't figure out how to make this test work on Py 2 and 3 :-(
    >>> d                       # doctest: +SKIP
    b'\xcf\x83\xe15~\xef\xb8\xbd\xf1T(P\xd6m\x80'

    >>> str(d)
    'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc'

    >>> len(d), len(str(d))
    (24, 32)

    >>> str(vmc_digest("", 24))
    'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc'

    >>> vmc_digest("", 17) 
    Traceback (most recent call last):
    ...
    ValueError: digest_size must be a multiple of 3

    >>> vmc_digest("", 66) 
    Traceback (most recent call last):
    ...
    ValueError: digest_size must be between 0 and 63 (bytes)


    SHA-512 is 2x faster than SHA1 on modern 64-bit platforms.
    However, few appliations require 512 bits (64 bytes) of keyspace.
    That larger size translates into proportionally larger key size
    requirements, with attendant performance implications.  By
    truncating the SHA-512 digest [1], users may obtain a tunable
    level of collision avoidance.

    The string returned by this function is Base 64 encoded with
    URL-safe characters [2], making it suitable for use with URLs or
    filesystem paths.  Base 64 encoding results in an output string
    that is 4/3 the size of the input.  If the length of the input
    string is not divisible by 3, the output is right-padded with
    equal signs (=), which have no information content. Therefore,
    this function requires that digest_size is evenly divisible by 3.
    (The resulting vmc_digest will be 4/3*digest_size bytes.)

    According to [3], the probability of a collision using b bits with
    m messages (sequences) is:

       P(b, m) = m^2 / 2^(b+1).

    Note that the collision probability depends on the number of
    messages, but not their size.  Solving for the number of messages:

       m(b, P) = sqrt(P * 2^(b+1))

    Solving for the number of *bits*:

       b(m, P) = log2(m^2/P) - 1

    For various values of m and P, the number of *bytes* required
    according to b(m,P) rounded to next multiple of 3 is:

    +-------+----------+----------+----------+----------+----------+----------+
    |   #m  | P<=1e-24 | P<=1e-21 | P<=1e-18 | P<=1e-15 | P<=1e-12 | P<=1e-09 |
    +-------+----------+----------+----------+----------+----------+----------+
    | 1e+06 |    15    |    12    |    12    |    9     |    9     |    9     |
    | 1e+09 |    15    |    15    |    12    |    12    |    9     |    9     |
    | 1e+12 |    15    |    15    |    15    |    12    |    12    |    9     |
    | 1e+15 |    18    |    15    |    15    |    15    |    12    |    12    |
    | 1e+18 |    18    |    18    |    15    |    15    |    15    |    12    |
    | 1e+21 |    21    |    18    |    18    |    15    |    15    |    15    |
    | 1e+24 |    21    |    21    |    18    |    18    |    15    |    15    |
    +-------+----------+----------+----------+----------+----------+----------+

    For example, given 1e+18 expected messages and a desired collision
    probability < 1e-15, we use digest_size = 15 (bytes).

    [1] http://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.180-4.pdf
    [2] https://tools.ietf.org/html/rfc3548#section-4
    [3] http://stackoverflow.com/a/4014407/342839
    [4] http://stackoverflow.com/a/22029380/342839
    [5] http://preshing.com/20110504/hash-collision-probabilities/
    [6] https://en.wikipedia.org/wiki/Birthday_problem

    """ 

    # TODO: Consider relaxing %3 constraint and stripping padding
    if digest_size % 3 != 0:
        raise ValueError("digest_size must be a multiple of 3")
    if not 0 <= digest_size <= 63:
        raise ValueError("digest_size must be between 0 and 63 (bytes)")

    sha512 = Digest(hashlib.sha512(data.encode(ENC)).digest())
    return sha512[:digest_size]
