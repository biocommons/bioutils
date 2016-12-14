import base64
import hashlib

import six


def truncated_sha512(binary_data, digest_size):
    """returns SHA-512 *binary* digest, truncated to digest_size
    *bytes*"""
    return hashlib.sha512(binary_data).digest()[:digest_size]


def truncated_digest(binary_data, digest_size):
    """Returns the URL-safe, Base 64 encoded string representation of a
    SHA-512 digest truncated to digest_size *bytes*

    >>> truncated_digest(b"", 18)
    'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc'

    >>> truncated_digest(b"", 17) 
    Traceback (most recent call last):
    ...
    ValueError: digest_size must be a multiple of 3

    >>> truncated_digest(b"", 66) 
    Traceback (most recent call last):
    ...
    ValueError: digest_size must be between 1 and 63 (bytes)


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
    (The resulting truncated_digest will be 4/3*digest_size bytes.)

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

    For 1e+18 messages and a desired collision probability < 1e-15, we
    use digest_size = 15.

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
    if not 1 < digest_size <= 63:
        raise ValueError("digest_size must be between 1 and 63 (bytes)")

    sha512t = truncated_sha512(binary_data=binary_data, digest_size=24)
    return base64.urlsafe_b64encode(sha512t).decode("ASCII")


if __name__ == "__main__":  # pragma: nocover
    import math
    import prettytable

    def B(P, m):
        """return number of *bytes* needed to achieve a collision probability
        P for m messages"""
        return math.ceil((math.log2(m / P) - 1) / 8 / 3) * 3

    m_bins = [1E6, 1E9, 1E12, 1E15, 1E18, 1E21, 1E24]
    P_bins = [1E-24, 1E-21, 1E-18, 1E-15, 1E-12, 1E-9]
    field_names = ["#m"] + ["P<={P}".format(P=P) for P in P_bins]
    pt = prettytable.PrettyTable(field_names=field_names)
    for n_m in m_bins:
        pt.add_row(["{:g}".format(n_m)] + [B(P, n_m) for P in P_bins])
    print(pt)
