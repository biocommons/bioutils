import base64

import six


@six.python_2_unicode_compatible
class Digest(bytes):
    """Represents a sliceable binary digest, with support for encoding and
    decoding using printable characters.

    Stringified representations use URL-safe base64 encodings. See
    https://tools.ietf.org/html/rfc4648#page-7

    >>> import hashlib

    >>> b = hashlib.sha512().digest()
    >>> len(b)
    64

    >>> d = Digest(b)
    >>> len(d)
    64
    >>> len(d.as_base64us())
    88
    >>> d.as_base64us()
    'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXcg_SpIdNs6c5H0NE8XYXysP-DGNKHfuwvY7kxvUdBeoGlODJ6-SfaPg=='
    >>> str(d)
    'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXcg_SpIdNs6c5H0NE8XYXysP-DGNKHfuwvY7kxvUdBeoGlODJ6-SfaPg=='

    >>> d24 = d[:24]
    >>> len(d24)
    24
    >>> len(d24.as_base64us())
    32
    >>> d24.as_base64us()
    'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc'

    >>> e = Digest.from_base64us(d.as_base64us())
    >>> e == d
    True

    >>> e24 = Digest.from_base64us(d24.as_base64us())
    >>> e24 == d24
    True

    """

    def __str__(self):
        """returns digest as URL-safe, base64 encoded digest, as Unicode

        """
        return self.as_base64us().decode("ASCII")


    # TODO: Consider requiring slice start == None or 0, and len % 3 == 0
    # Slicing %3 != 0 => strings will having suffix differences
    if six.PY2:                 # pragma: nocover
        def __getslice__(self, start, end):
            return Digest(bytes.__getslice__(self, start, end))
    else:                       # pragma: nocover
        def __getitem__(self, key):
            return Digest(bytes.__getitem__(self, key))


    def as_base64us(self):
        """returns digest as URL-safe, base64-encoded digest, as ASCII-encoded
        binary

        """
        return base64.urlsafe_b64encode(self)


    @staticmethod
    def from_base64us(d):
        """returns Digest object initialized from ASCII-encoded binary form of
        URL-safe, base64-encoded digest

        """

        return Digest(base64.urlsafe_b64decode(d))


if __name__ == "__main__":      # pragma: nocover
    import hashlib
    b = hashlib.sha512().digest()
    d = Digest(b)
    assert isinstance(d, Digest), "d isn't a Digest"
    d24 = d[:24]
    assert isinstance(d24, Digest), "d24 isn't a Digest"
    e = Digest.from_base64us(d.as_base64us())
    e24 = Digest.from_base64us(d24.as_base64us())
