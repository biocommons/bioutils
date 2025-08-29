"""Sliceable bytes subclass for binary digests with base64/base64url/hex encode/decode"""

import base64
import binascii

_enc = "ascii"


class Digest(bytes):
    """Represent a sliceable binary digest, with support for encoding and decoding using printable characters.

    Supported encoding and decodings are::
        * base64
        * base64url
        * hex (aka base16)

    The Base64 specification
    (https://tools.ietf.org/html/rfc4648#page-7) defines base64 and a
    URL-safe variant called base64url.

    "Stringified" Digest objects use URL-safe base64 encodings.


    >>> import hashlib

    >>> b = hashlib.sha512().digest()
    >>> len(b)
    64


    >>> d = Digest(b)  # creation
    >>> str(d)  # returns base64url
    'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXcg_SpIdNs6c5H0NE8XYXysP-DGNKHfuwvY7kxvUdBeoGlODJ6-SfaPg=='

    >>> d24 = d[:24]  # slice binary digest at first 24 bytes
    >>> str(d24)
    'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc'

    # encoding

    >>> d.as_base64url()
    'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXcg_SpIdNs6c5H0NE8XYXysP-DGNKHfuwvY7kxvUdBeoGlODJ6-SfaPg=='
    >>> d.as_hex()
    'cf83e1357eefb8bdf1542850d66d8007d620e4050b5715dc83f4a921d36ce9ce47d0d13c5d85f2b0ff8318d2877eec2f63b931bd47417a81a538327af927da3e'

    # decoding

    >>> d == Digest.from_base64(d.as_base64())
    True
    >>> d == Digest.from_base64url(d.as_base64url())
    True
    >>> d == Digest.from_hex(d.as_hex())
    True
    """

    def __str__(self) -> str:
        """Return digest as base64url string"""
        return self.as_base64url()

    # TODO: Consider requiring slice start == None or 0, and len % 3 == 0
    # Slicing %3 != 0 => strings will having suffix differences
    def __getitem__(self, key):  # noqa: ANN001 ANN204
        return Digest(bytes.__getitem__(self, key))

    # base64
    def as_base64(self) -> str:
        """Return Digest as a base64-encoded string.

        Returns:
            str: base64 encoding of Digest.

        """
        return base64.b64encode(self).decode(_enc)

    @staticmethod
    def from_base64(s: str):  # noqa: ANN205
        """Return Digest object initialized from a base64-encoded string.

        Args:
            s (str): A base64-encoded digest string.

        Returns:
            Digest: A Digest object initialized from s.

        """
        return Digest(base64.b64decode(s))

    # base64url
    def as_base64url(self) -> str:
        """Return Digest as URL-safe, base64-encoded string.

        Returns:
            str: URL-safe base64 encoding of Digest.

        """
        return base64.urlsafe_b64encode(self).decode(_enc)

    @staticmethod
    def from_base64url(s: str):  # noqa: ANN205
        """Returns Digest object initialized from a base64url string.

        Args:
            s (str): A base64url-encoded digest string.

        Returns:
            Digest: A Digest object initialized from s.

        """
        return Digest(base64.urlsafe_b64decode(s))

    # for backward compatibility with earlier versions
    # ("base64url" is the official name for the encoding)
    as_base64us = as_base64url
    from_base64us = from_base64url

    # hex
    def as_hex(self) -> str:
        """Returns Digest as hex string.

        Returns:
            str: A hex-encoding of Digest.

        """
        return binascii.hexlify(self).decode(_enc)

    @staticmethod
    def from_hex(s: str):  # noqa: ANN205
        """Returns Digest object initialized from hex string.

        Args:
            s (str): A hex-encoded digest string.

        Returns:
            Digest: A Digest object initialized from s.

        """
        return Digest(binascii.unhexlify(s))


if __name__ == "__main__":  # pragma: nocover
    import hashlib

    b = hashlib.sha512().digest()
    d = Digest(b)
    assert isinstance(d, Digest), "d isn't a Digest"  # noqa: S101
    d24 = d[:24]
    assert isinstance(d24, Digest), "d24 isn't a Digest"  # noqa: S101
    e = Digest.from_base64url(d.as_base64url())
    e24 = Digest.from_base64url(d24.as_base64url())
