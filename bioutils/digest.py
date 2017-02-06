import base64

import six


@six.python_2_unicode_compatible
class Digest(six.binary_type):
    """Represents a binary digest, with support for encoding and decoding using printable characters

    Stringified representations use URL-safe base64 encodings. See
    https://tools.ietf.org/html/rfc4648#page-7

    >>> import hashlib

    >>> d = Digest(hashlib.sha512().digest())
    >>> len(d)
    64
    >>> str(d)
    'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXcg_SpIdNs6c5H0NE8XYXysP-DGNKHfuwvY7kxvUdBeoGlODJ6-SfaPg=='

    >>> d24 = d[:24]
    >>> len(d24)
    24
    >>> str(d24)
    'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc'

    >>> e = Digest.from_base64us(d.as_base64us())
    >>> e == d
    True

    >>> e24 = Digest.from_base64us(d24.as_base64us())
    >>> e24 == d24
    True

    """


    def __str__(self):
        """returns digest as URL-safe, base64 encoded digest, as Unicode"""
        return self.as_base64us().decode("ASCII")

    def __getitem__(self, key):
        return Digest(bytes(self)[key])

    def as_base64us(self):
        """returns digest as URL-safe, base64-encoded digest, as ASCII-encoded binary"""
        return base64.urlsafe_b64encode(self)

    @staticmethod
    def from_base64us(d):
        """returns Digest object initialized from ASCII-encoded binary form of URL-safe, base64-encoded digest"""
        return Digest(base64.urlsafe_b64decode(d))
