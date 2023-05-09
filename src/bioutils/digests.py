# -*- coding: utf-8 -*-

import base64
import hashlib

from .sequences import normalize_sequence
from .vmc_digest import vmc_digest


def seq_seqhash(seq, normalize=True):
    """Converts sequence to 24-byte Truncated Digest.

    Args:
        seq (str): A sequence.
        normalize (bool, optional): Whether to normalize the sequence before conversion,
            i.e. to ensure representation as uppercase letters without whitespace or asterisks.
            Defaults to ``True``.

    Returns:
        str: 24-byte Truncated Digest representation of sequence.

    Examples:
        >>> seq_seqhash("")
        'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc'

        >>> seq_seqhash("ACGT")
        'aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2'

        >>> seq_seqhash("acgt")
        'aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2'

        >>> seq_seqhash("acgt", normalize=False)
        'eFwawHHdibaZBDcs9kW3gm31h1NNJcQe'
    """

    seq = normalize_sequence(seq) if normalize else seq
    return str(vmc_digest(seq, digest_size=24))


def seq_seguid(seq, normalize=True):
    """Converts sequence to seguid.

    This seguid is compatible with BioPython's seguid.

    Args:
        seq (str): A sequence.
        normalize (bool, optional): Whether to normalize the sequence before conversion,
            i.e. to ensure representation as uppercase letters without whitespace or asterisks.
            Defaults to ``True``.

    Returns:
        str: seguid representation of sequence.

    Examples:
        >>> seq_seguid('')
        '2jmj7l5rSw0yVb/vlWAYkK/YBwk'

        >>> seq_seguid('ACGT')
        'IQiZThf2zKn/I1KtqStlEdsHYDQ'

        >>> seq_seguid('acgt')
        'IQiZThf2zKn/I1KtqStlEdsHYDQ'

        >>> seq_seguid('acgt', normalize=False)
        'lII0AoG1/I8qKY271rgv5CFZtsU'
    """

    seq = normalize_sequence(seq) if normalize else seq
    bseq = seq.encode("ascii")
    return base64.b64encode(hashlib.sha1(bseq).digest()).decode("ascii").rstrip("=")


def seq_md5(seq, normalize=True):
    """Converts sequence to unicode md5 hex digest.

    Args:
        seq (str): A sequence.
        normalize (bool, optional): Whether to normalize the sequence before conversion,
            i.e. to ensure representation as uppercase letters without whitespace or asterisks.
            Defaults to ``True``.

    Returns:
        str: Unicode md5 hex digest representation of sequence.

    Examples:
        >>> seq_md5('')
        'd41d8cd98f00b204e9800998ecf8427e'

        >>> seq_md5('ACGT')
        'f1f8f4bf413b16ad135722aa4591043e'

        >>> seq_md5('ACGT*')
        'f1f8f4bf413b16ad135722aa4591043e'

        >>> seq_md5(' A C G T ')
        'f1f8f4bf413b16ad135722aa4591043e'

        >>> seq_md5('acgt')
        'f1f8f4bf413b16ad135722aa4591043e'

        >>> seq_md5('acgt', normalize=False)
        'db516c3913e179338b162b2476d1c23f'
    """

    seq = normalize_sequence(seq) if normalize else seq
    bseq = seq.encode("ascii")
    return hashlib.md5(bseq).hexdigest()


def seq_sha1(seq, normalize=True):
    """Converts sequence to unicode sha1 hexdigest.

    Args:
        seq (str): A sequence.
        normalize (bool, optional): Whether to normalize the sequence before conversion,
            i.e. to ensure representation as uppercase letters without whitespace or asterisks before encoding.
            Defaults to ``True``.

    Returns:
        str: Unicode sha1 hexdigest representation of sequence.

    Examples:
        >>> seq_sha1('')
        'da39a3ee5e6b4b0d3255bfef95601890afd80709'

        >>> seq_sha1('ACGT')
        '2108994e17f6cca9ff2352ada92b6511db076034'

        >>> seq_sha1('acgt')
        '2108994e17f6cca9ff2352ada92b6511db076034'

        >>> seq_sha1('acgt', normalize=False)
        '9482340281b5fc8f2a298dbbd6b82fe42159b6c5'
    """

    seq = normalize_sequence(seq) if normalize else seq
    bseq = seq.encode("ascii")
    return hashlib.sha1(bseq).hexdigest()


def seq_sha512(seq, normalize=True):
    """Converts sequence to unicode sha512 hexdigest.

    Args:
        seq (str): A sequence.
        normalize (bool, optional): Whether to normalize the sequence before conversion,
            i.e. to ensure representation as uppercase letters without whitespace or asterisks.
            Defaults to ``True``.

    Returns:
        str: Unicode sha512 hexdigest representation of sequence.

    Examples:
        >>> seq_sha512('')
        'cf83e1357eefb8bdf1542850d66d8007d620e4050b5715dc83f4a921d36ce9ce47d0d13c5d85f2b0ff8318d2877eec2f63b931bd47417a81a538327af927da3e'

        >>> seq_sha512('ACGT')
        '68a178f7c740c5c240aa67ba41843b119d3bf9f8b0f0ac36cf701d26672964efbd536d197f51ce634fc70634d1eefe575bec34c83247abc52010f6e2bbdb8253'

        >>> seq_sha512('acgt')
        '68a178f7c740c5c240aa67ba41843b119d3bf9f8b0f0ac36cf701d26672964efbd536d197f51ce634fc70634d1eefe575bec34c83247abc52010f6e2bbdb8253'

        >>> seq_sha512('acgt', normalize=False)
        '785c1ac071dd89b69904372cf645b7826df587534d25c41edb2862e54fb2940d697218f2883d2bf1a11cdaee658c7f7ab945a1cfd08eb26cbce57ee88790250a'
    """

    seq = normalize_sequence(seq) if normalize else seq
    bseq = seq.encode("ascii")
    return hashlib.sha512(bseq).hexdigest()


def seq_vmc_id(seq, normalize=True):
    """Converts sequence to VMC id.

    See https://github.com/ga4gh/vmc

    Args:
        seq (str): A sequence.
        normalize (bool, optional): Whether to normalize the sequence before conversion,
            i.e. to ensure representation as uppercase letters without whitespace or asterisks.
            Defaults to ``True``.

    Returns:
        str: VMC id representation of sequence.

    Examples:
        >>> seq_vmc_id("")
        'VMC:GS_z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc'

        >>> seq_vmc_id("ACGT")
        'VMC:GS_aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2'

        >>> seq_vmc_id("acgt")
        'VMC:GS_aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2'

        >>> seq_vmc_id("acgt", normalize=False)
        'VMC:GS_eFwawHHdibaZBDcs9kW3gm31h1NNJcQe'
    """

    return "{ir[namespace]}:{ir[accession]}".format(ir=seq_vmc_identifier(seq, normalize))


def seq_vmc_identifier(seq, normalize=True):
    """Converts sequence to VMC identifier (record).

    See https://github.com/ga4gh/vmc

    Args:
        seq (str): A sequence.
        normalize (bool, optional): Whether to normalize the sequence before conversion,
            i.e. to ensure representation as uppercase letters without whitespace or asterisks.
            Defaults to ``True``.

    Returns:
        str: VMC identifier (record) representation of sequnce.

    Examples:
        >>> seq_vmc_identifier("") == {'namespace': 'VMC', 'accession': 'GS_z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc'}
        True

        >>> seq_vmc_identifier("ACGT") == {'namespace': 'VMC', 'accession': 'GS_aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2'}
        True

        >>> seq_vmc_identifier("acgt") == {'namespace': 'VMC', 'accession': 'GS_aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2'}
        True

        >>> seq_vmc_identifier("acgt", normalize=False) == {'namespace': 'VMC', 'accession': 'GS_eFwawHHdibaZBDcs9kW3gm31h1NNJcQe'}
        True
    """

    seq = normalize_sequence(seq) if normalize else seq
    return {"namespace": "VMC", "accession": "GS_" + str(vmc_digest(seq))}
