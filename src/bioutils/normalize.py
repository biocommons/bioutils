# -*- coding: utf-8 -*-
"""Provides functionality for normalizing alleles, ensuring comparable representations.

"""

import copy
import enum
import logging
import math
from typing import Optional

import attr

_logger = logging.getLogger(__name__)
debug = False

NormalizationMode = enum.Enum(
    "NormalizationMode", "EXPAND LEFTSHUFFLE RIGHTSHUFFLE TRIMONLY VCF"
)
"""Enum passed to normalize to select the normalization mode.

Attributes:
    EXPAND: Normalize alleles to maximal extent both left and right.
    LEFTSHUFFLE: Normalize alleles to maximal extent left.
    RIGHTSHUFFLE: Normalize alleles to maximal extent right.
    TRIMONLY: Only trim the common prefix and suffix of alleles. Deprecated -- use `mode=None` with `trim=True` instead.
    VCF: Normalize with VCF. 
"""


def normalize(
    sequence,
    interval,
    alleles,
    mode: Optional[NormalizationMode] = NormalizationMode.EXPAND,
    bounds=None,
    anchor_length=0,
    trim: bool = True,
):
    """Normalizes the alleles that co-occur on sequence at interval, ensuring comparable representations.

    Normalization performs three operations:
    - trimming
    - shuffling
    - anchoring

    Args:
        sequence (str or iterable): The reference sequence; must support indexing and ``__getitem__``.
        interval (2-tuple of int): The location of alleles in the reference sequence as ``(start, end)``.
            Interbase coordinates.
        alleles (iterable of str): The sequences to be normalized. The first element
            corresponds to the reference sequence being unchanged and must be None.
        bounds (2-tuple of int, optional): Maximal extent of normalization left and right.
            Must be provided if sequence doesn't support ``__len__``. Defaults to ``(0, len(sequence))``.
        mode (NormalizationMode Enum or string, optional): A NormalizationMode Enum or the corresponding string.
            Defaults to ``EXPAND``. Set to None to skip shuffling. Does not affect trimming or anchoring.
        anchor (int, optional): number of flanking residues left and right. Defaults to ``0``.
        trim (bool): indicates whether to trim the common prefix and suffix of alleles. Defaults to True.
            Set to False to skip trimming. Does not affect shuffling or anchoring.

    Returns:
        tuple: ``(new_interval, [new_alleles])``

    Raises:
        ValueError: If normalization mode is VCF and `anchor_length` is nonzero.
        ValueError: If the interval start is greater than the end.
        ValueError: If the first (reference) allele is not `None`.
        ValueError: If there are not at least two distinct alleles.

    Examples:
        >>> sequence = "CCCCCCCCACACACACACTAGCAGCAGCA"
        >>> normalize(sequence, interval=(22,25), alleles=(None, "GC", "AGCAC"), mode='TRIMONLY')
        ((22, 24), ('AG', 'G', 'AGCA'))

        >>> normalize(sequence, interval=(22, 22), alleles=(None, 'AGC'), mode='RIGHTSHUFFLE')
        ((29, 29), ('', 'GCA'))

        >>> normalize(sequence, interval=(22, 22), alleles=(None, 'AGC'), mode='EXPAND')
        ((19, 29), ('AGCAGCAGCA', 'AGCAGCAGCAGCA'))
    """

    interval = _Interval(*interval)
    if interval.start > interval.end:
        raise ValueError("Interval start > end; must be start <= end")

    if bounds is None:
        bounds = _Interval(0, len(sequence))
    else:
        bounds = _Interval(*bounds)

    left_anchor = right_anchor = anchor_length

    if mode is not None and not isinstance(mode, NormalizationMode):
        mode = NormalizationMode[mode]  # e.g., mode="LEFTSHUFFLE" OK

    if mode == NormalizationMode.VCF:
        if anchor_length:
            raise ValueError(
                "May not provide non-zero anchor size with VCF normalization mode"
            )
        if not trim:
            raise ValueError("May not disable trimming with VCF normalization mode")
        mode = NormalizationMode.LEFTSHUFFLE
        left_anchor = 1
        right_anchor = 0

    if alleles[0] is not None:
        raise ValueError("First allele, the reference allele, must be None")
    alleles = list(alleles)  # in case tuple
    alleles[0] = sequence[interval.start : interval.end]

    if debug:
        _print_state(
            interval,
            bounds,
            sequence=sequence,
            alleles=alleles,
            comment="Starting state",
        )

    if trim:
        if len(set(alleles)) < 2:
            raise ValueError("Must have at least two distinct alleles to trim")

        # Trim: remove common suffix, prefix, and adjust interval to match
        l_trimmed, alleles = trim_left(alleles)
        interval.start += l_trimmed
        r_trimmed, alleles = trim_right(alleles)
        interval.end -= r_trimmed
        if debug:
            _print_state(
                interval,
                bounds,
                sequence=sequence,
                alleles=alleles,
                comment="After trimming",
            )

    lens = [len(a) for a in alleles]

    if mode == NormalizationMode.LEFTSHUFFLE:
        dist = roll_left(sequence, alleles, interval.start - 1, bounds.start)
        for i, a in enumerate(alleles):
            if lens[i]:
                adist = -dist % lens[i]
                alleles[i] = a[adist:] + a[:adist]
        interval.start -= dist
        interval.end -= dist

    elif mode == NormalizationMode.RIGHTSHUFFLE:
        dist = roll_right(sequence, alleles, interval.end, bounds.end - 1)
        for i, a in enumerate(alleles):
            if lens[i]:
                adist = dist % lens[i]
                alleles[i] = a[adist:] + a[:adist]
        interval.start += dist
        interval.end += dist

    elif mode == NormalizationMode.EXPAND:
        ldist = roll_left(sequence, alleles, interval.start - 1, bounds.start)
        rdist = roll_right(sequence, alleles, interval.end, bounds.end - 1)

        lseq = sequence[interval.start - ldist : interval.start]
        rseq = sequence[interval.end : interval.end + rdist]
        alleles = [lseq + a + rseq for a in alleles]

        interval.start -= ldist
        interval.end += rdist

    if debug:
        _print_state(
            interval,
            bounds,
            sequence=sequence,
            alleles=alleles,
            comment=f"After mode: {mode}",
        )

    # Add left and/or right flanking sequence
    if left_anchor or right_anchor:
        anchor_left = max(bounds.start, interval.start - left_anchor)
        anchor_right = min(bounds.end, interval.end + right_anchor)
        left_anchor_seq = sequence[anchor_left : interval.start]
        right_anchor_seq = sequence[interval.end : anchor_right]
        interval.start = anchor_left
        interval.end = anchor_right
        alleles = [left_anchor_seq + a + right_anchor_seq for a in alleles]
        if debug:
            _print_state(
                interval,
                bounds,
                sequence=sequence,
                alleles=alleles,
                comment="After anchoring",
            )

    return (interval.start, interval.end), tuple(alleles)


############################################################################
# INTERNAL


@attr.s
class _Interval:
    start = attr.ib(int)
    end = attr.ib(int)


# TODO: Rewrite trim_* code -- no need to modify array each time!
def trim_left(alleles):
    """Removes common prefix of given alleles.

    Args:
        alleles (list of str): A list of alleles.

    Returns:
        tuple: ``(number_trimmed, [new_alleles])``.

    Examples:
        >>> trim_left(["","AA"])
        (0, ['', 'AA'])

        >>> trim_left(["A","AA"])
        (1, ['', 'A'])

        >>> trim_left(["AT","AA"])
        (1, ['T', 'A'])

        >>> trim_left(["AA","AA"])
        (2, ['', ''])

        >>> trim_left(["CAG","CG"])
        (1, ['AG', 'G'])
    """

    trimmed = 0
    while all(len(a) > 0 for a in alleles):
        a0 = alleles[0]
        for a in alleles[1:]:
            if a0[0] != a[0]:
                return trimmed, alleles
        alleles = [a[1:] for a in alleles]
        trimmed += 1
    return (trimmed, alleles)


def trim_right(alleles):
    """Removes common suffix of given alleles.

    Args:
        alleles (list of str): A list of alleles.

    Returns:
        tuple: ``(number_trimmed, [new_alleles])``.

    Examples:
        >>> trim_right(["","AA"])
        (0, ['', 'AA'])

        >>> trim_right(["A","AA"])
        (1, ['', 'A'])

        >>> trim_right(["AT","AA"])
        (0, ['AT', 'AA'])

        >>> trim_right(["AA","AA"])
        (2, ['', ''])

        >>> trim_right(["CAG","CG"])
        (1, ['CA', 'C'])
    """

    trimmed = 0
    while all(len(a) > 0 for a in alleles):
        a0 = alleles[0]
        for a in alleles[1:]:
            if a0[-1] != a[-1]:
                return trimmed, alleles
        alleles = [a[:-1] for a in alleles]
        trimmed += 1
    return (trimmed, alleles)


def roll_left(sequence, alleles, ref_pos, bound):
    """Determines common distance all alleles can be rolled (circularly permuted) left
    within the reference sequence without altering it.

    Args:
        sequence (str): The reference sequence.
        alleles (list of str): The sequences to be normalized.
        ref_pos (int): The beginning index for rolling.
        bound (int): The lower bound index in the reference sequence for normalization, hence also for rolling.

    Returns:
        int: The distance that the alleles can be rolled.
    """

    # circularly permute sequence d steps, using modulo arithmetic
    lens = [len(a) for a in alleles]
    d = 0
    max_d = ref_pos - bound
    while d <= max_d and not any(
        a and a[-(d + 1) % lens[i]] != sequence[ref_pos - d]
        for i, a in enumerate(alleles)
    ):
        d += 1
    return d


def roll_right(sequence, alleles, ref_pos, bound):
    """Determines common distance all alleles can be rolled (circularly permuted) right
    within the reference sequence without altering it.

    Args:
        sequence (str): The reference sequence.
        alleles (list of str): The sequences to be normalized.
        ref_pos (int): The start index for rolling.
        bound (int): The upper bound index in the reference sequence for normalization, hence also for rolling.
    Returns:
        int: The distance that the alleles can be rolled
    """

    # circularly permute sequence d steps, using modulo arithmetic
    lens = [len(a) for a in alleles]
    d = 0
    max_d = bound - ref_pos
    while d <= max_d and not any(
        a and a[d % lens[i]] != sequence[ref_pos + d] for i, a in enumerate(alleles)
    ):
        d += 1
    return d


############################################################################
# Debugging

pfx = "        "


def _print_state(interval, bounds, sequence, alleles, comment):
    """ """
    line = pfx + "  " * interval.start + "^"
    if interval.end > interval.start:
        line += "-" * ((interval.end - interval.start - 1) * 2 + 1) + "^"
    print(line + f" [{interval.start},{interval.end}): {alleles}   | {comment}")
    return

    margin = 15
    leftseq = sequence[max(0, interval.start - margin) : interval.start]
    rightseq = sequence[interval.end : interval.end + margin]

    row_fmt = "{:>20.20s}{:>20.20s}{:^20.20s}{:<20.20s}{:<20.20s}"
    rows = [
        row_fmt.format(
            str(bounds.start),
            "",
            f"[{interval.start},{interval.end})",
            "",
            str(bounds.end),
        ),
        row_fmt.format("//", "|", "", "|", "//"),
        row_fmt.format("", leftseq, alleles[0], rightseq, ""),
    ] + [row_fmt.format("", "", a, "", "") for a in alleles[1:]]
    print("\n".join(rows))


def _print_seq_row(sequence):
    print(pfx + " ".join("0123456789" * math.ceil(len(sequence) / 10)))
    print(pfx + " " + " ".join(sequence))


if __name__ == "__main__":  # pragma: no cover
    from functools import partial

    if debug:
        _logger.setLevel("DEBUG")

    sequence = "CCCCCCCCACACACACACTAGCAGCAGCAT"
    #                        1                   2                   3
    #    0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0
    #     C C C C C C C C A C A C A C A C A C T A G C A G C A G C A T

    tests = [
        # {"interval": (5,5), "alleles": [None, "C"]},
        # {"interval": (5,6), "alleles": [None, "CC"]},
        # {"interval": (5,6), "alleles": [None, ""]},
        # {"interval": (13,13), "alleles": [None, "CA"]},
        # {"interval": (14,14), "alleles": [None, "AC"]},
        # {"interval": (13,15), "alleles": [None, ""]},
        {"interval": (22, 22), "alleles": [None, "AGC"]},
        {"interval": (22, 22), "alleles": [None, "AGCT"]},
        {"interval": (22, 22), "alleles": [None, "AGC", "AGCT"]},
        # {"interval": (22,25), "alleles": [None, ""]},
        # {"interval": (22,25), "alleles": [None, "", "AGC"]},
        # {"interval": (22,25), "alleles": [None, "", "AGCAGC"]},
    ]

    normalize_seq = partial(normalize, sequence=sequence)
    normalize_trim = partial(normalize_seq, mode=NormalizationMode.TRIMONLY)
    normalize_left = partial(normalize_seq, mode=NormalizationMode.LEFTSHUFFLE)
    normalize_right = partial(normalize_seq, mode=NormalizationMode.RIGHTSHUFFLE)
    normalize_expand = partial(normalize_seq, mode=NormalizationMode.EXPAND)
    normalize_vcf = partial(normalize_seq, mode=NormalizationMode.VCF)

    debug = True

    def test1(**kwargs):
        print(f"* {kwargs}")
        _print_seq_row(sequence)
        result = normalize_seq(**kwargs)
        kwargs["mode"] = str(kwargs["mode"])
        print(f"assert {result} == normalize_seq({kwargs})")

    for test in tests:
        print(
            "############################################################################"
        )
        for mode in ("EXPAND",):  # "LEFTSHUFFLE", "RIGHTSHUFFLE", "EXPAND"):
            for bm in (None,):
                if bm is None:
                    bounds = None
                else:
                    bounds = (test["interval"][0] - bm, test["interval"][1] + bm)
                test["bounds"] = bounds
                test1(mode=mode, **test)
