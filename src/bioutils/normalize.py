# -*- coding: utf-8 -*-
"""normalize an allele on a reference sequence

A normalized Allele is:
* Minimum length: the Allele cannot be written as a shorter Allele
  that generates the same resulting sequence; and,
* Right-aligned: the Allele cannot be written as another Allele with
  numerically greater coordinates that would generate the same
  resulting sequence.

"""

import copy
import enum
import math
import logging

import attr


_logger = logging.getLogger(__name__)
debug = False

NormalizationMode = enum.Enum("NormalizationMode", "TRIMONLY LEFTSHUFFLE RIGHTSHUFFLE EXPAND VCF")


def normalize(
        sequence,
        interval,
        alleles,
        mode = NormalizationMode.EXPAND,
        bounds = None,
        anchor_length = 0,
        ):
    """normalizes the given alleles, returning a new interval and alleles.

    interval: location of alleles as (start, end) tuple, interbase coordiates
    alleles[]: array of sequence strings; first element is ref sequence and may be None
    bounds: maximal extent of normalization left and right; default to
      (0, len(sequence)); must be provided if sequence doesn't support __len__
    sequence: the reference sequence
    mode: shuffle left (à la vcf), shuffle right (hgvs), extend (voca)
    anchor: 0 (# of bases left and right)

    Returns: (interval, alleles) tuple

    """
    
    interval = Interval(*interval)
    if interval.start > interval.end:
        raise ValueError("Interval start > end; must be start <= end")

    if bounds is None:
        bounds = Interval(0,len(sequence))
    else:
        bounds = Interval(*bounds)

    left_anchor = right_anchor = anchor_length

    if not isinstance(mode, NormalizationMode):
        mode = NormalizationMode[mode]  # e.g., mode="LEFTSHUFFLE" OK

    if mode == NormalizationMode.VCF:
        if left_anchor:
            _logger.warning("Overriding requested left anchor size for VCF mode")
        left_anchor = 1
        mode = NormalizationMode.LEFTSHUFFLE

    if alleles[0] is not None:
        raise ValueError("First allele, the reference allele, must be None")
    alleles = list(alleles)     # in case tuple
    alleles[0] = sequence[interval.start:interval.end]
    if len(set(alleles)) < 2:
        raise ValueError("Must have at least two distinct alleles to normalize")
        
    if debug:
        _print_state(interval, bounds, sequence=sequence, alleles=alleles, comment="Starting state")

    # Trim: remove common suffix, prefix, and advance start
    l_trimmed, alleles = trim_left(alleles)
    interval.start += l_trimmed
    r_trimmed, alleles = trim_right(alleles)
    interval.end -= r_trimmed
    if debug:
        _print_state(interval, bounds, sequence=sequence, alleles=alleles, comment="After trimming")

    lens = [len(a) for a in alleles]

    # Left shuffle/expand
    if mode in (NormalizationMode.LEFTSHUFFLE, NormalizationMode.EXPAND):
        dist = roll_left(sequence, alleles, interval.start - 1, bounds.start)
        _logger.debug(f"Rolled left {dist}")
        if mode == NormalizationMode.LEFTSHUFFLE:
            for i, a in enumerate(alleles):
                if lens[i]:
                    cdist = -dist % lens[i]
                    alleles[i] = a[cdist:] + a[:cdist]
            interval.start -= dist
            interval.end -= dist
        elif mode == NormalizationMode.EXPAND:
            addlseq = sequence[interval.start-dist:interval.start]
            alleles = [addlseq + a for a in alleles]
            interval.start -= dist
        if debug:
            _print_state(interval, bounds, sequence=sequence, alleles=alleles, comment=f"After left extend {mode}")

    # Right shuffle/expand
    if mode in (NormalizationMode.RIGHTSHUFFLE, NormalizationMode.EXPAND):
        dist = roll_right(sequence, alleles, interval.end, bounds.end - 1)
        _logger.debug(f"Rolled right {dist}")
        if mode == NormalizationMode.RIGHTSHUFFLE:
            for i, a in enumerate(alleles):
                if lens[i]:
                    cdist = dist % lens[i]
                    alleles[i] = a[cdist:] + a[:cdist]
            interval.start += dist
            interval.end += dist
        elif mode == NormalizationMode.EXPAND:
            addlseq = sequence[interval.end:interval.end+dist]
            alleles = [a + addlseq for a in alleles]
            interval.end += dist
        if debug:
            _print_state(interval, bounds, sequence=sequence, alleles=alleles, comment=f"After right extend {mode}")

    # Add left and/or right flanking sequence
    if left_anchor or right_anchor:
        anchor_left = max(bounds.start, interval.start - left_anchor)
        anchor_right = min(bounds.end, interval.end + right_anchor)
        left_anchor_seq = sequence[anchor_left:interval.start]
        right_anchor_seq = sequence[interval.end:anchor_right]
        interval.start = anchor_left
        interval.end = anchor_right
        alleles = [left_anchor_seq + a + right_anchor_seq for a in alleles]
        if debug:
            _print_state(interval, bounds, sequence=sequence, alleles=alleles, comment="After anchoring")

    return (interval.start, interval.end), tuple(alleles)




############################################################################
# INTERNAL

@attr.s
class Interval:
    start = attr.ib(int)
    end = attr.ib(int)


pfx="        "



# TODO: Rewrite trim_* code -- no need to modify array each time!
def trim_left(alleles):
    """remove common prefix from left of all alleles, returning
    (number_trimmed, [new_alleles])

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
    """remove common suffix from right of all alleles, returning
    (number_trimmed, [new_alleles])

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
    """circularly permute ('roll') alleles left

    """
    lens = [len(a) for a in alleles]
    d = 0
    max_d = ref_pos - bound
    while (d <= max_d
           and not any(a and a[(d-1) % lens[i]] != sequence[ref_pos - d]
                       for i, a in enumerate(alleles))):
        d += 1
    return d


def roll_right(sequence, alleles, ref_pos, bound):
    """circularly permute ('roll') alleles right

    """
    lens = [len(a) for a in alleles]
    d = 0
    max_d = bound - ref_pos
    while (d <= max_d
           and not any(a and a[d % lens[i]] != sequence[ref_pos + d]
                       for i, a in enumerate(alleles))):
        d += 1
    return d


def _print_state(interval, bounds, sequence, alleles, comment):
    line = pfx + "  "*interval.start + "^"
    if interval.end > interval.start:
        line += "-" * ((interval.end - interval.start - 1) * 2 + 1) + "^"
    print(line + f" [{interval.start},{interval.end}): {alleles}   | {comment}")
    return

    margin = 15
    leftseq = sequence[max(0, interval.start-margin):interval.start]
    rightseq = sequence[interval.end:interval.end+margin]

    row_fmt = "{:>20.20s}{:>20.20s}{:^20.20s}{:<20.20s}{:<20.20s}"
    rows = [
        row_fmt.format(str(bounds.start), "", f"[{interval.start},{interval.end})", "", str(bounds.end)),
        row_fmt.format("//", "|", "", "|", "//"),
        row_fmt.format("", leftseq, alleles[0], rightseq, "")
    ] + [
        row_fmt.format("", "", a, "", "")
        for a in alleles[1:]
    ]
    print("\n".join(rows))


def _print_seq_row(sequence):
    print(pfx + " ".join("0123456789"*math.ceil(len(sequence)/10)))
    print(pfx + " " + " ".join(sequence))



if __name__ == "__main__":    # pragma: no cover
    from functools import partial

    if debug:
        _logger.setLevel("DEBUG")

    sequence = "CCCCCCCCACACACACACTAGCAGCAGCA"
    #    0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9
    #     C C C C C C C C A C A C A C A C A C T A G C A G C A G C A

    tests = [
        {"interval": (5,5), "alleles": [None, "C"]},
        {"interval": (5,6), "alleles": [None, "CC"]},
        {"interval": (5,6), "alleles": [None, ""]},
        {"interval": (13,13), "alleles": [None, "CA"]},
        {"interval": (14,14), "alleles": [None, "AC"]},
        {"interval": (13,15), "alleles": [None, ""]},
        {"interval": (22,22), "alleles": [None, "AGC"]},
        {"interval": (22,25), "alleles": [None, ""]},
        {"interval": (22,25), "alleles": [None, "", "AGC"]},
        {"interval": (22,25), "alleles": [None, "", "AGCAGC"]},
        ]

    normalize_seq = partial(normalize, sequence=sequence)
    normalize_trim = partial(normalize_seq, mode=NormalizationMode.TRIMONLY)
    normalize_left = partial(normalize_seq, mode=NormalizationMode.LEFTSHUFFLE)
    normalize_right = partial(normalize_seq, mode=NormalizationMode.RIGHTSHUFFLE)
    normalize_expand = partial(normalize_seq, mode=NormalizationMode.EXPAND)
    normalize_vcf = partial(normalize_seq, mode=NormalizationMode.VCF)

    
    def test1(**kwargs):
        print(f"* {kwargs}")
        _print_seq_row(sequence)
        result = normalize_seq(**kwargs)
        kwargs["mode"] = "NormalizationMode." + kwargs["mode"].name
        print(f"assert {result} == normalize_seq({kwargs})")

    # for test in tests:
    #     print("############################################################################")
    #     for mode in NormalizationMode:
    #         test1(mode=mode, **test)
              
    #    for bm in (None, 2):
    #        if bm is None:
    #            bounds = None
    #        else:
    #            bounds = (test["interval"][0]-bm, test["interval"][1]+bm)
