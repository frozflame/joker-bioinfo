#!/usr/bin/env python3
# coding: utf-8

from __future__ import division, print_function

from collections import Counter


ambig_nucl = {
    'M': 'AC',
    'R': 'AG',
    'W': 'AT',
    'S': 'CG',
    'Y': 'CT',
    'K': 'GT',
    'V': 'ACG',
    'H': 'ACT',
    'D': 'AGT',
    'B': 'CGT',
    'X': 'GATC',
    'N': 'GATC',
}


ambig_nucl_gc_equiv = {
    'B': 0.6666666666666666,
    'D': 0.3333333333333333,
    'H': 0.3333333333333333,
    'K': 0.5,
    'M': 0.5,
    'N': 0.5,
    'R': 0.5,
    'S': 1.0,
    'V': 0.6666666666666666,
    'W': 0.0,
    'X': 0.5,
    'Y': 0.5,
}


def calc_gc_content(seq, percent=False):
    counter = Counter(seq.upper())
    count = counter['G'] + counter['C'] + 0.

    for x in counter:
        if x in ambig_nucl_gc_equiv:
            count += ambig_nucl_gc_equiv[x] * counter[x]

    ratio = count / len(seq)

    if percent:
        return int(ratio * 100)
    else:
        return ratio


def calc_n50_statistic(lenths):
    lenths = list(lenths)
    half = sum(lenths) / 2.
    for x in sorted(lenths, reverse=True):
        half -= x
        if half <= 0:
            return x


def guess_sequence_type(seq):
    counter = Counter(seq.lower())
    if sum(counter[k] for k in 'atgcn') >= len(seq) * .6:
        return 'nucl'
    else:
        return 'prot'
