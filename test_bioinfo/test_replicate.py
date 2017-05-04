#!/usr/bin/env python3
# coding: utf-8

from __future__ import division, print_function
from joker.bioinfo import revcompl


def test():
    fw = 'TTGATGGCTAAGAGTAAAATCTTAAAAAACACACTGGTTCTATATTTTCGTCAAGTTTTG'
    rc = 'CAAAACTTGACGAAAATATAGAACCAGTGTGTTTTTTAAGATTTTACTCTTAGCCATCAA'
    assert revcompl(fw, typ='dna') == rc


if __name__ == '__main__':
    test()
