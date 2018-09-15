#!/usr/bin/env python3
# coding: utf-8

from __future__ import division, print_function
from joker.bioinfo import revcompl, replication
from joker.bioinfo.replication import Transcoder


def test():
    fw = 'TTGATGGCTAAGAGTAAAATCTTAAAAAACACACTGGTTCTATATTTTCGTCAAGTTTTG'
    rc = 'CAAAACTTGACGAAAATATAGAACCAGTGTGTTTTTTAAGATTTTACTCTTAGCCATCAA'
    assert revcompl(fw, typ='dna') == rc
    predefined = [
        replication.compl_dna,
        replication.compl_rna,
        replication.revcompl_dna,
        replication.revcompl_rna,
    ]
    for p in predefined:
        r = Transcoder(**p).transcode(p['src']).decode()
        print(r)
        print(p)
        assert p['dst'] in {r, r[::-1]}


if __name__ == '__main__':
    test()
