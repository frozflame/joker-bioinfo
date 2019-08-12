#!/usr/bin/env python3
# coding: utf-8

from __future__ import division, print_function

from joker.bioinfo.io import fasta

p = '/Users/Hailong/Cloud/Imbark/code/jokerseries/joker-bioinfo/test_bioinfo/samp_single.fasta'
gen = fasta.read(p)
recs = list(gen)
