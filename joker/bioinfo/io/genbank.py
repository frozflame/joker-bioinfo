#!/usr/bin/env python3
# coding: utf-8

from __future__ import division, print_function

import textwrap

_sp = ' '
features_header = "FEATURES{}Location/Qualifiers".format(_sp * 12)


def fmt_location_line(feature_type, location):
    ft = feature_type
    line = _sp * 5 + ft + _sp * (16 - len(ft)) + location
    return [line]


def fmt_qualifier(key, value):
    splitted_lines = []
    whole_line = '/{}="{}"'.format(key, value)
    for line in textwrap.wrap(whole_line, 58):
        line = _sp * 21 + line
        splitted_lines.append(line)
    return splitted_lines


def fmt_feature(feature_type, location, qualifiers):
    lines = fmt_location_line(feature_type, location)
    for key, value in qualifiers.items():
        for line in fmt_qualifier(key, value):
            lines.append(line)
    return lines


def fmt_gene_cds(headpos, tailpos, qualifiers):
    gene_q = dict(gene=qualifiers.get('gene', '?'))

    if headpos < tailpos:
        location = '{}..{}'.format(headpos, tailpos)
    else:
        location = 'complement({}..{})'.format(tailpos, headpos)

    lines = fmt_feature('gene', location, gene_q)
    lines.extend(fmt_feature('CDS', location, qualifiers))
    return lines


def fmt_origin_sequence(seq):
    lines = ['ORIGIN']
    for i in range(0, len(seq), 60):
        line = '{0:>9} '.format(i + 1) + _sp.join(textwrap.wrap(seq[i:i + 60], 10))
        lines.append(line)
    lines.extend(['//', ''])
    return lines
