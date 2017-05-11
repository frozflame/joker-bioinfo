#!/usr/bin/env python3
# coding: utf-8

from __future__ import division, print_function

import collections
import itertools
import os
import re

from six import moves as six_moves


class SRecord(dict):
    """
    Friendly with interactive Python shell
    better displayed if it contains long string
    example:

        {'cmt': 'randSEQ', 'seq': 'TGCTTGGGGAATGTCT'~1000}

    """

    def __init__(self, *args, **kwargs):
        super(SRecord, self).__init__(*args, **kwargs)
        self.__dict__['attributes'] = {'cmt', 'seq', 'qual'}

    def __getattr__(self, key):
        if key in self.attributes:
            return self.get(key, None)
        else:
            # just to raise proper error
            return self.__getattribute__(key)

    def __setattr__(self, key, value):
        if key in self.attributes:
            self[key] = value
        else:
            self.__dict__[key] = value

    def __str__(self):
        title = '>{}'.format(self.cmt)
        return '{}\n{}'.format(title, self.seq)

    def __repr__(self):
        def fmt(s, length=32):
            # always use str.format to support non-string type s
            if len(repr(s)) <= length:
                return '{}'.format(repr(s))
            else:
                return '{}~{}'.format(repr(s[:length]), len(s))

        pairs = ((k, self[k]) for k in sorted(self))
        parts = ('{}: {}'.format(fmt(k), fmt(v)) for k, v in pairs)
        repr_ = '{{{}}}'.format(', '.join(parts))
        return repr_


def remove_whitespaces(string):
    return string[:0].join(string.split())


class FastAReader(object):
    def __init__(self, infile):
        self.cmt = ''
        self.count = 0
        self.queue = collections.deque()
        self.infile = infile

    def reset(self):
        self.cmt = ''
        self.count = 0
        self.queue.clear()
        self.infile.seek(0)

    def get_cmt(self):
        if self.cmt:
            return self.cmt
        else:
            return 'anonym.{:03}'.format(self.count)

    def set_cmt(self, cmtline):
        """
        Set self.cmt to new value, and return the old value 
        :param cmtline: a string like '>xxxxxx' 
        """
        cmt = self.get_cmt()
        # set state
        self.cmt = cmtline.lstrip().replace('>', '', 1)
        self.count += 1
        return cmt

    def pop_seq(self):
        seq = ''.join(self.queue)
        self.queue.clear()
        return seq

    def __iter__(self):
        self.reset()
        chunks = self.iterate_chunks()
        chunks = itertools.chain(chunks, ['>epilogue'])
        for chunk in chunks:
            if chunk.startswith('>'):
                seq = self.pop_seq()
                cmt = self.set_cmt(chunk)
                if seq or self.count > 1:
                    yield SRecord(cmt=cmt, seq=seq)
            else:
                self.queue.append(chunk)

    def iterate_chunks(self):
        """
        Yields strings.
        Each time, either 
            * pure sequence like 'AATCGCTTCGAGCAGGGCCC...CTGCC', or
            * cmt line like '>sp|O88942|NFAC1_MOUSE Nuclear ...'
        Leading and trailing spaces are removed.
        Spaces within sequences are removed.
        """
        cmtreg = re.compile(r'>[^>]*?\n')
        with self.infile:
            while True:
                chunk = self.infile.read(2 ** 20)
                headpos = 0
                if not chunk:
                    break
                for mat in cmtreg.finditer(chunk):
                    if mat.start() != headpos:
                        seq = chunk[headpos:mat.start()]
                        yield remove_whitespaces(seq)
                    yield mat.group().strip()
                    headpos = mat.end()
                seq = chunk[headpos:]
                yield remove_whitespaces(seq)


def read(infile):
    """
    Reading a FASTA file should NOT be complicated!

    :param infile: file-like object or string (file path)
    :return: a generator of SRecord objects

    Say we have

        >ORF00024
        ATCTGTCCTACTCCCGTC...TC
        >ORF00025
        GTCTGTCCTACTCCCGTC...TC

    Then we will get something equivolent to (but in generator form)

        [
            {
                'cmt': 'ORF00024',
                'seq': 'ATCTGTCCTACTCCCGTC...TC'
            },
            {
                'cmt': 'ORF00025',
                'seq': 'GTCTGTCCTACTCCCGTC...TC'
            }
        ]

    Nothing frustrates you. If you want to iterate through a multi-seq FASTA
    file a second time, `itertools.tee` may help you:

        recgen = fasta.read('contigs.fas')
        reciter, reciter1 = itertools.tee(recgen)

        for rec in reciter:
            print(rec.cmt, rec.seq)

        for rec in reciter1:
            print(rec.cmt, rec.seq)
    """
    if isinstance(infile, str):
        infile = open(infile)
    return iter(FastAReader(infile))


def read_one(infile):
    for rec in read(infile):
        return rec


def write(outfile, records, linesep=os.linesep, linewidth=70):
    """
    Reverse of `fasta.read`.

    :param outfile: a file-like object or outfile to the output FASTA file
    :param records: an iterable like [{'cmt': 'SEQ1', 'seq': 'ATCTC...T'}, ...]
    :param linesep: newline symbol
    :param linewidth: default 60
    :return: a generator
    """
    if isinstance(outfile, str):
        outfile = open(outfile, 'w')

    with outfile:
        # accept a single record
        if isinstance(records, dict):
            records = [records]

        for ix, rec in enumerate(records):
            cmt = rec.get('cmt')  # if empty, cstate will provide one
            seq = rec.get('seq')  # if empty, skip
            if not cmt:
                cmt = 'anonym.{:03}'.format(ix + 1)

            cmtline = '>{}{}'.format(cmt, linesep)
            outfile.write(cmtline)

            for i in six_moves.range(0, len(seq), linewidth):
                outfile.write(seq[i:i+linewidth])
                outfile.write(linesep)


def fix_cmt(cmt, prefix='', suffix=''):
    """
    Fix seqrecord comment

        >(prefix)Key(suffix) other_descriptions
        ATTCGGGGGTCTGGCTAG...
    """
    regex_key = re.compile(r'^\S+')
    return regex_key.sub(prefix + r'\g<0>' + suffix, cmt)


def fix_filename(name, prefix='', suffix=''):
    """
    Remove common fasta extensions and join with prefix & suffix

        (prefix)N2700.contigs<.fa>(suffix)
    """
    regex_ext = re.compile(r'\.(fa|fas|fasta|fna|ffn|faa|frn)$')
    return prefix + regex_ext.sub('', name) + suffix


def match_contig_name(names, keyword):
    keyword = re.escape(keyword)
    regex = re.compile(r'^(.*[^a-zA-Z0-9]+)?' + keyword + '([^a-zA-Z0-9]+.*)?$', re.I)
    for name in names:
        if regex.match(name):
            return name
