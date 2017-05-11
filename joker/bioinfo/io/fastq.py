#!/usr/bin/env python3
# encoding: utf-8

import os

from joker.bioinfo.io.fasta import SRecord


def read(infile):
    """
    Read a FASTQ file and yield dicts.

    Say we have

        @HWI-ST1426:113:HGTH7ADXX:1:1101:1358:2170 1:N:0:ATCACG
        ATATGAGGACAAACGATAATACCGCCGCCTTGGTTATCTAGGATCTCTTGAA...
        +
        BBBFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIFIIIIIIIIIIIIIIIIII...

    Then we will get

        {
            'cmt':  'HWI-ST1426:113:HGTH7ADXX:1:1101:1358:2170 1:N:0:ATCACG',
            'seq':  'ATATGAGGACAAACGATAATACCGCCGCCTTGGTTATCTAGGATCTCT...',
            'qual': 'BBBFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIFIIIIIIIIIIIIII...',
        }

    :param infile: a file object or path
    :return: a generator
    """
    if isinstance(infile, str):
        infile = open(infile)

    with infile:
        while True:
            cmt = infile.readline().strip()
            seq = infile.readline().strip()
            plus = infile.readline().strip()
            qual = infile.readline().strip()

            if not cmt:
                break
            if not cmt.startswith('@') or plus != '+':
                raise ValueError('fastq file <{}> is corrupted'.format(infile.path))
            yield SRecord(cmt=cmt[1:], seq=seq, qual=qual)


def read_one(infile):
    for rec in read(infile):
        return rec


def write(outfile, records, linesep=os.linesep):
    # `handle` is either a file object or a string

    if isinstance(outfile, str):
        outfile = open(outfile, 'w')

    with outfile:
        template = '@{cmt}{eol}{seq}{eol}+{qual}{eol}'
        for seqdict in records:
            block = template.format(eol=linesep, **seqdict)
            outfile.write(block)

# TODO
# <instrument>:<run number>:<flowcell ID>:<lane>:<tile>
# :<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<index sequence>
# http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm
