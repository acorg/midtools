#!/usr/bin/env python

from __future__ import division, print_function

from dark.reads import (
    addFASTACommandLineOptions, parseFASTACommandLineOptions)

from data.data import AlignedRead


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Extract reads from an alignment')

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)

    for read in reads:
        ar = AlignedRead(read)
        print(ar.read.toString('fasta'), end='')
