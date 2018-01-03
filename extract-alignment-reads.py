#!/usr/bin/env python

from __future__ import division, print_function
import sys

from dark.reads import (
    addFASTACommandLineOptions, parseFASTACommandLineOptions)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=('Extract reads from an alignment, ignoring the first '
                     'consensus sequence.'))

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)

    for read in reads:
        sequence = read.sequence

        # Scan for initial gaps.
        offset = 0
        for base in sequence:
            if base == '-':
                offset += 1
            else:
                break

        # Scan for final gaps.
        trailing = 0
        for base in sequence[::-1]:
            if base == '-':
                trailing += 1
            else:
                break

        # Make sure the read is not all gaps.
        assert offset + trailing < len(sequence)

        read.sequence = sequence[offset:len(sequence) - trailing]

        print(read.toString('fasta'), end='')
