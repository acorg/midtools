#!/usr/bin/env python

from __future__ import print_function
import sys

from data.mutate import mutateRead

from dark.reads import (
    addFASTACommandLineOptions, parseFASTACommandLineOptions)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Mutate reads.')

    parser.add_argument(
        '--rate', type=float, required=True,
        help='The per-base mutation rate to use')

    parser.add_argument(
        '--verbose', action='store_true', default=False,
        help=('Print (to stderr) the number of mutations made to each '
              'sequence.'))

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)
    rate = args.rate
    verbose = args.verbose

    for read in reads:
        count = len(mutateRead(read, rate))
        if verbose:
            print('%d mutation%s made in read (len %d) %s' % (
                count, '' if count == 1 else 's', len(read), read.id),
                  file=sys.stderr)
        print(read.toString('fasta'), end='')
