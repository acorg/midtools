#!/usr/bin/env python

from __future__ import print_function
import sys

from midtools.mutate import mutateRead
from midtools.utils import s

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

    parser.add_argument(
        '--idSuffix', default='',
        help=('Add this string to the end of the read ids. This is added '
              'after the string added by --editIds (if also used).'))

    parser.add_argument(
        '--editIds', action='store_true', default=False,
        help=('Add "-mutations:N" to the end of each read id, where N '
              'is the number of mutations introduced to the read.'))

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)
    rate = args.rate
    verbose = args.verbose
    editIds = args.editIds
    idSuffix = args.idSuffix
    format_ = 'fastq' if args.fastq else 'fasta'

    for read in reads:
        count = len(mutateRead(read, rate))
        if verbose:
            print('%d mutation%s made in read (len %d) %s' % (
                count, s(count), len(read), read.id), file=sys.stderr)
        read.id = (read.id +
                   (('-mutations:%d' % count) if editIds else '') +
                   idSuffix)
        print(read.toString(format_), end='')
