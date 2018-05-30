#!/usr/bin/env python

from __future__ import division, print_function

import sys
from collections import Counter
from math import log10

from data.data import addCommandLineOptions, parseCommandLineOptions


def baseCountsToStr(counts):
    """
    @param counts: A C{counter} instance.
    """
    total = sum(counts.values())
    return ', '.join([
        ('%s %d (%.2f%%)' % (base, counts[base], counts[base] / total * 100.0))
        for base in sorted(counts)])


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Analyze a set of aligned reads.')

    addCommandLineOptions(parser, 'significant-base-frequencies.html')

    parser.add_argument(
        '--verbose', action='store_true', default=False,
        help='Print verbose textual output showing read connections.')

    args = parser.parse_args()

    (genome, alignedReads, readCountAtLocation, baseCountAtLocation,
     readsAtLocation, _) = parseCommandLineOptions(args, False)

    print('Read %d aligned reads.' % len(alignedReads), file=sys.stderr)

    genomeLength = len(genome)
    genomeLengthWidth = int(log10(genomeLength)) + 1
    nucleotides = set('ACGT')

    for offset in range(genomeLength):
        counts = Counter()
        for read in readsAtLocation[offset]:
            base = read.base(offset)
            if base in nucleotides:
                counts[base] += 1
        print('Location %*d: base counts %s' % (
            genomeLengthWidth, offset + 1, baseCountsToStr(counts)))
