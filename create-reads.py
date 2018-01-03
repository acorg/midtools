#!/usr/bin/env python

from __future__ import print_function

import sys
from random import uniform, normalvariate
from math import log10

from dark.reads import Read
from data.mutate import mutateRead

from dark.reads import (
    addFASTACommandLineOptions, parseFASTACommandLineOptions)


def makeRead(genome, genomeLen, meanLength, sdLength,
             minReadLength, maxReadLength, id_, rate, circularGenome):
    """
    Make a read. Assumes a circular genome!
    """
    length = genomeLen + 1

    while (length > genomeLen or length <= 0 or
           length < minReadLength or length > maxReadLength):
        length = int(normalvariate(meanLength, sdLength) + 0.5)

    if circularGenome:
        offset = int(uniform(0.0, genomeLen))

        sequence = genome[offset:offset + length]

        # If we didn't get enough from the end of the genome, take whatever
        # else we need from its start.
        if len(sequence) < length:
            sequence += genome[0:length - len(sequence)]

        assert len(sequence) == length
    else:
        offset = int(uniform(0.0, genomeLen - length))
        sequence = genome[offset:offset + length]
        assert len(sequence) == length

    read = Read(id_, sequence)
    mutationOffsets = mutateRead(read, rate)
    return read, offset, mutationOffsets


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Mutate reads.')

    parser.add_argument(
        '--idPrefix', default='read-',
        help=('The prefix for the created read ids. The read number '
              'will be appended.'))

    parser.add_argument(
        '--count', type=int, default=100,
        help='The number of reads to create')

    parser.add_argument(
        '--minReadLength', type=int, default=10,
        help='The minimum length read to create')

    parser.add_argument(
        '--maxReadLength', type=int, default=None,
        help=('The maximum length read to create. Defaults to the genome '
              'length'))

    parser.add_argument(
        '--rate', type=float, default=0.0,
        help='The per-base mutation rate to use')

    parser.add_argument(
        '--meanLength', type=float, default=100.0,
        help='The mean read length')

    parser.add_argument(
        '--sdLength', type=float, default=10.0,
        help='The standard deviation of read length')

    parser.add_argument(
        '--verbose', action='store_true', default=False,
        help='Print (to stderr) information about the created reads.')

    parser.add_argument(
        '--circularGenome', action='store_true', default=False,
        help=('If specified, reads will wrap around the genome (currently not '
              'compatible with --alignReads).'))

    parser.add_argument(
        '--printGenome', action='store_true', default=False,
        help='If specified, print the genome as the first sequence.')

    parser.add_argument(
        '--alignReads', action='store_true', default=False,
        help='If specified, print the reads aligned to the genome.')

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = list(parseFASTACommandLineOptions(args))
    # There should only be one "read", the sequence we are to create other
    # reads from.
    assert len(reads) == 1, (
        'FASTA input contained %d sequence%s (expected just one).' % (
            len(reads), '' if len(reads) == 1 else 's'))
    genome = reads[0]
    genomeLen = len(genome)
    meanLength = args.meanLength

    if meanLength > genomeLen:
        raise ValueError('The mean read length (%d) is greater than the '
                         'genome length (%d)' % (int(meanLength), genomeLen))

    if meanLength <= 0:
        raise ValueError('The mean read length must be greater than zero')

    sdLength = args.sdLength

    if sdLength <= 0.0:
        raise ValueError('The read length standard deviation must be > 0.0')

    rate = args.rate

    if not (0.0 <= rate <= 1.0):
        raise ValueError('The read mutation rate must be in [0.0, 1.0]')

    minReadLength = args.minReadLength

    if minReadLength <= 0:
        raise ValueError('The minimum read length must be > 0')

    maxReadLength = args.maxReadLength

    if maxReadLength is None:
        maxReadLength = genomeLen
    elif maxReadLength <= 0:
        raise ValueError('The maximum read length must be > 0')

    alignReads = args.alignReads
    circularGenome = args.circularGenome

    if circularGenome and alignReads:
        raise ValueError('You cannot specify both --circularGenome and '
                         '--alignReads')

    idPrefix = args.idPrefix
    verbose = args.verbose
    genomeSequence = genome.sequence
    readCountWidth = int(log10(args.count)) + 1
    genomeLengthWidth = int(log10(genomeLen)) + 1

    if args.printGenome:
        print(genome.toString('fasta'), end='')

    for i in range(args.count):
        id_ = '%s%0*d' % (idPrefix, readCountWidth, i + 1)
        read, offset, mutationOffsets = makeRead(
            genomeSequence, genomeLen, meanLength, sdLength,
            minReadLength, maxReadLength, id_, rate, circularGenome)

        read.id = read.id + '-length-%0*d-offset-%0*d' % (
            genomeLengthWidth, len(read),
            genomeLengthWidth, offset)

        if mutationOffsets:
            read.id = read.id + '-mutations-at-%s' % (
                ','.join(map(str, sorted(mutationOffsets))))
        else:
            read.id = read.id + '-no-mutations'

        if verbose:
            print('Created read of length %d with %d mutations' %
                  (len(read), len(mutationOffsets)), file=sys.stderr)

        if alignReads:
            sequence = ('-' * offset) + read.sequence
            if len(sequence) < genomeLen:
                sequence += '-' * (genomeLen - len(sequence))
            read.sequence = sequence[:genomeLen]

        print(read.toString('fasta'), end='')
