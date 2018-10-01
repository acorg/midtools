#!/usr/bin/env python

from __future__ import print_function

import sys
from random import uniform, normalvariate
from math import log10

from dark.reads import (
    Read, addFASTACommandLineOptions, parseFASTACommandLineOptions)

from midtools.mutate import mutateRead
from midtools.utils import s


def makeRead(genome, meanLength, sdLength, minReadLength, maxReadLength,
             id_, rate, circularGenome):
    """
    Make a read, according to various parameters and constraints regarding its
    length.

    Note that when circularGenome is False, reads generated using this method
    will not in fact have a mean length of C{meanLength}. This is because they
    are sometimes truncated at the start and end of the genome.

    @param genome: The C{str} genome to base the read on.
    @param meanLength: The C{float} mean read length.
    @param sdLength: The C{float} standard deviation of the read lengths.
    @param minReadLength: The C{int} minimum read length.
    @param maxReadLength: The C{int} maximum read length.
    @param id_: The C{str} read id.
    @param rate: The per-base C{float} mutation rate.
    @param circularGenome: If C{True}, the genome will be treated as circular.
        Reads that would otherwise be truncated by running into the end of the
        genome will continue with bases from the start of the genome.
    """
    genomeLen = len(genome)
    length = -1

    while (0 >= length > genomeLen or
           length < minReadLength or
           length > maxReadLength):
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
        # For symmetry, we calculate an offset that allows the read to
        # overlap (by at least minReadLength bases) with the start or end
        # of the genome. If that happens, we truncate the read.
        offset = int(uniform(-(length - 1) + minReadLength,
                             genomeLen - minReadLength))

        if offset < 0:
            sequence = genome[:offset + length]
        else:
            sequence = genome[offset:offset + length]

    assert maxReadLength >= len(sequence) >= minReadLength, (
        'maxReadLength=%d, len(sequence)=%d, minReadLength=%d '
        'readLength=%d offset=%d' %
        (maxReadLength, len(sequence), minReadLength, length, offset))

    read = Read(id_, sequence)
    mutationOffsets = () if rate == 0.0 else mutateRead(read, rate)
    return read, offset, mutationOffsets


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Create DNA reads.')

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
        '--fastaReads', action='store_true', default=False,
        help='Make the reads be FASTA instead of FASTQ')

    parser.add_argument(
        '--qualityChar', default='I',
        help=('The quality character to use for all quality scores when '
              '--fastq is used'))

    parser.add_argument(
        '--circularGenome', action='store_true', default=False,
        help=('If specified, reads will wrap around the genome (currently not '
              'compatible with --alignReads).'))

    parser.add_argument(
        '--printGenome', action='store_true', default=False,
        help='If specified, print the genome as the first sequence.')

    parser.add_argument(
        '--alignReads', action='store_true', default=False,
        help=('If specified, print the reads aligned (with "-" characters) '
              'to the genome.'))

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = list(parseFASTACommandLineOptions(args))
    # There should only be one "read", the sequence we are to create other
    # reads from.
    assert len(reads) == 1, (
        'FASTA input contained %d sequence%s (expected just one).' % (
            len(reads), s(len(reads))))
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
        raise ValueError('The minimum read length must be positive')

    maxReadLength = args.maxReadLength

    if maxReadLength is None:
        maxReadLength = genomeLen
    elif maxReadLength <= 0:
        raise ValueError('The maximum read length must be positive')

    if minReadLength > maxReadLength:
        raise ValueError(
            'The minimum read length cannot exceed the maximum read length')

    alignReads = args.alignReads
    circularGenome = args.circularGenome

    if circularGenome and alignReads:
        raise ValueError(
            'You cannot specify both --circularGenome and --alignReads')

    idPrefix = args.idPrefix
    verbose = args.verbose
    genomeSequence = genome.sequence
    readCountWidth = int(log10(args.count)) + 1
    genomeLengthWidth = int(log10(genomeLen)) + 1

    if args.printGenome:
        print(genome.toString('fasta'), end='')

    fastq, format_ = (False, 'fasta') if args.fastaReads else (True, 'fastq')
    qualityChar = args.qualityChar

    for i in range(args.count):
        id_ = '%s%0*d' % (idPrefix, readCountWidth, i + 1)
        read, offset, mutationOffsets = makeRead(
            genomeSequence, meanLength, sdLength,
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

        if fastq:
            read.quality = qualityChar * len(read.sequence)

        print(read.toString(format_), end='')
