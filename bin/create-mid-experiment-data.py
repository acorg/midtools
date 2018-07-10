#!/usr/bin/env python

from __future__ import division, print_function

import argparse
import sys
from os import mkdir
from os.path import exists, join
from random import choice
from six.moves import shlex_quote as quote

from dark.process import Executor


def main(args, logfp):
    """
    Create genomes and reads for a multiple infection detection experiment.

    @param args: A namespace instance, as returned by parse_args
    @param logfp: A file object to write log information to.
    """
    print('Invocation arguments', args, file=logfp)
    qOutputDir = quote(args.outputDir)
    genome1 = join(qOutputDir, 'genome-1.fasta')
    genome2 = join(qOutputDir, 'genome-2.fasta')
    reads1 = join(qOutputDir, 'reads-1.fasta')
    reads2 = join(qOutputDir, 'reads-2.fasta')
    reads12 = join(qOutputDir, 'reads-12.fasta')

    executor = Executor(args.dryRun)

    if args.genome1Filename:
        executor.execute('ln -s %s %s' %
                         (quote(args.genome1Filename), genome1))
    else:
        if args.genomeLength < 1:
            print('Random initial genome length must be > 0.', file=sys.stderr)
            sys.exit(3)
        print('Writing random starting genome of length %d to %s' %
              (args.genomeLength, genome1), file=logfp)
        if not args.dryRun:
            with open(genome1, 'w') as fp:
                print('>genome-1\n%s' % (
                    ''.join([choice('ACGT')
                             for _ in range(args.genomeLength)])), file=fp)

    if args.genome2Filename:
        executor.execute('ln -s %s %s' %
                         (quote(args.genome2Filename), genome2))
    else:
        # Make a second genome using the given mutation rate.
        executor.execute(
            'mutate-reads.py --rate %s < %s | '
            'filter-fasta.py --quiet '
            '--idLambda \'lambda _: "genome-2"\'> %s' %
            (args.genome2MutationRate, genome1, genome2))

    cmdPrefix = (
        'create-reads.py --maxReadLength %d --minReadLength %d '
        '--meanLength %d --sdLength %d --rate %f ' %
        (args.maxReadLength, args.minReadLength, args.meanReadLength,
         args.sdReadLength, args.readMutationRate))

    for info in [
            {
                'reads': reads1,
                'fasta': genome1,
                'number': 1,
                'count': args.genome1ReadCount or args.readCount,
            },
            {
                'reads': reads2,
                'fasta': genome2,
                'number': 2,
                'count': args.genome2ReadCount or args.readCount,
            }]:
        executor.execute(cmdPrefix + ('--idPrefix genome-%(number)d-read- '
                         '--count %(count)d < %(fasta)s > %(reads)s' % info))

    executor.execute('cat %s %s > %s' % (reads1, reads2, reads12))

    print('\n'.join(executor.log), file=logfp)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Create data for a dual infection experiment')

    parser.add_argument(
        '--outputDir', required=True,
        help='The output directory. Must not already exist.')

    parser.add_argument(
        '--genome1Filename',
        help=('The FASTA file containing the first genome to create reads '
              'from. If not specified, a random genome will be created.'))

    parser.add_argument(
        '--genome2Filename',
        help=('The FASTA file containing the second genome to create reads '
              'from. If not specified, a genome will be created from the '
              'first genome.'))

    parser.add_argument(
        '--genomeLength', type=int, default=100,
        help=('If any random genomes need to be made, this will be their '
              'length.'))

    parser.add_argument(
        '--dryRun', action='store_true', default=False,
        help='If specified, simply print the actions that would be taken.')

    parser.add_argument(
        '--force', action='store_true', default=False,
        help=('If specified, overwrite the contents of the --outputDir '
              'directory'))

    parser.add_argument(
        '--genome2MutationRate', type=float, default=0.075,
        help='The per-base mutation rate to use to create the second genome '
        '(if it is not provided by --genome2File).')

    parser.add_argument(
        '--readCount', type=int, default=100,
        help=('The number of reads to create for both genomes (only used '
              'if --genome1ReadCount or --genome2ReadCount are not given).'))

    parser.add_argument(
        '--genome1ReadCount', type=int,
        help='The number of reads to create for genome 1.')

    parser.add_argument(
        '--genome2ReadCount', type=int,
        help='The number of reads to create for genome 2.')

    parser.add_argument(
        '--minReadLength', type=int, default=10,
        help='The minimum length read to create.')

    parser.add_argument(
        '--maxReadLength', type=int, default=100,
        help='The maximum length read to create.')

    parser.add_argument(
        '--readMutationRate', type=float, default=0.01,
        help='The per-base mutation rate to use to create reads.')

    parser.add_argument(
        '--meanReadLength', type=float, default=100.0,
        help='The mean length of created reads.')

    parser.add_argument(
        '--sdReadLength', type=float, default=10.0,
        help='The standard deviation of the length of created reads.')

    args = parser.parse_args()

    if exists(args.outputDir):
        if not args.dryRun and not args.force:
            print('Output directory %r already exists. Exiting.' %
                  args.outputDir, file=sys.stderr)
            sys.exit(2)
    else:
        if not args.dryRun:
            mkdir(args.outputDir)

    with open(join(args.outputDir, 'LOG'), 'w') as logfp:
        main(args, logfp)
