#!/usr/bin/env python

from __future__ import division, print_function

import sys
import plotly
import plotly.graph_objs as go
import re
import colorlover
from math import log10
from tempfile import mkdtemp
from os.path import exists
from os import mkdir

from dark.fasta import FastaReads

from data.analysis import ReadAnalysis
from data.data import parseCommandLineOptions
from data.utils import baseCountsToStr, nucleotidesToStr


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Find which reads agree and disagree with one another.')

    parser.add_argument(
        '--verbose', action='store_true', default=False,
        help='Print verbose textual output showing read connections.')

    parser.add_argument(
        '--saveReducedFASTA', default=False, action='store_true',
        help=('If given write out a FASTA file of the original input but '
              'with just the signifcant locations'))

    parser.add_argument(
        '--outputDir',
        help=('The directory to write result files to.'))

    parser.add_argument(
        '--agreementThreshold', type=float, default=0.9,
        help=('Only reads with agreeing nucleotides at at least this fraction '
              'of the significant sites they have in common will be '
              'considered connected (this is for the second phase of adding '
              'reads to a component.'))

    parser.add_argument(
        '--fastaFile', metavar='FILENAME', required=True,
        help='The name of the FASTA input file.')

    parser.add_argument(
        '--referenceGenome', metavar='FILENAME', action='append',
        help=('The name of a FASTA file containing reference genomes (may '
              'be repeated).'))

    parser.add_argument(
        '--minReads', type=int, default=5,
        help=('The minimum number of reads that must cover a location for it '
              'to be considered significant.'))

    parser.add_argument(
        '--homogeneousCutoff', type=float, default=0.9,
        help=('If the most common nucleotide at a location occurs more than '
              'this fraction of the time (i.e., amongst all reads that cover '
              'the location) then the locaion will be considered homogeneous '
              'and therefore uninteresting.'))

    parser.add_argument(
        '--trim', type=int, default=0,
        help=('The number of bases to trim from the start and end of each '
              'read (to try to remove possible damage).'))

    args = parser.parse_args()
    verbose = args.verbose

    (genomeLength, alignedReads, readCountAtOffset,
     baseCountAtOffset, readsAtOffset,
     significantOffsets) = parseCommandLineOptions(args, True)

    if verbose:
        print('Read %d aligned reads of length %d. '
              'Found %d significant locations.' %
              (len(alignedReads), genomeLength, len(significantOffsets)))

    if not significantOffsets:
        print('Exiting due to finding no significant locations.')
        sys.exit(2)

    if args.outputDir:
        outputDir = args.outputDir
        if not exists(outputDir):
            mkdir(outputDir)
    else:
        outputDir = mkdtemp()
        print('Writing output files to %s' % outputDir)

    ra = ReadAnalysis(
        args.fastaFile, outputDir, genomeLength, alignedReads,
        readCountAtOffset, baseCountAtOffset, readsAtOffset,
        significantOffsets, args.agreementThreshold, verbose,
        args.referenceGenome)

    if verbose:
        print('Saving significant offsets')
    ra.saveSignificantOffsets()

    if verbose:
        print('Saving component FASTA')
    ra.saveComponentFasta()

    if args.saveReducedFASTA:
        if verbose:
            print('Saving reduced FASTA')
        ra.saveReducedFasta()

    if verbose:
        print('Writing analysis summary')
    ra.summarize()

    if verbose:
        print('Saving component consensuses')
    ra.saveComponentConsensuses()

    if verbose:
        print('Saving closest consensuses to references')
    ra.saveClosestReferenceConsensuses()
