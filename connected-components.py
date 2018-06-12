#!/usr/bin/env python

from __future__ import division, print_function

from itertools import chain
from data.analysis import ReadAnalysis


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Find which reads agree and disagree with one another.')

    parser.add_argument(
        '--alignmentFile', metavar='FILENAME', action='append', nargs='+',
        help='The name of a SAM/BAM alignment file (may be repeated).')

    parser.add_argument(
        '--referenceGenome', metavar='FILENAME', action='append', nargs='+',
        help=('The name of a FASTA file containing reference genomes that '
              'were used to create the alignment files (may be repeated).'))

    parser.add_argument(
        '--referenceId', metavar='NAME', action='append', nargs='*',
        help=('The sequence id whose alignment should be analyzed (may '
              'be repeated). All ids must be present in --referenceGenome '
              'file. One of the SAM/BAM files given using --alignmentFile '
              'should have an alignment against the given argument. If '
              'omitted, all references that are aligned to in the given '
              'BAM/SAM files will be analyzed.'))

    parser.add_argument(
        '--outputDir',
        help='The directory to save result files to.')

    parser.add_argument(
        '--minReads', type=int, default=ReadAnalysis.DEFAULT_MIN_READS,
        help=('The minimum number of reads that must cover a location for it '
              'to be considered significant.'))

    parser.add_argument(
        '--homogeneousCutoff', type=float,
        default=ReadAnalysis.DEFAULT_HOMOGENEOUS_CUTOFF,
        help=('If the most common nucleotide at a location occurs more than '
              'this fraction of the time (i.e., amongst all reads that cover '
              'the location) then the locaion will be considered homogeneous '
              'and therefore uninteresting.'))

    parser.add_argument(
        '--agreementThreshold', type=float,
        default=ReadAnalysis.DEFAULT_AGREEMENT_THRESHOLD,
        help=('Only reads with agreeing nucleotides at at least this fraction '
              'of the significant sites they have in common will be '
              'considered connected (this is for the second phase of adding '
              'reads to a component.'))

    parser.add_argument(
        '--saveReducedFASTA', default=False, action='store_true',
        help=('If given, write out a FASTA file of the original input but '
              'with just the signifcant locations.'))

    parser.add_argument(
        '--verbose', action='store_true', default=False,
        help='Print verbose output showing read connections.')

    args = parser.parse_args()

    ReadAnalysis(
        alignmentFiles=list(chain.from_iterable(args.alignmentFile)),
        referenceGenomeFiles=list(chain.from_iterable(args.referenceGenome)),
        referenceIds=list(chain.from_iterable(args.referenceId)),
        outputDir=args.outputDir,
        minReads=args.minReads,
        homogeneousCutoff=args.homogeneousCutoff,
        agreementThreshold=args.agreementThreshold,
        saveReducedFASTA=args.saveReducedFASTA,
        verbose=args.verbose).run()
