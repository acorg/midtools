#!/usr/bin/env python

from __future__ import division, print_function

from itertools import chain

from midtools.clusterAnalysis import ClusterAnalysis
from midtools.options import addAnalysisCommandLineOptions


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Find which reads agree and disagree with one another.')

    addAnalysisCommandLineOptions(parser)

    parser.add_argument(
        '--cutoff', type=float, default=ClusterAnalysis.DEFAULT_CUTOFF,
        help=('Clustering will be stopped once the minimum distance between '
              'remaining clusters exceeds this value.'))

    args = parser.parse_args()

    referenceIds = (list(chain.from_iterable(args.referenceId))
                    if args.referenceId else None)
    ClusterAnalysis(
        alignmentFiles=list(chain.from_iterable(args.alignmentFile)),
        referenceGenomeFiles=list(chain.from_iterable(args.referenceGenome)),
        referenceIds=referenceIds,
        cutoff=args.cutoff,
        outputDir=args.outputDir,
        minReads=args.minReads,
        homogeneousCutoff=args.homogeneousCutoff,
        saveReducedFASTA=args.saveReducedFASTA,
        plotSAM=args.plotSAM,
        verbose=args.verbose).run()
