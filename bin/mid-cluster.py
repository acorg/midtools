#!/usr/bin/env python

from itertools import chain

from midtools.clusterAnalysis import ClusterAnalysis
from midtools.options import addAnalysisCommandLineOptions


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Make consensus sequences by aligning reads to references "
            "and finding which reads agree and disagree with one another "
            "using clustering."
        ),
    )

    addAnalysisCommandLineOptions(parser)

    parser.add_argument(
        "--maxClusterDist",
        type=float,
        default=ClusterAnalysis.DEFAULT_MAX_CLUSTER_DIST,
        help=(
            "Clustering will be stopped once the minimum distance between "
            "remaining clusters exceeds this value."
        ),
    )

    parser.add_argument(
        "--alternateNucleotideMinFreq",
        type=float,
        default=ClusterAnalysis.ALTERNATE_NUCLEOTIDE_MIN_FREQ_DEF,
        help=(
            "The (0.0 to 1.0) frequency that an alternative nucleotide "
            "(i.e., not the one chosen for a consensus) must have in order "
            "to be selected for the alternate consensus."
        ),
    )

    args = parser.parse_args()

    referenceIds = (
        list(chain.from_iterable(args.referenceId)) if args.referenceId else None
    )

    ClusterAnalysis(
        args.sampleName,
        list(chain.from_iterable(args.alignmentFile)),
        list(chain.from_iterable(args.referenceGenome)),
        args.outputDir,
        referenceIds=referenceIds,
        minReads=args.minReads,
        maxClusterDist=args.maxClusterDist,
        homogeneousCutoff=args.homogeneousCutoff,
        plotSAM=args.plotSAM,
        alternateNucleotideMinFreq=args.alternateNucleotideMinFreq,
        saveReducedFASTA=args.saveReducedFASTA,
        verbose=args.verbose,
    ).run()
