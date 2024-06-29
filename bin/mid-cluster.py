#!/usr/bin/env python

from itertools import chain

from midtools.analysis import ReadAnalysis
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

    parser.add_argument(
        "--minCCIdentity",
        type=float,
        default=ClusterAnalysis.MIN_CC_IDENTITY_DEFAULT,
        help=(
            "The minimum nucleotide identity fraction [0.0, 1.0] that a consistent "
            "component must have with a reference in order to contribute to the "
            "consensus being made against the reference."
        ),
    )

    parser.add_argument(
        "--noCoverageStrategy",
        default="N",
        choices=("N", "reference"),
        help=(
            "The approach to use when making a consensus if there are no reads "
            "covering a site. A value of 'N' means to use an ambigous N nucleotide "
            "code, whereas a value fo 'reference' means to take the base from the "
            "reference sequence."
        ),
    )

    args = parser.parse_args()

    referenceIds = (
        list(chain.from_iterable(args.referenceId)) if args.referenceId else None
    )

    analysis = ReadAnalysis(
        args.sampleName,
        list(chain.from_iterable(args.alignmentFile)),
        list(chain.from_iterable(args.referenceGenome)),
        args.outputDir,
        referenceIds=referenceIds,
        minReads=args.minReads,
        homogeneousCutoff=args.homogeneousCutoff,
        plotSAM=args.plotSAM,
        saveReducedFASTA=args.saveReducedFASTA,
        verbose=args.verbose,
    )

    clusterAnalysis = ClusterAnalysis(
        analysis,
        maxClusterDist=args.maxClusterDist,
        alternateNucleotideMinFreq=args.alternateNucleotideMinFreq,
        minCCIdentity=args.minCCIdentity,
        noCoverageStrategy=args.noCoverageStrategy,
    )

    analysis.run(clusterAnalysis.analyzeReference)
