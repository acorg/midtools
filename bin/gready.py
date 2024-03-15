#!/usr/bin/env python

from itertools import chain

from midtools.greadyAnalysis import GreadyAnalysis
from midtools.options import addAnalysisCommandLineOptions


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Find which reads agree and disagree with one another.",
    )

    addAnalysisCommandLineOptions(parser)

    parser.add_argument(
        "--cutoff",
        type=float,
        default=0.5,
        help=("Reads with a score less than this will not be put into the consensus."),
    )

    args = parser.parse_args()

    referenceIds = (
        list(chain.from_iterable(args.referenceId)) if args.referenceId else None
    )
    GreadyAnalysis(
        alignmentFiles=list(chain.from_iterable(args.alignmentFile)),
        referenceGenomeFiles=list(chain.from_iterable(args.referenceGenome)),
        referenceIds=referenceIds,
        cutoff=args.cutoff,
        outputDir=args.outputDir,
        minReads=args.minReads,
        homogeneousCutoff=args.homogeneousCutoff,
        saveReducedFASTA=args.saveReducedFASTA,
        plotSAM=args.plotSAM,
        verbose=args.verbose,
    ).run()
