#!/usr/bin/env python

from itertools import chain
from midtools.connectedComponentAnalysis import ConnectedComponentAnalysis
from midtools.options import addAnalysisCommandLineOptions


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Find which reads agree and disagree with one another.",
    )

    addAnalysisCommandLineOptions(parser)

    parser.add_argument(
        "--agreementThreshold",
        type=float,
        default=ConnectedComponentAnalysis.DEFAULT_AGREEMENT_THRESHOLD,
        help=(
            "Only reads with agreeing nucleotides at at least this fraction "
            "of the significant sites they have in common will be "
            "considered connected (this is for the second phase of adding "
            "reads to a component."
        ),
    )

    args = parser.parse_args()

    referenceIds = (
        list(chain.from_iterable(args.referenceId)) if args.referenceId else None
    )
    ConnectedComponentAnalysis(
        alignmentFiles=list(chain.from_iterable(args.alignmentFile)),
        referenceGenomeFiles=list(chain.from_iterable(args.referenceGenome)),
        referenceIds=referenceIds,
        outputDir=args.outputDir,
        minReads=args.minReads,
        homogeneousCutoff=args.homogeneousCutoff,
        agreementThreshold=args.agreementThreshold,
        saveReducedFASTA=args.saveReducedFASTA,
        plotSAM=args.plotSAM,
        verbose=args.verbose,
    ).run()
