#!/usr/bin/env python

from __future__ import division, print_function

from midtools.options import addCommandLineOptions, parseCommandLineOptions
from midtools.plotting import plotCoverageAndSignificantLocations


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Plot read coverage and significant locations",
    )

    parser.add_argument("--title", help="The overall title plot")

    addCommandLineOptions(parser, "coverage-and-significant-locations.html")
    args = parser.parse_args()

    (
        genomeLength,
        alignedReads,
        paddedSAM,
        readCountAtOffset,
        baseCountAtOffset,
        readsAtOffset,
        significantOffsets,
    ) = parseCommandLineOptions(args)

    print(
        "Read %d aligned reads. Found %d significant locations."
        % (len(alignedReads), len(significantOffsets))
    )

    plotCoverageAndSignificantLocations(
        readCountAtOffset,
        genomeLength,
        significantOffsets,
        args.outfile,
        title=args.title,
        show=args.show,
    )
