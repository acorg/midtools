#!/usr/bin/env python

import sys

from midtools.options import addCommandLineOptions, parseCommandLineOptions
from midtools.plotting import plotBaseFrequencies


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Show significant genome location nucleotide "
            "frequencies for a set of aligned reads."
        ),
    )

    addCommandLineOptions(parser, "significant-base-frequencies.html")

    parser.add_argument(
        "--sampleName",
        help=("The name of the sample, to appear in the --valuesFile (if given)."),
    )

    parser.add_argument("--title", help="The plot title.")

    parser.add_argument(
        "--sortOn",
        choices=("max", "entropy"),
        help=(
            "If specified, locations will be sorted according to either the "
            "maximum nucleotide frequency or the nucleotide entropy at the "
            "location."
        ),
    )

    parser.add_argument(
        "--histogram",
        action="store_true",
        default=False,
        help=(
            "If specified and --sortOn is used, the values (according to "
            "--sortOn) will be shown in a histogram."
        ),
    )

    parser.add_argument(
        "--valuesFile",
        help=(
            "The filename to write the raw max or entropy scores to (as a "
            "JSON list). This option can only be used if --sortOn is used."
        ),
    )

    parser.add_argument(
        "--titleFontSize",
        type=int,
        default=12,
        help="The font size for the plot title.",
    )

    parser.add_argument(
        "--axisFontSize",
        type=int,
        default=12,
        help="The font size for the axis titles.",
    )

    args = parser.parse_args()

    if args.valuesFile and args.sortOn is None:
        print(
            "You have specified a file to write sorted location values "
            "to (using --valuesFile) but have not specified what to "
            'sort on using --sortOn. Choices are "max" or "entropy".',
            file=sys.stderr,
        )
        sys.exit(1)

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
        "Read %d aligned reads of length %d. "
        "Found %d significant locations."
        % (len(alignedReads), genomeLength, len(significantOffsets))
    )

    plotBaseFrequencies(
        significantOffsets,
        baseCountAtOffset,
        readCountAtOffset,
        args.outfile,
        valuesFile=args.valuesFile,
        title=args.title,
        sampleName=args.sampleName,
        minReads=args.minReads,
        homogeneousCutoff=args.homogeneousCutoff,
        sortOn=args.sortOn,
        histogram=args.histogram,
        show=args.show,
        titleFontSize=args.titleFontSize,
        axisFontSize=args.axisFontSize,
    )
