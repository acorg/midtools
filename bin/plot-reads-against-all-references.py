#!/usr/bin/env python

import argparse

from dark.sam import SAMFilter
from midtools.plotting import plotAllReferencesSAM


def makeParser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Make a plot aligning all reads to all references in a SAM/BAM file."
        ),
    )

    parser.add_argument(
        "--sampleName", default="sample", help="The name of the sample."
    )
    parser.add_argument("--samFile", required=True, help="The SAM/BAM file to read.")
    parser.add_argument(
        "--outputFile", required=True, help="The output (HTML) file for the plot."
    )
    parser.add_argument(
        "--hbv",
        action="store_true",
        help="Make an HBV-specific plot, categorized by genotype.",
    )

    return parser


def main() -> None:
    args = makeParser().parse_args()

    plotAllReferencesSAM(
        SAMFilter(args.samFile),
        args.outputFile,
        hbv=args.hbv,
        sampleName=args.sampleName,
        alignmentFile=args.samFile,
        jitter=0.5,
    )


if __name__ == "__main__":
    main()
