#!/usr/bin/env python

import sys
from collections import Counter
from math import log10

from midtools.options import addCommandLineOptions, parseCommandLineOptions
from midtools.utils import baseCountsToStr


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Analyze a set of aligned reads.",
    )

    addCommandLineOptions(parser)

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print verbose textual output showing read connections.",
    )

    args = parser.parse_args()

    (
        genomeLength,
        alignedReads,
        paddedSAM,
        readCountAtOffset,
        baseCountAtOffset,
        readsAtOffset,
        _,
    ) = parseCommandLineOptions(args, False)

    print("Read %d aligned reads." % len(alignedReads), file=sys.stderr)

    genomeLengthWidth = int(log10(genomeLength)) + 1
    nucleotides = set("ACGT")

    for offset in range(genomeLength):
        counts = Counter()
        for read in readsAtOffset[offset]:
            base = read.base(offset)
            if base in nucleotides:
                counts[base] += 1
        print(
            "Location %*d: base counts %s"
            % (genomeLengthWidth, offset + 1, baseCountsToStr(counts))
        )
