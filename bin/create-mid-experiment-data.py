#!/usr/bin/env python

import argparse
import sys
from random import choice
from pathlib import Path

from dark.fasta import FastaReads
from dark.process import Executor
from dark.reads import Reads

from midtools.mutate import mutateRead


def main(args, outDir, logfp):
    """
    Create genomes and reads for a multiple infection detection experiment.

    @param args: An argparse C{Namespace} instance, as returned by C{parse_args}.
    @param outDir: A C{Path} to the output directly.
    @param logfp: A file object to write log information to.
    """
    print("Invocation arguments", args, file=logfp)

    # Use genome1 values if the genome2 values are not specified.
    if args.genome2ReadCount is None:
        args.genome2ReadCount = args.genome1ReadCount
        print(
            f"Genome 2 read count defaulting to genome 1 value "
            f"({args.genome2ReadCount}).",
            file=sys.stderr,
        )

    if args.genome2ReadMutationRate is None:
        args.genome2ReadMutationRate = args.genome1ReadMutationRate
        print(
            f"Genome 2 read mutation rate defaulting to genome 1 value "
            f"({args.genome2ReadMutationRate:.3f}).",
            file=sys.stderr,
        )

    genome1 = outDir / "genome-1.fasta"
    genome2 = outDir / "genome-2.fasta"
    genome2locations = outDir / "genome-2.locations"

    reads1 = outDir / "reads-1.fastq"
    reads2 = outDir / "reads-2.fastq"
    reads12 = outDir / "reads-12.fastq"

    executor = Executor()
    execute = executor.execute

    # Make genome 1.
    if args.genome1Filename:
        execute(f"ln -s {args.genome1Filename!r} {str(genome1)!r}")
    else:
        if args.genomeLength < 1:
            print("Random initial genome length must be > 0.", file=sys.stderr)
            sys.exit(1)
        print(
            f"Writing random starting genome of length {args.genomeLength} to "
            f"{str(genome1)!r}",
            file=logfp,
        )
        sequence = "".join([choice("ACGT") for _ in range(args.genomeLength)])
        with open(genome1, "w") as fp:
            print(">genome-1\n%s" % sequence, file=fp)

    # Make genome 2.
    if args.genome2Filename:
        execute(f"ln -s {args.genome2Filename!r} {str(genome2)!r}")
    else:
        # Make a second genome (from the first genome) using the given
        # mutation rate. Print its mutated locations to a file.
        (genome1read,) = list(FastaReads(genome1))
        offsets = mutateRead(genome1read, args.genome2MutationRate)
        with open(genome2locations, "w") as fp:
            print("\n".join(str(offset + 1) for offset in sorted(offsets)), file=fp)
        genome1read.id = "genome-2"
        Reads([genome1read]).save(genome2)

    # Make reads.
    cmdPrefix = (
        f"create-reads.py "
        f"--maxReadLength {args.maxReadLength} "
        f"--minReadLength {args.minReadLength} "
        f"--meanLength {args.meanReadLength} "
        f"--sdLength {args.sdReadLength} "
    )

    for info in [
        {
            "reads": str(reads1),
            "fasta": str(genome1),
            "number": 1,
            "count": args.genome1ReadCount,
            "rate": args.genome1ReadMutationRate,
        },
        {
            "reads": str(reads2),
            "fasta": str(genome2),
            "number": 2,
            "count": args.genome2ReadCount,
            "rate": args.genome2ReadMutationRate,
        },
    ]:
        execute(
            cmdPrefix
            + (
                "--idPrefix genome-%(number)d-read- "
                "--count %(count)d < %(fasta)s > %(reads)s" % info
            )
        )

    execute(f"cat {str(reads1)!r} {str(reads2)!r} > {str(reads12)!r}")

    print("\n".join(executor.log), file=logfp)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Create data for a double-infection experiment.",
    )

    parser.add_argument("--outDir", required=True, help="The output directory.")

    parser.add_argument(
        "--genome1Filename",
        help=(
            "The FASTA file containing the first genome to create reads "
            "from. If not specified, a random genome will be created."
        ),
    )

    parser.add_argument(
        "--genome2Filename",
        help=(
            "The FASTA file containing the second genome to create reads "
            "from. If not specified, a genome will be created from the "
            "first genome using a mutation rate given by "
            "--genome2MutationRate."
        ),
    )

    parser.add_argument(
        "--genome2MutationRate",
        type=float,
        default=0.05,
        help=(
            "The per-base mutation rate to use to create the second genome from "
            "genome 1 (if that genome is not provided by --genome2File)."
        ),
    )

    parser.add_argument(
        "--genome1ReadMutationRate",
        type=float,
        help="The per-base mutation rate to use when creating reads from genome 1.",
    )

    parser.add_argument(
        "--genome2ReadMutationRate",
        type=float,
        help=(
            "The per-base mutation rate to use when creating reads from the second "
            "genome. If not given, the --genome1ReadMutationRate value will be used."
        ),
    )

    parser.add_argument(
        "--genome1ReadCount",
        default=100,
        type=int,
        help="The number of reads to create for genome 1.",
    )

    parser.add_argument(
        "--genome2ReadCount",
        type=int,
        help=(
            "The number of reads to create for genome 2. If not "
            "given, the value of --genome1ReadCount will be used"
        ),
    )

    parser.add_argument(
        "--genomeLength",
        type=int,
        default=100,
        help="If any random genomes need to be made, this will be their length.",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="If specified, overwrite the contents of the --outDir directory",
    )

    parser.add_argument(
        "--minReadLength",
        type=int,
        default=30,
        help="The minimum length read to create.",
    )

    parser.add_argument(
        "--maxReadLength",
        type=int,
        default=100,
        help="The maximum length read to create.",
    )
    parser.add_argument(
        "--meanReadLength",
        type=float,
        default=80.0,
        help="The mean length of created reads.",
    )

    parser.add_argument(
        "--sdReadLength",
        type=float,
        default=10.0,
        help="The standard deviation of the length of created reads.",
    )

    parser.add_argument(
        "--qualityChar",
        default="I",
        help="The quality character to use for all read quality scores.",
    )

    args = parser.parse_args()

    outDir = Path(args.outDir)

    if outDir.exists():
        if not args.force:
            print(
                f"Output directory {str(outDir)!r} already exists. Exiting.",
                file=sys.stderr,
            )
            sys.exit(1)
    else:
        outDir.mkdir()

    with open(outDir / "LOG", "w") as logfp:
        main(args, outDir, logfp)
