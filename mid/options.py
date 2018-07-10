from dark.sam import PaddedSAM

from mid.data import gatherData, findSignificantOffsets
from mid.read import AlignedRead


def addCommandLineOptions(parser, outfileDefaultName=None):
    """
    Add standard command-line options to an argument parser.

    @param parser: An C{ArgumentParser} instance.
    @param outfileDefaultName: The C{str} output file to use as a default
        in case the user does not give one on the command line.
    """
    parser.add_argument(
        '--samFile', metavar='FILENAME', required=True,
        help='The name of the SAM/BAM input file.')

    parser.add_argument(
        '--referenceId', metavar='SEQUENCE-ID',
        help='The id of the reference to use in the SAM/BAM input file.')

    parser.add_argument(
        '--outFile', default=outfileDefaultName,
        help='The filename to store the resulting HTML.')

    parser.add_argument(
        '--minReads', type=int, default=5,
        help=('The minimum number of reads that must cover a location for it '
              'to be considered significant.'))

    parser.add_argument(
        '--homogeneousCutoff', type=float, default=0.9,
        help=('If the most common nucleotide at a location occurs more than '
              'this fraction of the time (i.e., amongst all reads that cover '
              'the location) then the locaion will be considered homogeneous '
              'and therefore uninteresting.'))

    parser.add_argument(
        '--show', action='store_true', default=False,
        help='If specified, show the figure interactively.')


def parseCommandLineOptions(args, returnSignificantOffsets=True):
    """
    Deal with the various command-line options added to the ArgumentParser
    instance by addCommandLineOptions.

    @param args: The result of calling C{parse_args} on an C{ArgumentParser}
        instance (the one that was passed to C{addCommandLineOptions}, unless
        we're testing).
    @param returnSignificantOffsets: If C{True} also return a list of the
        significant offsets (else that element of the return value will be
        C{None}).
    @return: A C{tuple}: (genomeLength, alignedReads, padddedSAM,
        readCountAtOffset, baseCountAtOffset, readsAtOffset,
        significantOffsets).
    """
    genomeLength = None
    alignedReads = []
    paddedSAM = PaddedSAM(args.samFile)
    sam = paddedSAM.samfile

    referenceId = args.referenceId

    if referenceId:
        # Make sure the reference id is present in the SAM file.
        tid = sam.get_tid(referenceId)
        if tid == -1:
            raise ValueError('Reference %s not in %s alignment file.' %
                             (referenceId, args.samFile))
        else:
            genomeLength = sam.lengths[tid]
    else:
        # If there is just one reference, use it. Otherwise exit.
        if sam.nreferences == 1:
            referenceId = sam.references[0]
        else:
            # The error message will finish with ' ().' if for some reason
            # the SAM file has zero references. But that probably cannot
            # happen so I'm not worrying about it.
            raise ValueError(
                'If you do not specify a reference sequence with '
                '--referenceId, the SAM/BAM file must contain exactly one '
                'reference. But %s contains %d (%s).' %
                (args.samFile, sam.nreferences,
                 ', '.join(sorted(sam.references))))

    for query in paddedSAM.queries(referenceName=referenceId):
        if genomeLength is None:
            genomeLength = len(query)
        alignedReads.append(AlignedRead(query.id, query.sequence))

    readCountAtOffset, baseCountAtOffset, readsAtOffset = gatherData(
        genomeLength, alignedReads)

    if returnSignificantOffsets:
        significantOffsets = list(findSignificantOffsets(
            baseCountAtOffset, readCountAtOffset, args.minReads,
            args.homogeneousCutoff))
        for read in alignedReads:
            read.setSignificantOffsets(significantOffsets)
    else:
        significantOffsets = None

    return (genomeLength, alignedReads, paddedSAM, readCountAtOffset,
            baseCountAtOffset, readsAtOffset, significantOffsets)
