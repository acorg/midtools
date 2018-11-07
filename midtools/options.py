from dark.sam import SAMFilter, PaddedSAM

from midtools.data import gatherData, findSignificantOffsets
from midtools.read import AlignedRead


def addCommonOptions(parser):
    """
    Add standard command-line options to an argument parser.

    @param parser: An C{ArgumentParser} instance.
    """
    parser.add_argument(
        '--minReads', type=int, default=5,
        help=('The minimum number of reads that must cover a location for it '
              'to be considered significant.'))

    parser.add_argument(
        '--homogeneousCutoff', type=float, default=0.9,
        help=('If the most common nucleotide at a location occurs more than '
              'this fraction of the time (i.e., amongst all reads that cover '
              'the location) then the location will be considered homogeneous '
              'and therefore uninteresting.'))


def addCommandLineOptions(parser, outfileDefaultName=None):
    """
    Add standard command-line options to an argument parser.

    @param parser: An C{ArgumentParser} instance.
    @param outfileDefaultName: The C{str} output file to use as a default
        in case the user does not give one on the command line.
    """

    addCommonOptions(parser)
    SAMFilter.addFilteringOptions(parser)

    parser.add_argument(
        '--outfile', default=outfileDefaultName,
        help='The filename to store the resulting HTML.')

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
    samFilter = SAMFilter.parseFilteringOptions(args)

    if samFilter.referenceIds and len(samFilter.referenceIds) > 1:
        raise ValueError('Only one reference id can be given.')

    referenceLengths = samFilter.referenceLengths()

    if len(referenceLengths) == 1:
        referenceId, genomeLength = referenceLengths.popitem()
    else:
        raise ValueError(
            'If you do not specify a reference sequence with '
            '--referenceId, the SAM/BAM file must contain exactly one '
            'reference. But %s contains %d (%s).' %
            (args.samfile, len(referenceLengths)))

    paddedSAM = PaddedSAM(samFilter)

    for query in paddedSAM.queries():
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


def addAnalysisCommandLineOptions(parser):
    """
    Add command-line options used in a read analysis.
    """

    addCommonOptions(parser)

    parser.add_argument(
        '--referenceGenome', metavar='FILENAME', action='append', nargs='+',
        required=True,
        help=('The name of a FASTA file containing reference genomes that '
              'were used to create the alignment files (may be repeated).'))

    parser.add_argument(
        '--alignmentFile', metavar='FILENAME', action='append', nargs='+',
        required=True,
        help='The name of a SAM/BAM alignment file (may be repeated).')

    parser.add_argument(
        '--referenceId', metavar='NAME', action='append', nargs='*',
        help=('The sequence id whose alignment should be analyzed (may '
              'be repeated). All ids must be present in --referenceGenome '
              'file. One of the SAM/BAM files given using --alignmentFile '
              'should have an alignment against the given argument. If '
              'omitted, all references that are aligned to in the given '
              'BAM/SAM files will be analyzed.'))

    parser.add_argument(
        '--outputDir',
        help='The directory to save result files to.')

    parser.add_argument(
        '--saveReducedFASTA', default=False, action='store_true',
        help=('If given, write out a FASTA file of the original input but '
              'with just the signifcant locations.'))

    parser.add_argument(
        '--plotSAM', default=False, action='store_true',
        help=('If given, save plots showing where reads are aligned to on '
              'the genome along with their alignment scores.'))

    parser.add_argument(
        '--verbose', type=int, default=0,
        help=('The integer verbosity level (0 = no output, 1 = some output, '
              'etc).'))
