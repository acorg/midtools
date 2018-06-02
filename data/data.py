from __future__ import division

from collections import Counter

from dark.fasta import FastaReads

from data.read import AlignedRead


def gatherData(genomeLength, alignedReads):
    """
    Analyze the aligned reads.

    @param genomeLength: The C{int} length of the genome the reads were
        aligned to.
    @param alignedReads: A C{list} of C{AlignedRead} instances.
    @return: A tuple of C{list}s (readCountAtOffset, baseCountAtOffset,
        readsAtOffset), each indexed from zero to the genom length.
    """
    readCountAtOffset = []
    baseCountAtOffset = []
    readsAtOffset = []

    nucleotides = set('ACGT')

    for offset in range(genomeLength):
        reads = set()
        counts = Counter()
        for read in alignedReads:
            base = read.base(offset)
            if base in nucleotides:
                counts[base] += 1
                reads.add(read)
        baseCountAtOffset.append(counts)
        readCountAtOffset.append(sum(counts.values()))
        readsAtOffset.append(reads)

    return readCountAtOffset, baseCountAtOffset, readsAtOffset


def findSignificantOffsets(baseCountAtOffset, readCountAtOffset,
                           minReads, homogeneousCutoff):
    """
    Find the genome offsets that have significant base variability.

    @param baseCountAtOffset: A C{list} of C{Counter} instances giving
        the count of each nucleotide at each genome offset.
    @param readCountAtOffset: A C{list} of C{int} counts of the total
        number of reads at each genome offset (i.e., just the sum of the
        values in C{baseCountAtOffset})
    @param minReads: The C{int} minimum number of reads that must cover
        a offset for it to be considered significant.
    @param homogeneousCutoff: A C{float} frequency. If the most common
        nucleotide at a offset occurs *more than* this fraction of the time
        (i.e., amongst all reads that cover the offset) then the locaion
        will be considered homogeneous and therefore uninteresting.
    @return: A generator that yields 0-based significant offsets.
    """
    for offset, (readCount, counts) in enumerate(
            zip(readCountAtOffset, baseCountAtOffset)):
        if (readCount >= minReads and
                max(counts.values()) / readCount <= homogeneousCutoff):
            yield offset


def addCommandLineOptions(parser, outfileDefaultName=None):
    """
    Add standard command-line options to an argument parser.

    @param parser: An C{ArgumentParser} instance.
    @param outfileDefaultName: The C{str} output file to use as a default
        in case the user does not give one on the command line.
    """
    parser.add_argument(
        '--fastaFile', metavar='FILENAME', required=True,
        help='The name of the FASTA input file.')

    parser.add_argument(
        '--outFile', default=outfileDefaultName,
        help='The filename to store the resulting HTML.')

    parser.add_argument(
        '--minReads', type=int, default=5,
        help=('The minimum number of reads that must cover a location for it '
              'to be considered significant.'))

    parser.add_argument(
        '--trim', type=int, default=0,
        help=('The number of bases to trim from the start and end of each '
              'read (to try to remove possible damage).'))

    parser.add_argument(
        '--homogeneousCutoff', type=float, default=0.9,
        help=('If the most common nucleotide at a location occurs more than '
              'this fraction of the time (i.e., amongst all reads that cover '
              'the location) then the locaion will be considered homogeneous '
              'and therefore uninteresting.'))

    parser.add_argument(
        '--show', action='store_true', default=False,
        help='If specified, show the figure interactively.')


def parseCommandLineOptions(args, findSigOffsets):
    """
    Deal with the various command-line options added to the ArgumentParser
    instance by addCommandLineOptions.

    @param args: The result of calling C{parse_args} on an C{ArgumentParser}
        instance (the one that was passed to C{addCommandLineOptions}, unless
        we're testing).
    @param findSigOffsets: If C{True} also return a list of the significant
        offsets (else that element of the return value will be C{None}).
    @return: A 5-tuple: (genomeLength, alignedReads, readCountAtOffset,
        baseCountAtOffset, readsAtOffset, significantOffsets)
    """
    genomeLength = None
    reads = FastaReads(args.fastaFile)
    alignedReads = []
    trim = args.trim
    for count, read in enumerate(reads, start=1):
        if count == 1:
            genomeLength = len(read)
        else:
            assert len(read) == genomeLength, (
                'Read number %d with id %r in %s had unexpected length (%d '
                'instead of %d)' %
                (count, read.id, args.fastaFile, len(read), genomeLength))
        try:
            ar = AlignedRead(read.id, read.sequence)
        except ValueError as e:
            # Ignore reads that are all gaps.
            assert str(e) == 'Read is all gaps.'
        else:
            if trim:
                if ar.trim(trim):
                    alignedReads.append(ar)
            else:
                alignedReads.append(ar)

    readCountAtOffset, baseCountAtOffset, readsAtOffset = gatherData(
        genomeLength, alignedReads)

    if findSigOffsets:
        significantOffsets = list(findSignificantOffsets(
            baseCountAtOffset, readCountAtOffset, args.minReads,
            args.homogeneousCutoff))
        for read in alignedReads:
            read.setSignificantOffsets(significantOffsets)
    else:
        significantOffsets = None

    return (genomeLength, alignedReads, readCountAtOffset,
            baseCountAtOffset, readsAtOffset, significantOffsets)
