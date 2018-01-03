from collections import Counter

from dark.reads import (
    addFASTACommandLineOptions, parseFASTACommandLineOptions)
from dark.fasta import FastaReads


class AlignedRead(object):
    """
    Hold information about a read that has been aligned to a consensus.

    @param read: A C{dark.reads.Read} instance.
    """
    def __init__(self, read):
        sequence = read.sequence.upper()

        # Scan for initial gaps.
        offset = 0
        for base in sequence:
            if base == '-':
                offset += 1
            else:
                break

        # Scan for final gaps.
        trailing = 0
        for base in sequence[::-1]:
            if base == '-':
                trailing += 1
            else:
                break

        # Make sure the read is not all gaps.
        assert offset + trailing < len(sequence)

        read.sequence = sequence[offset:len(sequence) - trailing]
        self.read = read
        self.offset = offset

    def __str__(self):
        return '<AlignedRead: (offset %4d, len %2d) %s>' % (
            self.offset, len(self.read), self.read.id)

    def base(self, n):
        """
        Get the nucleotide base at a given offset.

        @param n: An C{int} offset on the genome.
        @return: The C{str} nucleotide, or C{None} if the read does not cover
            the genome at that offset.
        """
        offset = self.offset
        if n >= offset and n < offset + len(self.read):
            return self.read.sequence[n - offset]

    def trim(self, n):
        """
        Trim bases from the start and end of the read.

        @param n: The C{int} number of bases to remove from each end.
        @return: A C{bool} to indicate whether the trimming was performed or
            not (due to the read being too short).
        """
        assert n >= 0, ('Trim amount (%d) cannot be negative.' % n)
        if 2 * n < len(self.read):
            self.read = self.read[n:len(self.read) - n]
            self.offset += n
            return True
        else:
            return False


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
    addFASTACommandLineOptions(parser)

    parser.add_argument(
        '--genomeFile', required=True,
        help=('The filename of the FASTA file containing the supposed genome '
              'that the reads came from.'))

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
        '--homogeneousCutoff', type=float, default=0.85,
        help=('If the most common nucleotide at a location occurs more than '
              'this fraction of the time (i.e., amongst all reads that cover '
              'the location) then the locaion will be considered homogeneous '
              'and therefore uninteresting.'))

    parser.add_argument(
        '--show', action='store_true', default=False,
        help='If specified, show the figure interactively.')


def parseCommandLineOptions(args):
    """
    Deal with the various command-line options added to the ArgumentParser
    instance by addCommandLineOptions.

    @param args: The result of calling C{parse_args} on an C{ArgumentParser}
        instance (the one that was passed to C{addCommandLineOptions}, unless
        we're testing).
    @return: A 5-tuple: (genomeLength, alignedReads, readCountAtOffset,
        baseCountAtOffset, readsAtOffset)
    """
    genomeReads = list(FastaReads(args.genomeFile))
    assert len(genomeReads) == 1, (
        'The genome file %s has %d sequences in it. Expected just one.' %
        len(genomeReads))
    genome = genomeReads[0]
    genomeLength = len(genome)

    reads = parseFASTACommandLineOptions(args)
    alignedReads = []
    trim = args.trim
    for count, read in enumerate(reads, start=1):
        assert len(read) == genomeLength, (
            'Read number %d with id %r in %s had unexpected length (%d '
            'instead of %d)' %
            (count, read.id, args.fastaFile, len(read), genomeLength))
        ar = AlignedRead(read)
        if trim:
            if ar.trim(trim):
                alignedReads.append(ar)
        else:
            alignedReads.append(ar)

    readCountAtOffset, baseCountAtOffset, readsAtOffset = gatherData(
        genomeLength, alignedReads)

    return (genome, alignedReads, readCountAtOffset, baseCountAtOffset,
            readsAtOffset)
