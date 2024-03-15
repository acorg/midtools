from __future__ import division

from collections import Counter


def gatherData(genomeLength, alignedReads):
    """
    Analyze the aligned reads.

    @param genomeLength: The C{int} length of the genome the reads were
        aligned to.
    @param alignedReads: A C{list} of C{AlignedRead} instances.
    @return: A tuple of C{list}s (readCountAtOffset, baseCountAtOffset,
        readsAtOffset), each indexed from zero to the genome length.
    """
    readCountAtOffset = []
    baseCountAtOffset = []
    readsAtOffset = []

    nucleotides = set("ACGT")

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


def findSignificantOffsets(
    baseCountAtOffset, readCountAtOffset, minReads, homogeneousCutoff
):
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
        zip(readCountAtOffset, baseCountAtOffset)
    ):
        if (
            readCount >= minReads
            and max(counts.values()) / readCount <= homogeneousCutoff
        ):
            yield offset
