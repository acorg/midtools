from __future__ import annotations

from collections import Counter
from typing import Iterator, Optional

from midtools.read import AlignedRead


class OffsetBases:
    """
    Maintain the count of nucleotide bases at an offset.
    """

    def __init__(self) -> None:
        self._counts: Counter[str] = Counter()
        self._commonest: Optional[set[str]] = None
        self._clean = True

    def __eq__(self, other: object) -> bool:
        """
        Are two instances equal?

        @param other: Another C{OffsetBases} instance.
        @return: A Boolean indicating equality.
        """
        if not isinstance(other, OffsetBases):
            return NotImplemented

        return self._counts == other._counts

    def __str__(self) -> str:
        return self.baseCountsToStr()

    def __repr__(self):
        return "<%s %s>" % (self.__class__.__name__, self.baseCountsToStr())

    def incorporateBase(self, base: str) -> None:
        """
        Incorporate a new instance of a base at an offset.

        @param base: A C{str} base.
        """
        self._counts[base] += 1
        self._clean = False

    def unincorporateBase(self, base: str) -> None:
        """
        Remove an instance of a base at an offset.

        @param base: A C{str} base.
        """
        self._counts[base] -= 1
        # Sanity check.
        assert self._counts[base] >= 0
        self._clean = False

    @property
    def commonest(self) -> Optional[set[str]]:
        """
        Find the commonest bases.

        @return: Either C{None} if no bases have been added or else a C{set}
            of C{str} bases that are most common.
        """
        if not self._clean:
            self._recalculate()
        return self._commonest

    def _recalculate(self) -> None:
        """
        Re-calculate the set of commonest bases.
        """
        orderedCounts = self._counts.most_common()
        maxCount = orderedCounts[0][1]
        self._commonest = set(x[0] for x in orderedCounts if x[1] == maxCount)
        self._clean = True

    def baseCountsToStr(self) -> str:
        """
        Convert base counts to a string.

        @return: A C{str} representation of nucleotide counts.
        """
        # Avoid circular import.
        from midtools.utils import baseCountsToStr

        return baseCountsToStr(self._counts)

    def merge(self, other: OffsetBases) -> None:
        """
        Merge in the counts from another instance.

        @param other: An C{OffsetBases} base.
        """
        self._counts += other._counts
        self._clean = False

    @staticmethod
    def multiplicativeDistance(a: OffsetBases, b: OffsetBases) -> float:
        """
        Measure the multiplicative distance from one set of offsets to another,
        as 1.0 minus the sum of the multiplied probabilities of their common
        nucleotides.

        @param a: An C{OffsetBases} instance.
        @param b: An C{OffsetBases} instance.
        @raise ZeroDivisionError: if C{a} or C{b} have no nucleotides (neither
            of which should be possible in normal operation).
        @return: The C{float} [0.0, 1.0] distance between C{a} and C{b}.
        """
        aCounts = a._counts
        bCounts = b._counts
        aTotal = sum(aCounts.values())
        bTotal = sum(bCounts.values())
        commonBases = set(aCounts) & set(bCounts)

        # Let a ZeroDivisionError occur if aTotal or bTotal is zero.
        return 1.0 - sum(
            (aCounts[base] / aTotal) * (bCounts[base] / bTotal) for base in commonBases
        )

    @staticmethod
    def homogeneousDistance(a: OffsetBases, b: OffsetBases) -> float:
        """
        Measure the homogeneous distance from one set of nucleotides to
        another.

        @param a: An C{OffsetBases} instance.
        @param b: An C{OffsetBases} instance.
        @raise ZeroDivisionError: if C{a} and C{b} have no nucleotides (neither
            of which should be possible in normal operation).
        @return: The C{float} [0.0, 0.75] distance between C{a} and C{b}.
        """
        aCounts = a._counts
        bCounts = b._counts
        denom = sum(aCounts.values()) + sum(bCounts.values())
        allBases = set(aCounts) | set(bCounts)

        # Let a ZeroDivisionError occur if denom is zero.
        return 1.0 - max(((aCounts[base] + bCounts[base]) / denom) for base in allBases)

    @staticmethod
    def highestFrequenciesMultiple(a: OffsetBases, b: OffsetBases) -> Optional[float]:
        """
        How much does the most frequent nucleotide occur more than the second
        most?

        @param a: An C{OffsetBases} instance.
        @param b: An C{OffsetBases} instance.
        @return: A C{float} made from the divsion of the count of the most
            frequent nucleotide with the count of the second most in the
            combined counts. If there is only one nucleotide, return C{None}.
        """
        counts = a._counts + b._counts
        if len(counts) == 1:
            # There is only one nucleotide.
            return None
        else:
            orderedCounts = counts.most_common()
            return orderedCounts[0][1] / orderedCounts[1][1]


def analyzeOffets(
    genomeLength: int, alignedReads: list[AlignedRead]
) -> tuple[list[int], list[Counter[str]], list[set[AlignedRead]]]:
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
        counts: Counter[str] = Counter()
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
    baseCountAtOffset: list[Counter[str]],
    readCountAtOffset: list[int],
    minReads: int,
    homogeneousCutoff: float,
) -> Iterator[int]:
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
