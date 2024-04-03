from collections import Counter

from midtools.utils import baseCountsToStr


class OffsetBases:
    """
    Maintain the count of nucleotide bases at an offset.
    """

    def __init__(self):
        self._counts = Counter()
        self._commonest = None
        self._clean = True

    def __eq__(self, other):
        """
        Are two instances equal?

        @param other: Another C{OffsetBases} instance.
        @return: A Boolean indicating equality.
        """
        return self._counts == other._counts

    def __str__(self):
        return self.baseCountsToStr()

    def __repr__(self):
        return "<%s %s>" % (self.__class__.__name__, self.baseCountsToStr())

    def incorporateBase(self, base):
        """
        Incorporate a new instance of a base at an offset.

        @param base: A C{str} base.
        """
        self._counts[base] += 1
        self._clean = False

    def unincorporateBase(self, base):
        """
        Remove an instance of a base at an offset.

        @param base: A C{str} base.
        """
        self._counts[base] -= 1
        # Sanity check.
        assert self._counts[base] >= 0
        self._clean = False

    @property
    def commonest(self):
        """
        Find the commonest bases.

        @return: Either C{None} if no bases have been added or else a C{set}
            of C{str} bases that are most common.
        """
        if not self._clean:
            self._recalculate()
        return self._commonest

    def _recalculate(self):
        """
        Re-calculate the set of commonest bases.
        """
        orderedCounts = self._counts.most_common()
        maxCount = orderedCounts[0][1]
        self._commonest = set(x[0] for x in orderedCounts if x[1] == maxCount)
        self._clean = True

    def baseCountsToStr(self):
        """
        Convert base counts to a string.

        @return: A C{str} representation of nucleotide counts.
        """
        return baseCountsToStr(self._counts)

    def merge(self, other):
        """
        Merge in the counts from another instance.

        @param other: An C{OffsetBases} base.
        """
        self._counts += other._counts
        self._clean = False

    @staticmethod
    def multiplicativeDistance(a, b):
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
    def homogeneousDistance(a, b):
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
    def highestFrequenciesMultiple(a, b):
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
