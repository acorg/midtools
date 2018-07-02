from collections import OrderedDict

from dark.reads import Read


class AlignedRead(Read):
    """
    Hold information about a read that has been aligned to a consensus.

    @param read: A C{dark.reads.Read} instance.
    """
    def __init__(self, id_, sequence):
        self.significantOffsets = OrderedDict()
        self._originalLength = len(sequence)

        # Scan the sequence for initial gaps.
        offset = 0
        for base in sequence:
            if base == '-':
                offset += 1
            else:
                break

        if offset == len(sequence):
            raise ValueError('Read is all gaps.')

        # Scan for final gaps.
        trailing = 0
        for base in sequence[::-1]:
            if base == '-':
                trailing += 1
            else:
                break

        # Make sure the read is not all gaps.
        assert offset + trailing < len(sequence)
        self.offset = offset

        Read.__init__(
            self, id_, sequence[offset:len(sequence) - trailing].upper())

    def __str__(self):
        if self.significantOffsets:
            bases = ', bases %s, offsets (total %d): %s' % (
                ''.join(self.significantOffsets.values()),
                len(self.significantOffsets),
                ','.join(map(str, self.significantOffsets)))
        else:
            bases = ''

        return '<AlignedRead: (offset %4d, len %2d%s) %s>' % (
            self.offset, len(self), bases, self.id)

    def __lt__(self, other):
        """
        Sort order is according to number of signifcant offsets (decreasing),
        then the offset where the alignment begins, then the normal
        C{dark.reads.Read} sort order.
        """
        t1 = (-len(self.significantOffsets), self.offset)
        t2 = (-len(other.significantOffsets), other.offset)
        if t1 == t2:
            return Read.__lt__(self, other)
        else:
            return t1 < t2

    def agreesWith(self, other, agreementFraction):
        """
        Two reads agree if they have identical bases in a sufficiently high
        fraction of their shared significant offsets.

        @param other: Another C{AlignedRead} instance.
        @param agreementFraction: A [0..1] C{float} fraction. Agreement is true
            if the fraction of identical bases is at least this high.
        @return: C{True} if the reads agree, C{False} if not.
        """
        sharedCount = identicalCount = 0
        getOtherBase = other.significantOffsets.get
        for offset, base in self.significantOffsets.items():
            otherBase = getOtherBase(offset)
            if otherBase:
                sharedCount += 1
                identicalCount += (otherBase == base)
        if sharedCount:
            return (identicalCount / sharedCount) >= agreementFraction
        else:
            return True

    def setSignificantOffsets(self, significantOffsets):
        """
        Find the base at each of the significant offsets covered by this read.

        @param significantOffsets: A C{list} of C{int} offsets.
        """
        newSignificantOffsets = OrderedDict()
        for offset in significantOffsets:
            base = self.base(offset)
            if base is not None:
                # Note that we cannot break out of this loop early if base
                # is None because some reads have embedded gaps ('-'), for
                # which self.base returns None. So we have to continue on
                # to higher offsets.
                newSignificantOffsets[offset] = base
        self.significantOffsets = newSignificantOffsets

    def base(self, n):
        """
        Get the nucleotide base at a given offset.

        @param n: An C{int} offset on the genome.
        @return: The C{str} nucleotide, or C{None} if the read does not cover
            the genome at that offset.
        """
        offset = self.offset
        if n >= offset and n < offset + len(self):
            b = self.sequence[n - offset]
            return None if b == '-' else b

    def trim(self, n):
        """
        Trim bases from the start and end of the read.

        @param n: The C{int} number of bases to remove from each end.
        @return: A C{bool} to indicate whether the trimming was performed or
            not (due to the read being too short).
        """
        assert n >= 0, ('Trim amount (%d) cannot be negative.' % n)
        if 2 * n < len(self):
            self.sequence = self.sequence[n:len(self) - n]
            self.offset += n
            return True
        else:
            return False

    def toString(self):
        """
        Make a FASTA string for the read, including its original gaps.
        """
        return '>%s\n%s%s%s\n' % (
            self.id,
            '-' * self.offset,
            self.sequence,
            '-' * (self._originalLength - self.offset - len(self.sequence)))
