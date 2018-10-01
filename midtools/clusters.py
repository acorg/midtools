from __future__ import division

from collections import defaultdict
from itertools import count

from midtools.distances import DistanceCache
from midtools.offsets import OffsetBases
from midtools.utils import nucleotidesToStr, s


class ReadCluster(object):
    """
    Maintain a cluster of reads.
    """

    MIN_COMMONEST_MULTIPLE = 10.0

    def __init__(self):
        self.nucleotides = defaultdict(OffsetBases)
        self.reads = set()

    def __len__(self):
        return len(self.reads)

    def __str__(self):
        result = ['Cluster with %d read%s:' %
                  (len(self.reads), s(len(self.reads)))]
        for read in self.reads:
            result.append('  %s' % read)
        result.append(nucleotidesToStr(self.nucleotides, prefix='  '))
        return '\n'.join(result)

    def add(self, read):
        """
        Add a read to this cluster.

        @param read: An C{alignedRead} instance.
        """
        assert read not in self.reads
        self.reads.add(read)

        nucleotides = self.nucleotides
        for offset, base in read.significantOffsets.items():
            nucleotides[offset].incorporateBase(base)

    def merge(self, other):
        """
        Merge another cluster into this one.

        @param other: An C{ReadCluster} instance.
        """
        # Sanity check that the two clusters are distinct.
        assert not (self.reads & other.reads)
        self.reads.update(other.reads)

        for offset, offsetBases in other.nucleotides.items():
            self.nucleotides[offset].merge(offsetBases)

    @staticmethod
    def commonOffsetsMaxFraction(a, b):
        """
        Compute the maximum fraction of the number of common offsets
        (between two clusters) and the number of offsets covered by each
        cluster.

        @param a: A C{ReadCluster} instance.
        @param b: A C{ReadCluster} instance.
        @raise ZeroDivisionError: if C{a} and C{b} have no offsets (neither of
            which should be possible in normal operation).
        @return: The C{float} maximum fraction.
        """
        common = set(a.nucleotides) & set(b.nucleotides)
        return len(common) / min(len(a.nucleotides), len(b.nucleotides))

    @staticmethod
    def commonNucleotidesMultiplicativeDistance(a, b):
        """
        Measure the distance from one cluster to another, as the sum of the
        multiplied probabilities of nucleotides.

        @param a: A C{ReadCluster} instance.
        @param b: A C{ReadCluster} instance.
        @raise ZeroDivisionError: if C{a} or C{b} has no offsets (neither of
            which should be possible in normal operation).
        @return: The C{float} [0.0, 1.0] distance between C{a} and C{b}.
        """
        aNucleotides = a.nucleotides
        bNucleotides = b.nucleotides
        commonOffsets = set(aNucleotides) & set(bNucleotides)

        if commonOffsets:
            similarity = sum(
                1.0 - min(
                    OffsetBases.multiplicativeDistance(
                        aNucleotides[offset], bNucleotides[offset]),
                    OffsetBases.homogeneousDistance(
                        aNucleotides[offset], bNucleotides[offset]))
                for offset in commonOffsets)
            return 1.0 - (similarity / len(commonOffsets))
        else:
            return 1.0

    @staticmethod
    def commonNucleotidesAgreementDistance(a, b):
        """
        Measure the distance from one cluster to another, according to how
        often the intersection of the commonest nucleotides (at each site)
        is non-empty. The distance is 1.0 minus the fraction of common sites
        with a non-empty intersection.

        @param a: A C{ReadCluster} instance.
        @param b: A C{ReadCluster} instance.
        @return: The C{float} [0.0, 1.0] distance between C{a} and C{b}.
        """
        aNucleotides = a.nucleotides
        bNucleotides = b.nucleotides
        commonOffsets = set(aNucleotides) & set(bNucleotides)

        if commonOffsets:
            matching = 0
            for offset in commonOffsets:
                aNucleotidesAtOffset = aNucleotides[offset]
                bNucleotidesAtOffset = bNucleotides[offset]
                if (aNucleotidesAtOffset.commonest &
                        bNucleotidesAtOffset.commonest):
                    matching += 1
                else:
                    multiple = OffsetBases.highestFrequenciesMultiple(
                        aNucleotidesAtOffset, bNucleotidesAtOffset)
                    # Sanity: the multiple cannot be None because that
                    # would mean only one nucleotide is present, and that
                    # case is dealt with by the first part of this if/then.
                    assert multiple is not None
                    if multiple >= ReadCluster.MIN_COMMONEST_MULTIPLE:
                        matching += 1

            return 1.0 - (matching / len(commonOffsets))
        else:
            return 1.0


class ReadClusters(object):
    """
    Maintain clusters of reads.
    """

    COMMON_OFFSETS_MAX_FRACTION_MIN = 0.9

    def __init__(self):
        self._count = count()
        self.readClusters = defaultdict(ReadCluster)
        # self.distanceCache = DistanceCache(self.multiplicativeDistance)
        self.distanceCache = DistanceCache(
            self.commonNucleotidesAgreementDistance)

    def __len__(self):
        """
        Our length is the number of read clusters.

        @return: An C{int} length.
        """
        return len(self.readClusters)

    def add(self, read):
        """
        Add a single read as a new cluster.

        @param read: An C{alignedRead} instance.
        @return: The C{int} number of the new cluster.
        """
        count = next(self._count)
        self.readClusters[count].add(read)
        return count

    def commonNucleotidesAgreementDistance(self, a, b):
        """
        Measure the distance from one cluster to another based on the fraction
        of shared offsets where the commonest nucleotide(s) for the offset
        has a non-empty intersection. Hard to explain...

        @param a: An C{int} cluster number.
        @param b: An C{int} cluster number.
        @return: The C{float} [0.0, 1.0] distance between C{a} and C{b}.
        """
        a = self.readClusters[a]
        b = self.readClusters[b]

        return 1.0 - (
            (1.0 - ReadCluster.commonNucleotidesAgreementDistance(a, b))
            *
            max(self.COMMON_OFFSETS_MAX_FRACTION_MIN,
                ReadCluster.commonOffsetsMaxFraction(a, b)))

    def multiplicativeDistance(self, a, b):
        """
        Calculate an inter-cluster distance based on the maximum fraction of
        offsets (for the two clusters) present in common and the
        (multiplicative) degree to which they match.
        """
        a = self.readClusters[a]
        b = self.readClusters[b]

        # The following is 1.0 minus the product of two similarities,
        # giving a distance. We don't let the similarity penalty for
        # covered offset fraction be greater than 0.75.
        return 1.0 - (
            (1.0 - ReadCluster.commonNucleotidesMultiplicativeDistance(a, b))
            *
            max(self.COMMON_OFFSETS_MAX_FRACTION_MIN,
                ReadCluster.commonOffsetsMaxFraction(a, b)))

    def mergeDescription(self, a, b, distance):
        """
        Make a textual description of a cluster merge.

        @param a: An C{int} cluster number.
        @param b: An C{int} cluster number.
        @param distance: The C{float} [0.0, 1.0] distance between the clusters.
        @return: A C{str} side-by-side descriptions of clusters C{a} and C{b}.
        """
        cluster1 = self.readClusters[a]
        cluster2 = self.readClusters[b]

        result1 = []
        result2 = []
        matches = []
        sharedCount = matchCount = 0

        allOffsets = sorted(
            set(cluster1.nucleotides) | set(cluster2.nucleotides))

        for offset in allOffsets:

            inCount = 0

            if offset in cluster1.nucleotides:
                result1.append(cluster1.nucleotides[offset].baseCountsToStr())
                inCount += 1
            else:
                result1.append('-')

            if offset in cluster2.nucleotides:
                result2.append(cluster2.nucleotides[offset].baseCountsToStr())
                inCount += 1
            else:
                result2.append('-')

            if inCount == 2:
                sharedCount += 1
                if (cluster1.nucleotides[offset].commonest &
                        cluster2.nucleotides[offset].commonest):
                    matches.append('*')
                    matchCount += 1
                else:
                    multiple = OffsetBases.highestFrequenciesMultiple(
                        cluster1.nucleotides[offset],
                        cluster2.nucleotides[offset])
                    # Sanity: the multiple cannot be None because that
                    # would mean only one nucleotide is present, and that
                    # case is dealt with by the first part of this if/then.
                    assert multiple is not None
                    if multiple >= ReadCluster.MIN_COMMONEST_MULTIPLE:
                        matchCount += 1
                        matches.append('+')
                    else:
                        matches.append('')
            else:
                matches.append('')

        result1Width = max(len(line) for line in result1)
        result2Width = max(len(line) for line in result2)

        return '\n'.join(
            [
                ('Merging clusters %d and %d with distance %.2f' %
                 (a, b, distance)),
                ('Cluster %d has %d read%s, covering %d offset%s' %
                 (a,
                  len(cluster1.reads), s(len(cluster1.reads), ),
                  len(cluster1.nucleotides), s(len(cluster1.nucleotides)))),
                ('Cluster %d has %d read%s, covering %d offset%s' %
                 (b,
                  len(cluster2.reads), s(len(cluster2.reads)),
                  len(cluster2.nucleotides), s(len(cluster2.nucleotides)))),
                ('%d matches out of %d shared offsets' %
                 (matchCount, sharedCount)),
            ] + [
                '  %d: %*s    %*s    %s' %
                (offset, result1Width, line1, result2Width, line2, match) for
                (offset, line1, line2, match) in
                zip(allOffsets, result1, result2, matches)
            ]
        )

    def mergeDescriptionWithOffsetScores(self, a, b, distance):
        """
        Make a textual description of a cluster merge, including per-offset
        score information.

        @param a: An C{int} cluster number.
        @param b: An C{int} cluster number.
        @param distance: The C{float} [0.0, 1.0] distance between the clusters.
        @return: A C{str} side-by-side descriptions of clusters C{a} and C{b}.
        """
        cluster1 = self.readClusters[a]
        cluster2 = self.readClusters[b]

        result1 = []
        result2 = []
        offsetScores = []
        matches = []
        sharedCount = matchCount = 0

        allOffsets = sorted(
            set(cluster1.nucleotides) | set(cluster2.nucleotides))

        for offset in allOffsets:

            inCount = 0

            if offset in cluster1.nucleotides:
                result1.append(cluster1.nucleotides[offset].baseCountsToStr())
                inCount += 1
            else:
                result1.append('-')

            if offset in cluster2.nucleotides:
                result2.append(cluster2.nucleotides[offset].baseCountsToStr())
                inCount += 1
            else:
                result2.append('-')

            if inCount == 2:
                sharedCount += 1
                if (cluster1.nucleotides[offset].commonest &
                        cluster2.nucleotides[offset].commonest):
                    matches.append('*')
                    matchCount += 1
                else:
                    matches.append('')

                offsetScores.append(
                    '%.3f' % min(
                        OffsetBases.multiplicativeDistance(
                            cluster1.nucleotides[offset],
                            cluster2.nucleotides[offset]),
                        OffsetBases.homogeneousDistance(
                            cluster1.nucleotides[offset],
                            cluster2.nucleotides[offset])))
            else:
                matches.append('')
                offsetScores.append('')

        result1Width = max(len(line) for line in result1)
        result2Width = max(len(line) for line in result2)
        offsetScoresWidth = max(len(line) for line in offsetScores)

        return '\n'.join(
            [
                ('Merging clusters %d and %d with distance %.2f' %
                 (a, b, distance)),
                ('Cluster %d has %d read%s, covering %d offset%s' %
                 (a,
                  len(cluster1.reads), s(len(cluster1.reads), ),
                  len(cluster1.nucleotides), s(len(cluster1.nucleotides)))),
                ('Cluster %d has %d read%s, covering %d offset%s' %
                 (b,
                  len(cluster2.reads), s(len(cluster2.reads)),
                  len(cluster2.nucleotides), s(len(cluster2.nucleotides)))),
                ('%d matches out of %d shared offsets' %
                 (matchCount, sharedCount)),
            ] + [
                '  %d: %*s    %*s    %*s    %s' %
                (offset, result1Width, line1, result2Width, line2,
                 offsetScoresWidth, offsetScore, match) for
                (offset, line1, line2, offsetScore, match) in
                zip(allOffsets, result1, result2, offsetScores, matches)
            ]
        )

    def analyze(self, cutoff, fp=None):
        """
        Perform the cluster analysis, up to a given distance cut off.

        @param cutoff: The C{float} distance at which clustering will be
            stopped.
        @param fp: A file-like object to write information to, or C{None}
            if no output should be produced.
        @return: A generator that yields C{ReadCluster} instances, being the
            clusters that were too distant (according to C{cutoff} to be
            merged.
        """
        if fp:
            print('Starting cluster analysis on %d reads' % len(self), file=fp)

        for a in self.readClusters:
            self.distanceCache.add(a)

        while True:
            lowestDistance = self.distanceCache.lowestDistance()
            if lowestDistance is None or lowestDistance > cutoff:
                break

            a, b = self.distanceCache.pop()
            self.distanceCache.remove(a)
            self.distanceCache.remove(b)

            if fp:
                print(self.mergeDescription(a, b, lowestDistance), file=fp)

            # Merge.
            count = next(self._count)
            self.readClusters[a].merge(self.readClusters[b])
            self.readClusters[count] = self.readClusters[a]
            self.distanceCache.add(count)

            # Remove the merged clusters.
            del self.readClusters[a]
            del self.readClusters[b]

            if fp:
                print('New cluster (%d) has %d reads and covers %d offsets. '
                      'Remaining cluster count: %d.' %
                      (count, len(self.readClusters[count].reads),
                       len(self.readClusters[count].nucleotides), len(self)),
                      file=fp)

        if fp:
            if lowestDistance is None:
                print('Distance priority queue is empty, all reads '
                      'are now in one cluster.', file=fp)
            else:
                print('The lowest inter-cluster distance has risen '
                      'to %.2f, which exceeds the %.2f cutoff value. '
                      '%d unmerged clusters remain.' %
                      (lowestDistance, cutoff, len(self.readClusters)),
                      file=fp)

                # Print the remaining inter-cluster distances.
                print('\nThe remaining cluster differences that were too '
                      'distant to qualify are:', file=fp)

                while True:
                    lowestDistance = self.distanceCache.lowestDistance()
                    if lowestDistance is None:
                        break

                    a, b = self.distanceCache.pop()
                    # self.distanceCache.remove(a)
                    # self.distanceCache.remove(b)

                    print(self.mergeDescription(a, b, lowestDistance), file=fp)

        if fp:
            print('End of cluster analysis.', file=fp)

        return self.readClusters.values()
