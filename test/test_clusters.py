from unittest import TestCase

from midtools.clusters import ReadCluster, ReadClusters
from midtools.offsets import OffsetBases
from midtools.read import AlignedRead


class TestReadCluster(TestCase):
    """
    Test the ReadCluster class.
    """

    def testEmpty(self):
        """
        An empty instance must have length zero.
        """
        self.assertEqual(0, len(ReadCluster()))

    def testAddOne(self):
        """
        Adding a read must result in the expected nucleotides and read
        being stored.
        """
        read = AlignedRead("id", "---ACGT--")
        read.setSignificantOffsets([3])
        rc = ReadCluster()
        rc.add(read)
        expectedBases = OffsetBases()
        expectedBases.incorporateBase("A")
        self.assertEqual({read}, rc.reads)
        self.assertEqual({3: expectedBases}, rc.nucleotides)
        self.assertEqual(1, len(rc))

    def testAddTwo(self):
        """
        Adding two reads must result in the expected nucleotides and reads
        being stored.
        """
        read1 = AlignedRead("id1", "---ACCC--")
        read1.setSignificantOffsets([3, 4])

        read2 = AlignedRead("id2", "---TG--")
        read2.setSignificantOffsets([3])

        rc = ReadCluster()
        rc.add(read1)
        rc.add(read2)

        self.assertEqual({read1, read2}, rc.reads)

        expectedBases3 = OffsetBases()
        expectedBases3.incorporateBase("A")
        expectedBases3.incorporateBase("T")
        expectedBases4 = OffsetBases()
        expectedBases4.incorporateBase("C")

        self.assertEqual({3: expectedBases3, 4: expectedBases4}, rc.nucleotides)

    def testMerge(self):
        """
        Merging two clusters must result in the expected nucleotides and reads
        being stored.
        """
        read1 = AlignedRead("id1", "---ACCC--")
        read1.setSignificantOffsets([3, 4])
        read2 = AlignedRead("id2", "---TG--")
        read2.setSignificantOffsets([3])
        rc1 = ReadCluster()
        rc1.add(read1)
        rc1.add(read2)

        read3 = AlignedRead("id3", "---GGCC--")
        read3.setSignificantOffsets([3, 4, 6])
        rc2 = ReadCluster()
        rc2.add(read3)

        expectedBases3 = OffsetBases()
        expectedBases3.incorporateBase("A")
        expectedBases3.incorporateBase("T")
        expectedBases3.incorporateBase("G")

        expectedBases4 = OffsetBases()
        expectedBases4.incorporateBase("C")
        expectedBases4.incorporateBase("G")

        expectedBases6 = OffsetBases()
        expectedBases6.incorporateBase("C")

        rc1.merge(rc2)

        self.assertEqual(
            {
                3: expectedBases3,
                4: expectedBases4,
                6: expectedBases6,
            },
            rc1.nucleotides,
        )

    def testDistanceOne(self):
        """
        The distance between two clusters must be zero when the nucleotides
        they have at their signifcant offsets do not match.
        """
        read1 = AlignedRead("id1", "---ACCCG-")
        read1.setSignificantOffsets([3, 4, 5, 6, 7])
        rc1 = ReadCluster()
        rc1.add(read1)

        read2 = AlignedRead("id2", "---TG--A")
        read2.setSignificantOffsets([3, 4, 7])
        rc2 = ReadCluster()
        rc2.add(read2)

        self.assertEqual(1.0, ReadCluster.commonNucleotidesAgreementDistance(rc1, rc2))

    def testDistanceZero(self):
        """
        The distance between two clusters must be 0.0 when they match at all
        common significant offsets.
        """
        read1 = AlignedRead("id1", "---TCCCG-")
        read1.setSignificantOffsets([3, 4, 5, 6, 7])
        rc1 = ReadCluster()
        rc1.add(read1)

        read2 = AlignedRead("id2", "---TC--A")
        read2.setSignificantOffsets([3, 4, 7])
        read3 = AlignedRead("id3", "---TG--G")
        read3.setSignificantOffsets([7])
        rc2 = ReadCluster()
        rc2.add(read2)
        rc2.add(read3)

        # All three common offsets match.
        self.assertEqual(0.0, ReadCluster.commonNucleotidesAgreementDistance(rc1, rc2))

    def testNonZeroDistance(self):
        """
        The distance between two clusters must be as expected.
        """
        read1 = AlignedRead("id1", "---TCCCG-")
        read1.setSignificantOffsets([3, 4, 5, 6, 7])
        rc1 = ReadCluster()
        rc1.add(read1)

        read2 = AlignedRead("id2", "---TG--A")
        read2.setSignificantOffsets([3, 4, 7])
        rc2 = ReadCluster()
        rc2.add(read2)

        # Two of the three common offsets do not match.
        self.assertAlmostEqual(
            2 / 3, ReadCluster.commonNucleotidesAgreementDistance(rc1, rc2)
        )

    def testCommonOffsetsMaxFractionOne(self):
        """
        The commonOffsetsMaxFraction method must return 1.0 when all of
        one of the offsets in one cluster are in the intersection of the
        offsets of another cluster.
        """
        read1 = AlignedRead("id1", "---TCCCG-")
        read1.setSignificantOffsets([3, 4])
        read2 = AlignedRead("id2", "---TCCCG-")
        read2.setSignificantOffsets([6, 7])
        rc1 = ReadCluster()
        rc1.add(read1)
        rc1.add(read2)

        read3 = AlignedRead("id3", "---TG--A")
        read3.setSignificantOffsets([3, 4, 7])
        rc2 = ReadCluster()
        rc2.add(read3)

        self.assertAlmostEqual(1.0, ReadCluster.commonOffsetsMaxFraction(rc1, rc2))

    def testCommonOffsetsMaxFractionZero(self):
        """
        The commonOffsetsMaxFraction method must return 0.0 when two
        clusters have no offsets in common.
        """
        read1 = AlignedRead("id1", "---TCCCG-")
        read1.setSignificantOffsets([5, 6])
        rc1 = ReadCluster()
        rc1.add(read1)

        read2 = AlignedRead("id2", "---TG--A")
        read2.setSignificantOffsets([3, 4, 7])
        rc2 = ReadCluster()
        rc2.add(read2)

        self.assertAlmostEqual(0.0, ReadCluster.commonOffsetsMaxFraction(rc1, rc2))

    def testCommonOffsetsMaxFractionOneHalf(self):
        """
        The commonOffsetsMaxFraction method must return 0.5 when half
        the offsets of one cluster are in common with the offsets of another
        cluster.
        """
        read1 = AlignedRead("id1", "---TCCCG-")
        read1.setSignificantOffsets([5, 6])
        rc1 = ReadCluster()
        rc1.add(read1)

        read2 = AlignedRead("id2", "---TG--A")
        read2.setSignificantOffsets([3, 4, 6, 7])
        rc2 = ReadCluster()
        rc2.add(read2)

        self.assertAlmostEqual(0.0, ReadCluster.commonOffsetsMaxFraction(rc1, rc2))

    def testMultiplicativeDistanceOne(self):
        """
        The commonNucleotidesMultiplicativeDistance method must return 1.0
        when two clusters have no offsets in common.
        """
        read1 = AlignedRead("id1", "-----TC-")
        read1.setSignificantOffsets([5, 6])
        rc1 = ReadCluster()
        rc1.add(read1)

        read2 = AlignedRead("id2", "---TG--")
        read2.setSignificantOffsets([3, 4])
        rc2 = ReadCluster()
        rc2.add(read2)

        self.assertAlmostEqual(
            1.0, ReadCluster.commonNucleotidesMultiplicativeDistance(rc1, rc2)
        )

    def testMultiplicativeDistanceZero(self):
        """
        The commonNucleotidesMultiplicativeDistance method must return 0.0
        when two clusters agree completely on nucleotides at their common
        offsets.
        """
        read1 = AlignedRead("id1", "-----TC-")
        read1.setSignificantOffsets([5, 6])
        rc1 = ReadCluster()
        rc1.add(read1)

        read2 = AlignedRead("id2", "-----TC-")
        read2.setSignificantOffsets([5, 6])
        read3 = AlignedRead("id3", "-----TC-")
        read3.setSignificantOffsets([5, 6])
        rc2 = ReadCluster()
        rc2.add(read2)
        rc2.add(read3)

        self.assertAlmostEqual(
            0.0, ReadCluster.commonNucleotidesMultiplicativeDistance(rc1, rc2)
        )

    def testMultiplicativeDistanceOneQuarter(self):
        """
        The commonNucleotidesMultiplicativeDistance method must return 0.25
        when the ratio of two bases in one cluster is 50:50. You need to look
        at that method to see what's going on. It's more complicated than it
        used to be because it uses a min of a multiplicative distance and a
        homogeneous one.
        """
        read1 = AlignedRead("id1", "---TCTC-")
        read1.setSignificantOffsets([3, 4, 5, 6])
        rc1 = ReadCluster()
        rc1.add(read1)

        read2 = AlignedRead("id2", "---TG--")
        read2.setSignificantOffsets([3, 4])
        rc2 = ReadCluster()
        rc2.add(read2)

        self.assertAlmostEqual(
            0.25, ReadCluster.commonNucleotidesMultiplicativeDistance(rc1, rc2)
        )


class TestReadClusters(TestCase):
    """
    Test the ReadClusters class.
    """

    def testNonZeroCommonNucleotidesAgreementDistance(self):
        """
        The distance between two clusters must be as expected, as calculated by
        the commonNucleotidesAgreementDistance method.
        """
        read1 = AlignedRead("id1", "---TCCCG-")
        read1.setSignificantOffsets([3, 4, 5, 6, 7])

        read2 = AlignedRead("id2", "---TG--A")
        read2.setSignificantOffsets([3, 4, 7])

        rc = ReadClusters()
        cluster1 = rc.add(read1)
        cluster2 = rc.add(read2)

        # Two of the three common offsets do not match.
        self.assertAlmostEqual(
            2 / 3, rc.commonNucleotidesAgreementDistance(cluster1, cluster2)
        )

    def testMultiplicativeDistanceOne(self):
        """
        The multiplicative distance between two clusters must be 1.0 when
        they have no offsets in common.
        """
        read1 = AlignedRead("id1", "-----CCG-")
        read1.setSignificantOffsets([5, 6, 7])

        read2 = AlignedRead("id2", "---TG--")
        read2.setSignificantOffsets([3, 4])

        rc = ReadClusters()
        cluster1 = rc.add(read1)
        cluster2 = rc.add(read2)

        self.assertAlmostEqual(1.0, rc.multiplicativeDistance(cluster1, cluster2))

    def testMultiplicativeDistanceOneTenth(self):
        """
        The multiplicative distance between two clusters must be 0.1 when
        they are identical and the minimum fraction of common offsets
        (ReadClusters.COMMON_OFFSETS_MAX_FRACTION_MIN) is 0.9. That fraction
        is used in the following because 3 of the first cluster's sites and
        3 of the second's are in common, and those fractions are 3/4 and 3/6
        which are both less than the 0.9 value of
        ReadClusters.COMMON_OFFSETS_MAX_FRACTION_MIN so it it used to scale
        the commonNucleotidesAgreementDistance distance.
        """
        # Check that the above comment is possibly correct :-)
        self.assertEqual(ReadClusters.COMMON_OFFSETS_MAX_FRACTION_MIN, 0.9)

        read1 = AlignedRead("id1", "-----CCGT")
        read1.setSignificantOffsets([5, 6, 7, 8])

        read2 = AlignedRead("id2", "--TTTCCG-")
        read2.setSignificantOffsets([2, 3, 4, 5, 6, 7])

        rc = ReadClusters()
        cluster1 = rc.add(read1)
        cluster2 = rc.add(read2)

        self.assertAlmostEqual(0.1, rc.multiplicativeDistance(cluster1, cluster2))

    def testMultiplicativeDistanceOneTenthLowOffsetCoverage(self):
        """
        The multiplicative distance between two clusters must be 0.1 when
        they are identical but the maximum fraction of common offsets is
        one half (in which case the 0.9 minimum offset coverage fraction
        will be applied, as explained in the docstring for the
        testMultiplicativeDistanceOneQuarter test above).
        """
        read1 = AlignedRead("id1", "-----CCGTTT")
        read1.setSignificantOffsets([5, 6, 7, 8, 9, 10])

        read2 = AlignedRead("id2", "--TTTCCG-")
        read2.setSignificantOffsets([2, 3, 4, 5, 6, 7])

        rc = ReadClusters()
        cluster1 = rc.add(read1)
        cluster2 = rc.add(read2)

        self.assertAlmostEqual(0.1, rc.multiplicativeDistance(cluster1, cluster2))

    def testMultiplicativeDistanceTwoOffsetsInCommonOneAgreesOneDoesNot(self):
        """
        The explanation for the numeric result here can only be seen by
        """
        read1 = AlignedRead("id1", "-----CAGT")
        read1.setSignificantOffsets([5, 6, 7, 8])

        read2 = AlignedRead("id2", "--TTTCC--")
        read2.setSignificantOffsets([2, 3, 4, 5, 6])

        rc = ReadClusters()
        cluster1 = rc.add(read1)
        cluster2 = rc.add(read2)

        self.assertAlmostEqual(
            1.0 - ((1.0 - 0.25) * 0.9), rc.multiplicativeDistance(cluster1, cluster2)
        )

    def testMultiplicativeDistanceZero(self):
        """
        The multiplicative distance between two clusters must be 0.0 when
        they are identical.
        """
        read1 = AlignedRead("id1", "-----CCG-")
        read1.setSignificantOffsets([5, 6, 7])

        read2 = AlignedRead("id2", "-----CCG-")
        read2.setSignificantOffsets([5, 6, 7])

        rc = ReadClusters()
        cluster1 = rc.add(read1)
        cluster2 = rc.add(read2)

        self.assertAlmostEqual(0.0, rc.multiplicativeDistance(cluster1, cluster2))

    def testFullMerge(self):
        """
        A merge of three reads must result in a cluster with the expected
        nucleotide counts.
        """
        read1 = AlignedRead("id1", "---TCCCG-")
        read1.setSignificantOffsets([3, 4, 5, 6, 7])

        read2 = AlignedRead("id2", "---AC--A")
        read2.setSignificantOffsets([3, 4, 7])

        read3 = AlignedRead("id3", "---TG--G")
        read3.setSignificantOffsets([7])

        rc = ReadClusters()
        rc.add(read1)
        rc.add(read2)
        rc.add(read3)

        (cluster,) = rc.analyze(0.7)

        expectedBases3 = OffsetBases()
        expectedBases3.incorporateBase("T")
        expectedBases3.incorporateBase("A")

        expectedBases4 = OffsetBases()
        expectedBases4.incorporateBase("C")
        expectedBases4.incorporateBase("C")

        expectedBases5 = OffsetBases()
        expectedBases5.incorporateBase("C")

        expectedBases6 = OffsetBases()
        expectedBases6.incorporateBase("C")

        expectedBases7 = OffsetBases()
        expectedBases7.incorporateBase("G")
        expectedBases7.incorporateBase("A")
        expectedBases7.incorporateBase("G")

        self.assertEqual(
            {
                3: expectedBases3,
                4: expectedBases4,
                5: expectedBases5,
                6: expectedBases6,
                7: expectedBases7,
            },
            cluster.nucleotides,
        )
