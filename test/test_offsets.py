from unittest import TestCase

from midtools.offsets import OffsetBases

multiplicativeDistance = OffsetBases.multiplicativeDistance
homogeneousDistance = OffsetBases.homogeneousDistance
highestFrequenciesMultiple = OffsetBases.highestFrequenciesMultiple


class TestOffsetBases(TestCase):
    """
    Test the OffsetBases class.
    """

    def testInitialized(self):
        """
        The class must have the expected content when initialized.
        """
        ob = OffsetBases()
        ob.incorporateBase("a")
        self.assertEqual({"a"}, ob.commonest)

    def testAddNewInstanceOfCommonest(self):
        """
        Adding another copy of a base that is already the commonest must
        leave the commonest set unchanged.
        """
        ob = OffsetBases()
        ob.incorporateBase("a")
        ob.incorporateBase("a")
        self.assertEqual({"a"}, ob.commonest)

    def testDraw(self):
        """
        Adding a base that creates a draw must result in the commonest set
        being added to.
        """
        ob = OffsetBases()
        ob.incorporateBase("a")
        ob.incorporateBase("g")
        self.assertEqual({"a", "g"}, ob.commonest)

    def testCreateDrawThenResolveIt(self):
        """
        Adding a base that creates a draw and then bases that remove the draw
        must result in the commonest set being as expected.
        """
        ob = OffsetBases()
        ob.incorporateBase("a")
        ob.incorporateBase("g")
        ob.incorporateBase("t")
        ob.incorporateBase("t")
        self.assertEqual({"t"}, ob.commonest)

    def testMerge(self):
        """
        Merging another instance must result in the commonest set being
        as expected.
        """
        ob1 = OffsetBases()
        ob1.incorporateBase("a")
        ob1.incorporateBase("g")
        ob1.incorporateBase("t")
        ob1.incorporateBase("t")

        ob2 = OffsetBases()
        ob2.incorporateBase("g")
        ob2.incorporateBase("c")
        ob2.incorporateBase("c")

        ob1.merge(ob2)

        self.assertEqual({"c", "g", "t"}, ob1.commonest)


class TestMultiplicativeDistance(TestCase):
    """
    Test the multiplicativeDistance function.
    """

    def testZero(self):
        """
        The multiplicativeDistance method must return zero when the bases at an
        offset completely agree.
        """
        ob1 = OffsetBases()
        ob1.incorporateBase("a")
        ob1.incorporateBase("a")

        ob2 = OffsetBases()
        ob2.incorporateBase("a")

        self.assertEqual(0.0, multiplicativeDistance(ob1, ob2))

    def testOne(self):
        """
        The multiplicativeDistance method must return one when the bases at an
        offset completely disagree.
        """
        ob1 = OffsetBases()
        ob1.incorporateBase("t")
        ob1.incorporateBase("c")

        ob2 = OffsetBases()
        ob2.incorporateBase("a")

        self.assertEqual(1.0, multiplicativeDistance(ob1, ob2))

    def testOneHalf(self):
        """
        The multiplicativeDistance method must return 0.5 when the two bases
        at an offset are balanced (because (0.5 * 0.5) + (0.5 * 0.5) = 0.25).
        """
        ob1 = OffsetBases()
        ob1.incorporateBase("t")
        ob1.incorporateBase("a")

        ob2 = OffsetBases()
        ob2.incorporateBase("t")
        ob2.incorporateBase("a")

        self.assertEqual(0.5, multiplicativeDistance(ob1, ob2))

    def testSixSixteenths(self):
        """
        The multiplicativeDistance method must return 10/16 when both
        have AATC because (0.5 * 0.5) + (0.25 * 0.25) + (0.25 * 0.25) =
        4/16 + 1/16 + 1/16 = 6/16 (and 1 - 6/16 = 10/16).
        """
        ob1 = OffsetBases()
        ob1.incorporateBase("a")
        ob1.incorporateBase("a")
        ob1.incorporateBase("c")
        ob1.incorporateBase("t")

        ob2 = OffsetBases()
        ob2.incorporateBase("a")
        ob2.incorporateBase("a")
        ob2.incorporateBase("c")
        ob2.incorporateBase("t")

        self.assertEqual(10 / 16, multiplicativeDistance(ob1, ob2))

    def testFourNinths(self):
        """
        The multiplicativeDistance method must return 4/9 when both
        have AAT because (2/3 * 2/3) + (1/3 * 1/3) = 5/9 (and 1 - 5/9 = 4/9).
        """
        ob1 = OffsetBases()
        ob1.incorporateBase("a")
        ob1.incorporateBase("a")
        ob1.incorporateBase("t")

        ob2 = OffsetBases()
        ob2.incorporateBase("a")
        ob2.incorporateBase("a")
        ob2.incorporateBase("t")

        self.assertAlmostEqual(4 / 9, multiplicativeDistance(ob1, ob2))


class TestHomogeneousDistance(TestCase):
    """
    Test the homogeneousDistance function.
    """

    def testZero(self):
        """
        The homogeneousDistance method must return zero when both there is no
        discrepancy in the base.
        """
        ob1 = OffsetBases()
        ob1.incorporateBase("a")
        ob1.incorporateBase("a")

        ob2 = OffsetBases()
        ob2.incorporateBase("a")

        self.assertAlmostEqual(0.0, homogeneousDistance(ob1, ob2))

    def testOneHalf(self):
        """
        The homogeneousDistance method must return 0.5 when the nucleotides
        are equally divided between two choices.
        """
        ob1 = OffsetBases()
        ob1.incorporateBase("a")
        ob1.incorporateBase("a")
        ob1.incorporateBase("g")

        ob2 = OffsetBases()
        ob2.incorporateBase("g")

        self.assertAlmostEqual(0.5, homogeneousDistance(ob1, ob2))

    def testPointOne(self):
        """
        The homogeneousDistance method must return 0.1 when the maximum
        nucleotide fraction is 0.9.
        """
        ob1 = OffsetBases()
        for _ in range(40):
            ob1.incorporateBase("a")
        for _ in range(3):
            ob1.incorporateBase("g")

        ob2 = OffsetBases()
        for _ in range(50):
            ob2.incorporateBase("a")
        for _ in range(7):
            ob2.incorporateBase("g")

        self.assertAlmostEqual(0.1, homogeneousDistance(ob1, ob2))

    def testThreeQuarters(self):
        """
        The homogeneousDistance method must return 0.75 when all nucleotides
        are equally represented. This is its maximum distance.
        """
        ob1 = OffsetBases()
        ob1.incorporateBase("a")
        ob1.incorporateBase("c")
        ob1.incorporateBase("g")

        ob2 = OffsetBases()
        ob2.incorporateBase("t")

        self.assertAlmostEqual(0.75, homogeneousDistance(ob1, ob2))


class TestHighestFrequenciesMultiple(TestCase):
    """
    Test the highestFrequenciesMultiple function.
    """

    def testOneNucleotide(self):
        """
        When two offsets both have just a single identical nucleotide, the
        highestFrequenciesMultiple method must return None.
        """
        ob1 = OffsetBases()
        ob1.incorporateBase("a")

        ob2 = OffsetBases()
        ob2.incorporateBase("a")

        self.assertIs(None, highestFrequenciesMultiple(ob1, ob2))

    def testTwo(self):
        """
        When the frequency of the most common nucleotide in the sum of two
        OffsetBases instance is twice the second most common, the
        highestFrequenciesMultiple method must return 2.0.
        """
        ob1 = OffsetBases()
        ob1.incorporateBase("a")
        ob1.incorporateBase("a")
        ob1.incorporateBase("a")
        ob1.incorporateBase("a")
        ob1.incorporateBase("a")
        ob1.incorporateBase("a")
        ob1.incorporateBase("c")

        ob2 = OffsetBases()
        ob2.incorporateBase("a")
        ob2.incorporateBase("a")
        ob2.incorporateBase("a")
        ob2.incorporateBase("a")
        ob2.incorporateBase("c")
        ob2.incorporateBase("c")
        ob2.incorporateBase("c")
        ob2.incorporateBase("c")
        ob2.incorporateBase("g")
        ob2.incorporateBase("t")

        # The total count has 10 x 'a' and 5 x 'c'.
        self.assertEqual(2.0, highestFrequenciesMultiple(ob1, ob2))

    def testDraw(self):
        """
        When the frequency of the two most common nucleotides in the sum of two
        OffsetBases instance is the same, the highestFrequenciesMultiple method
        must return 1.0.
        """
        ob1 = OffsetBases()
        ob1.incorporateBase("a")
        ob1.incorporateBase("a")
        ob1.incorporateBase("a")
        ob1.incorporateBase("c")

        ob2 = OffsetBases()
        ob2.incorporateBase("a")
        ob2.incorporateBase("c")
        ob2.incorporateBase("c")
        ob2.incorporateBase("c")
        ob2.incorporateBase("g")
        ob2.incorporateBase("t")

        # The total count has 4 x 'a' and 4 x 'c'.
        self.assertEqual(1.0, highestFrequenciesMultiple(ob1, ob2))
