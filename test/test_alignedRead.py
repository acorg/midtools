from unittest import TestCase

from midtools.read import AlignedRead


class TestAlignedRead(TestCase):
    """
    Tests for the AlignedRead class.
    """

    def testTrim(self):
        """
        The trim function must work as expected.
        """
        ar = AlignedRead("id", "---ACGTACGT--")
        self.assertEqual(8, len(ar))
        self.assertTrue(ar.trim(2))
        self.assertEqual("GTAC", ar.sequence)

    def testTrimZero(self):
        """
        The trim function must work as expected when the trim quantity is 0.
        """
        ar = AlignedRead("id", "---ACGTACGT--")
        self.assertTrue(ar.trim(0))
        self.assertEqual("ACGTACGT", ar.sequence)

    def testTrimWithNegativeAmount(self):
        """
        The trim function must raise an AssertionError if the amount to trim
        is negative.
        """
        ar = AlignedRead("id", "---ACGTACGT--")
        error = r"^Trim amount \(-4\) cannot be negative\.$"
        self.assertRaisesRegex(AssertionError, error, ar.trim, -4)

    def testTrimWithReadTooShort(self):
        """
        The trim function must return False if the read is too short.
        """
        ar = AlignedRead("id", "---ACGTACGT--")
        self.assertFalse(ar.trim(4))

    def testSetSignificantOffsets(self):
        """
        Adding some significant offsets must work as expected.
        """
        ar = AlignedRead("id", "---ACGTACGT--")
        ar.setSignificantOffsets([2, 3, 4])
        self.assertEqual({3: "A", 4: "C"}, ar.significantOffsets)

    def testBaseNone(self):
        """
        The base method must return C{None} for an offset the read does not
        have.
        """
        ar = AlignedRead("id", "---ACGTACGT--")
        self.assertEqual(None, ar.base(0))

    def testBase(self):
        """
        The base method must return the expected nucleotide.
        """
        ar = AlignedRead("id", "---ACGTGCGT--")
        self.assertEqual("A", ar.base(3))

    def testToPaddedString(self):
        """
        The toPaddedString method must return the expected string.
        """
        ar = AlignedRead("id", "---ACGTGCGT--")
        self.assertEqual(">id\n---ACGTGCGT--\n", ar.toPaddedString())
