from unittest import TestCase

from dark.reads import Read

from data.data import AlignedRead


class TestAlignedRead(TestCase):
    """
    Tests for the AlignedRead class.
    """
    def testTrim(self):
        """
        The trim function must work as expected.
        """
        ar = AlignedRead(Read('id', '---ACGTACGT--'))
        self.assertEqual(8, len(ar.read))
        self.assertTrue(ar.trim(2))
        self.assertEqual('GTAC', ar.read.sequence)

    def testTrimZero(self):
        """
        The trim function must work as expected when the trim quantity is 0.
        """
        ar = AlignedRead(Read('id', '---ACGTACGT--'))
        self.assertTrue(ar.trim(0))
        self.assertEqual('ACGTACGT', ar.read.sequence)

    def testTrimWithNegativeAmount(self):
        """
        The trim function must raise an AssertionError if the amount to trim
        is negative.
        """
        ar = AlignedRead(Read('id', '---ACGTACGT--'))
        error = '^Trim amount \(-4\) cannot be negative\.$'
        self.assertRaisesRegex(AssertionError, error, ar.trim, -4)

    def testTrimWithReadTooShort(self):
        """
        The trim function must return False if the read is too short.
        """
        ar = AlignedRead(Read('id', '---ACGTACGT--'))
        self.assertFalse(ar.trim(4))
