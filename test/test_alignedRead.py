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
        ar = AlignedRead('id', '---ACGTACGT--')
        self.assertEqual(8, len(ar))
        self.assertTrue(ar.trim(2))
        self.assertEqual('GTAC', ar.sequence)

    def testTrimZero(self):
        """
        The trim function must work as expected when the trim quantity is 0.
        """
        ar = AlignedRead('id', '---ACGTACGT--')
        self.assertTrue(ar.trim(0))
        self.assertEqual('ACGTACGT', ar.sequence)

    def testTrimWithNegativeAmount(self):
        """
        The trim function must raise an AssertionError if the amount to trim
        is negative.
        """
        ar = AlignedRead('id', '---ACGTACGT--')
        error = '^Trim amount \(-4\) cannot be negative\.$'
        self.assertRaisesRegex(AssertionError, error, ar.trim, -4)

    def testTrimWithReadTooShort(self):
        """
        The trim function must return False if the read is too short.
        """
        ar = AlignedRead('id', '---ACGTACGT--')
        self.assertFalse(ar.trim(4))
