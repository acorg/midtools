from unittest import TestCase

from dark.reads import Read

from midtools.mutate import mutateRead


class TestMutateRead(TestCase):
    """
    Tests for the mutate read function.
    """

    def testRateOne(self):
        """
        If the mutation rate is 1.0 all bases must be mutated.
        """
        read = Read("id", "ACGTACGT")
        offsets = mutateRead(read, 1.0, "Z")
        self.assertEqual("ZZZZZZZZ", read.sequence)
        self.assertEqual(8, len(offsets))

    def testRateZero(self):
        """
        If the mutation rate is 0.0 no bases can be mutated.
        """
        read = Read("id", "ACGTACGT")
        offsets = mutateRead(read, 0.0, "Z")
        self.assertEqual("ACGTACGT", read.sequence)
        self.assertEqual(0, len(offsets))
