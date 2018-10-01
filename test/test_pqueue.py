from unittest import TestCase

from six import assertRaisesRegex

from midtools.pqueue import PriorityQueue


class TestPriorityQueue(TestCase):
    """
    Test the PriorityQueue class.
    """
    def testLowestPriorityRaisesOnEmpty(self):
        """
        The lowestPriority method must raise a KeyError if the queue is empty.
        """
        pq = PriorityQueue()
        error = 'peek on an empty priority queue'
        assertRaisesRegex(self, KeyError, error, pq.lowestPriority)

    def testPopRaisesOnEmpty(self):
        """
        The pop method must raise a KeyError if the queue is empty.
        """
        pq = PriorityQueue()
        error = 'pop from an empty priority queue'
        assertRaisesRegex(self, KeyError, error, pq.pop)

    def testEmptyLength(self):
        """
        An empty queue must have a zero length.
        """
        pq = PriorityQueue()
        self.assertEqual(0, len(pq))

    def testLengthAfterAdd(self):
        """
        An queue with one thing in it must have a length of one.
        """
        pq = PriorityQueue()
        pq.add(3)
        self.assertEqual(1, len(pq))

    def testLengthAfterDoubleAdd(self):
        """
        An queue with an item that is added twice must have a length of one.
        """
        pq = PriorityQueue()
        pq.add(3)
        pq.add(3)
        self.assertEqual(1, len(pq))

    def testContains(self):
        """
        The __contains__ function must work as expected.
        """
        pq = PriorityQueue()
        pq.add('hey')
        self.assertTrue('hey' in pq)
        self.assertFalse('hi' in pq)

    def testBooleanWhenEmpty(self):
        """
        The queue must test False when empty.
        """
        pq = PriorityQueue()
        self.assertFalse(pq)

    def testBooleanWhenNotEmpty(self):
        """
        The queue must test False when empty.
        """
        pq = PriorityQueue()
        pq.add(4)
        self.assertTrue(pq)

    def testPop(self):
        """
        The queue pop method must return the lowest priority item.
        """
        pq = PriorityQueue()
        pq.add('hey', 9)
        pq.add('you', 7)
        pq.add('two', 8)
        self.assertEqual('you', pq.pop())

    def testRemove(self):
        """
        The remove method must work as expected.
        """
        pq = PriorityQueue()
        pq.add('hey')
        pq.remove('hey')
        self.assertFalse(pq)
        self.assertEqual(0, len(pq))

    def testLowestPriority(self):
        """
        The lowestPriority method must return the lowest priority without
        disturbing the length of the queue.
        """
        pq = PriorityQueue()
        pq.add('hey', 9)
        pq.add('you', 7)
        pq.add('two', 8)
        self.assertEqual(7, pq.lowestPriority())
        self.assertEqual(3, len(pq))
