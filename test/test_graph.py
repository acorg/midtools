from unittest import TestCase

from dark.reads import Read
from data.data import AlignedRead
from data.graph import Node, connectedComponents, componentOffsets


class TestConnectedComponents(TestCase):
    """
    Tests for the connectedComponents function.
    """
    def testEmpty(self):
        """
        If an empty set of nodes is passed to connectedComponents it must
        return an empty generator.
        """
        self.assertEqual([], list(connectedComponents([])))

    def testOneNode(self):
        """
        A set of one node should result in one connected component with just
        that one node in it.
        """
        read = AlignedRead('id', 'ACGT')
        node = Node(read)
        self.assertEqual(
            [
                {node},
            ],
            list(connectedComponents([node])))

    def testTwoUnconnectedNodes(self):
        """
        A list of two unconnected nodes should result in two connected
        components with one node in each.
        """
        read1 = AlignedRead('id1', 'ACGT')
        node1 = Node(read1)
        read2 = AlignedRead('id2', 'ACGT')
        node2 = Node(read2)
        result = list(connectedComponents([node1, node2]))
        self.assertTrue(
            result == [{node1}, {node2}] or result == [{node2}, {node1}])

    def testTwoConnectedNodes(self):
        """
        A list of two connected nodes should result in one connected
        component containing both nodes.
        """
        read1 = AlignedRead('id1', 'ACGT')
        node1 = Node(read1)
        read2 = AlignedRead('id2', 'ACGT')
        node2 = Node(read2)
        node1.add(node2)
        node2.add(node1)
        result = list(connectedComponents([node1, node2]))
        self.assertEqual(result, [{node1, node2}])

    def testTwoConnectedNodesOneUnconnected(self):
        """
        A list of two connected nodes and one that's unconnected should result
        in one connected component containing both nodes and one singleton.
        """
        read1 = AlignedRead('id1', 'ACGT')
        node1 = Node(read1)
        read2 = AlignedRead('id2', 'ACGT')
        node2 = Node(read2)
        node1.add(node2)
        node2.add(node1)
        read3 = AlignedRead('id3', 'ACGT')
        node3 = Node(read3)
        result = list(connectedComponents([node1, node2, node3]))
        self.assertTrue(
            result == [{node1, node2}, {node3}] or
            result == [{node3}, {node1, node2}])

    def testRepeatTwoConnectedNodesOneUnconnected(self):
        """
        A list of two connected nodes and one that's unconnected should result
        in one connected component containing both nodes and one singleton.
        This must be able to be repeated (i.e., the code in
        connectedComponents must not disrupt the original graph data).
        """
        read1 = AlignedRead('id1', 'ACGT')
        node1 = Node(read1)
        read2 = AlignedRead('id2', 'ACGT')
        node2 = Node(read2)
        node1.add(node2)
        node2.add(node1)
        read3 = AlignedRead('id3', 'ACGT')
        node3 = Node(read3)

        for _ in range(5):
            result = list(connectedComponents([node1, node2, node3]))
            self.assertTrue(
                result == [{node1, node2}, {node3}] or
                result == [{node3}, {node1, node2}])


class TestComponentOffsets(TestCase):
    """
    Tests for the componentOffsets function.
    """
    def testNoGaps(self):
        """
        A read with no leading gaps should have offsets starting from zero.
        """
        read = AlignedRead('id', 'ACGT')
        node = Node(read)
        self.assertEqual({0, 1, 2, 3}, componentOffsets({node}))

    def testLeadingGaps(self):
        """
        A read with 3 leading gaps should have offsets starting from 3.
        """
        read = AlignedRead('id', '---ACGT')
        node = Node(read)
        self.assertEqual({3, 4, 5, 6}, componentOffsets({node}))

    def testTrailingGaps(self):
        """
        A read with trailing gaps should have offsets starting from 0.
        """
        read = AlignedRead('id', 'ACGT----')
        node = Node(read)
        self.assertEqual({0, 1, 2, 3}, componentOffsets({node}))

    def testLeadingAndTrailingGaps(self):
        """
        A read with 2 leading gaps and some trailing gaps should have offsets
        starting from 2.
        """
        read = AlignedRead('id', '--ACGT-------')
        node = Node(read)
        self.assertEqual({2, 3, 4, 5}, componentOffsets({node}))
