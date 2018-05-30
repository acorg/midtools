class Node(set):
    """
    Hold an aligned read and the set of other aligned reads it is connected to
    in a graph.

    @param read: An C{AlignedRead} instance.
    """
    def __init__(self, read):
        self.read = read
        set.__init__(self)

    def __str__(self):
        return '<Node %s Connected to %d>' % (self.read.__str__(), len(self))

    def __hash__(self):
        return self.read.__hash__()

    def __lt__(self, other):
        return (self.read.offset, self.read.read.id) < (
            other.read.offset, other.read.read.id)

    def neighbors(self):
        return set(self)


def connectedComponents(nodes):
    """
    Find connected components in a set of nodes.
    """
    nodes = set(nodes)

    while nodes:
        n = nodes.pop()

        # The connected component we are building.
        component = {n}

        # Maintain a queue of things yet to be processed for this component.
        componentQueue = [n]

        # Iterate the queue. When it is empty, we have finished visiting a
        # group of connected nodes.
        while componentQueue:

            # Consume the next item from the component queue.
            n = componentQueue.pop(0)
            neighbors = n.neighbors()

            # Remove the neighbors we already visited.
            neighbors.difference_update(component)

            # Remove the remaining nodes from the global set.
            nodes.difference_update(neighbors)

            # Add them to the group of connected nodes.
            component.update(neighbors)

            # Add them to the component queue, so we visit them in the next
            # iterations.
            componentQueue.extend(neighbors)

        yield component


def componentOffsets(component):
    """
    Get all genome offsets spanned by the reads in a connected component.
    """
    offsets = set()
    for node in component:
        read = node.read
        offsets.update(range(read.offset, read.offset + len(read.read)))
    return offsets
