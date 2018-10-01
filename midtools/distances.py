from collections import defaultdict

from midtools.pqueue import PriorityQueue


def _key(a, b):
    return (a, b) if a <= b else (b, a)


class DistanceCache(object):
    """
    Maintain a set of distances between objects, with lazy evaluation and
    removal from the set.

    @param distFunc: A function that computes the distance between two objects.
    """
    def __init__(self, distFunc):
        self._distFunc = distFunc
        self._distances = defaultdict(dict)
        self._pq = PriorityQueue()

    def distance(self, a, b):
        """
        Find the distance between a pair of objects.

        @param a: An immutable object.
        @param b: An immutable object.
        @return: The distance between C{a} and C{b}, according to the distance
            function passed to __init__.
        """
        return self._distances[a][b]

    def add(self, a):
        """
        Add an object.

        @param a: An immutable object.
        """
        assert a not in self._distances

        if self._distances:
            for b in list(self._distances):
                d = self._distFunc(a, b)
                self._distances[b][a] = self._distances[a][b] = d
                self._pq.add(_key(a, b), d)
        else:
            # This is the first element, so it has no distances to
            # anything.  Mention it to create its distance dictionary so it
            # will be found when subsequent elements are added.
            self._distances[a]

    def lowestDistance(self):
        """
        Get the lowest distance between any two clusters.

        @return: A C{float} distance.
        """
        try:
            return self._pq.lowestPriority()
        except KeyError:
            return None

    def pop(self):
        """
        Pop the lowest distance cluster pair.

        @raise KeyError: If the distance priority queue is empty.
        @return: A 2-C{tuple} of C{int} cluster numbers.
        """
        return self._pq.pop()

    def __contains__(self, pair):
        """
        Test if a pair has a computed distance (useful for testing).

        @param pair: A 2-tuple of objects.
        @return: A Boolean indicating membership.
        """
        return ((pair[0], pair[1]) in self._distances or
                (pair[1], pair[0]) in self._distances)

    def remove(self, a):
        """
        Remove an object.

        @param a: An object.
        """
        errorCount = 0
        for b in self._distances:
            if b != a:
                try:
                    self._pq.remove(_key(a, b))
                except KeyError:
                    # We allow one KeyError since 'a' has likely just been
                    # popped as part of the lowest scoring pain.
                    errorCount += 1
                    if errorCount > 1:
                        raise
                del self._distances[b][a]

        del self._distances[a]
