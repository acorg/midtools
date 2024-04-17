from typing import Any, Callable, Optional
from collections import defaultdict

from midtools.pqueue import PriorityQueue


def _key(a, b) -> tuple[int, int]:
    return (a, b) if a <= b else (b, a)


class DistanceCache:
    """
    Maintain a set of distances between objects, with lazy evaluation and
    removal from the set.

    @param distFunc: A function that computes the distance between two objects.
    """

    def __init__(self, distFunc: Callable[[Any, Any], float]) -> None:
        self._distFunc = distFunc
        self._distances: dict[Any, dict[Any, float]] = defaultdict(dict)
        self._pq = PriorityQueue()

    def distance(self, a: Any, b: Any) -> float:
        """
        Find the distance between a pair of objects.

        @param a: An immutable object.
        @param b: An immutable object.
        @return: The distance between C{a} and C{b}, according to the distance
            function passed to __init__.
        """
        return self._distances[a][b]

    def add(self, a: Any) -> None:
        """
        Add an object.

        @param a: An immutable object.
        """
        assert a not in self._distances

        if self._distances:
            for b in list(self._distances):
                distance = self._distFunc(a, b)
                self._distances[b][a] = self._distances[a][b] = distance
                self._pq.add(_key(a, b), distance)
        else:
            # This is the first element, so it has no distances to
            # anything.  Mention it to create its distance dictionary so it
            # will be found when subsequent elements are added.
            self._distances[a]

    def lowestDistance(self) -> Optional[float]:
        """
        Get the lowest distance between any two clusters.

        @return: A C{float} distance.
        """
        try:
            return self._pq.lowestPriority()
        except KeyError:
            return None

    def pop(self) -> tuple[Any, Any]:
        """
        Pop the lowest distance cluster pair.

        @raise KeyError: If the distance priority queue is empty.
        @return: A 2-C{tuple} of C{int} cluster numbers.
        """
        return self._pq.pop()

    def __contains__(self, pair: tuple[Any, Any]) -> bool:
        """
        Test if a pair has a computed distance (useful for testing).

        @param pair: A 2-tuple of objects.
        @return: A C{bool} indicating membership.
        """
        a, b = pair
        try:
            self._distances[a][b]
        except KeyError:
            return False
        else:
            return True

    def remove(self, a: Any) -> None:
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
                    # popped as part of the lowest scoring pair.
                    errorCount += 1
                    if errorCount > 1:
                        raise
                del self._distances[b][a]

        del self._distances[a]
