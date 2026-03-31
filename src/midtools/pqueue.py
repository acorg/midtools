import itertools
from typing import Any
from collections import UserList
from collections.abc import Hashable

from heapq import heappush, heappop

_REMOVED = object()


class PriorityQueue(UserList):
    """
    Maintain a priority queue.

    See https://docs.python.org/3.6/library/heapq.html
    """

    def __init__(self) -> None:
        super().__init__(self)
        self._entries: dict[Hashable, list[int | Any]] = {}
        self._counter = itertools.count()
        self._validCount = 0

    def __contains__(self, task: Hashable) -> bool:
        "Is a task in the queue?"
        return task in self._entries

    def __len__(self) -> int:
        "How many valid items are in the queue?"
        return self._validCount

    def add(self, task: Hashable, priority: int = 0):
        "Add a new task or update the priority of an existing task"
        if task in self._entries:
            self.remove(task)
        count = next(self._counter)
        entry = [priority, count, task]
        self._entries[task] = entry
        heappush(self.data, entry)
        self._validCount += 1

    def remove(self, task: Hashable):
        "Mark an existing task as removed. Raise KeyError if not found."
        entry = self._entries.pop(task)
        entry[-1] = _REMOVED
        self._validCount -= 1

    def pop(self, index=0) -> Any:
        "Remove and return the lowest priority task. Raise KeyError if empty."
        while self.data:
            priority, count, task = heappop(self.data)
            if task is not _REMOVED:
                del self._entries[task]
                self._validCount -= 1
                return task
        raise KeyError("pop called on an empty priority queue")

    def lowestPriority(self) -> int:
        "Return the lowest priority item. Raise KeyError if the queue is empty."
        while self.data:
            priority, _, task = self.data[0]
            if task is _REMOVED:
                heappop(self.data)
            else:
                return priority
        raise KeyError("lowestPriority called on an empty priority queue")
