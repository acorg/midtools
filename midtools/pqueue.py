import itertools

from heapq import heappush, heappop

REMOVED = None


class PriorityQueue(list):
    """
    Maintain a priority queue.

    See https://docs.python.org/3.6/library/heapq.html
    """

    def __init__(self):
        list.__init__(self)
        self._entries = {}
        self._counter = itertools.count()
        self._validCount = 0

    def __contains__(self, task):
        "Is a task in the queue?"
        return task in self._entries

    def __len__(self):
        "How many valid items are in the queue?"
        return self._validCount

    def add(self, task, priority=0):
        "Add a new task or update the priority of an existing task"
        if task in self._entries:
            self.remove(task)
        count = next(self._counter)
        entry = [priority, count, task]
        self._entries[task] = entry
        heappush(self, entry)
        self._validCount += 1

    def remove(self, task):
        "Mark an existing task as removed. Raise KeyError if not found."
        entry = self._entries.pop(task)
        entry[-1] = REMOVED
        self._validCount -= 1

    def pop(self):
        "Remove and return the lowest priority task. Raise KeyError if empty."
        while self:
            priority, count, task = heappop(self)
            if task is not REMOVED:
                del self._entries[task]
                self._validCount -= 1
                return task
        raise KeyError("pop from an empty priority queue")

    def lowestPriority(self):
        "return the lowest priority. Raise KeyError if empty"
        while self:
            priority, _, task = self[0]
            if task is REMOVED:
                heappop(self)
            else:
                return priority
        raise KeyError("peek on an empty priority queue")
