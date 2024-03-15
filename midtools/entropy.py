from math import log
from sklearn.metrics.cluster import entropy

LOG2 = log(2.0)


def entropy2(labels):
    return entropy(labels) / LOG2


MAX_ENTROPY = entropy2(["a", "b", "c", "d"])
