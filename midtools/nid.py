from sklearn.metrics import mutual_info_score
from sklearn.metrics.cluster import entropy

from typing import Iterable

# from numpy import seterr
# seterr(divide='raise')


def normalized_information_distance(c1: Iterable[str], c2: Iterable[str]) -> float:
    """
    Calculate Normalized Information Distance

    Taken from Vinh, Epps, and Bailey, (2010). Information Theoretic Measures
    for Clusterings Comparison: Variants, Properties, Normalization and
    Correction for Chance, JMLR
    <http://jmlr.csail.mit.edu/papers/volume11/vinh10a/vinh10a.pdf>
    """
    denom = max(entropy(c1), entropy(c2))

    # The clusterings are identical (so return 1.0) if both have zero entropy.
    return (1.0 - (mutual_info_score(c1, c2) / denom)) if denom else 1.0
