from math import log
from typing import Iterable

from sklearn.utils._array_api import (
    _max_precision_float_dtype,
    get_namespace_and_device,
)

LOG2 = log(2.0)


def entropy(labels):
    """Calculate the entropy for a labeling.

    Parameters
    ----------
    labels : array-like of shape (n_samples,), dtype=int
        The labels.

    Returns
    -------
    entropy : float
       The entropy for a labeling.

    Notes
    -----
    The logarithm used is the natural logarithm (base-e).

    This code is a copy of the 'entropy' function from
    sklearn/metrics/cluster/_supervised.py
    in versions prior to 1.8. It was renamed to _entropy in 1.8 and deprecated
    because it was undocumented. But I (Terry) was using it. So I took a copy on
    2026-03-31 to use it here.
    """
    xp, is_array_api_compliant, device_ = get_namespace_and_device(labels)
    labels_len = labels.shape[0] if is_array_api_compliant else len(labels)

    if labels_len == 0:
        return 1.0

    pi = xp.astype(xp.unique_counts(labels)[1], _max_precision_float_dtype(xp, device_))

    # single cluster => zero entropy
    if pi.size == 1:
        return 0.0

    pi_sum = xp.sum(pi)
    # log(a / b) should be calculated as log(a) - log(b) for
    # possible loss of precision
    # Always convert the result as a Python scalar (on CPU) instead of a device
    # specific scalar array.
    return float(-xp.sum((pi / pi_sum) * (xp.log(pi) - log(pi_sum))))


def entropy2(labels: Iterable[str]) -> float:
    return entropy(labels) / LOG2


MAX_ENTROPY = entropy2(["a", "b", "c", "d"])
