#!/usr/bin/env python3
"""Pure numpy implementation of Transfer Entropy for binned spike trains.

Replaces pyinform which ships x86_64-only binaries (broken on Apple Silicon).
Implements TE(X->Y) = sum p(y_t+1, y_t^k, x_t^k) * log2(p(y_t+1|y_t^k,x_t^k) / p(y_t+1|y_t^k))

Based on Schreiber (2000) "Measuring information transfer", Phys Rev Lett.
"""

import numpy as np


def transfer_entropy(x: np.ndarray, y: np.ndarray, k: int = 1) -> float:
    """Compute TE from X to Y with history length k.

    Args:
        x: 1D binary array (binned spike train), shape (n_bins,)
        y: 1D binary array (binned spike train), shape (n_bins,)
        k: history length (default: 1)

    Returns:
        Transfer entropy value in bits
    """
    if len(x) != len(y):
        raise ValueError("x and y must have same length")
    if k < 1:
        raise ValueError("k must be >= 1")
    if len(x) < k + 2:
        return 0.0

    n = len(x)
    # State: (x_{t-k+1}, ..., x_t, y_{t-k+1}, ..., y_t)
    # We need y_{t+1} as the target
    # So we create state vectors for t = k to n-2

    # Count occurrences of each (past_state, y_next) combination
    # Past state = tuple of (x_t-k+1...x_t) concatenated with (y_t-k+1...y_t)
    states = {}
    for t in range(k, n - 1):
        past = tuple(x[t - k + 1:t + 1]) + tuple(y[t - k + 1:t + 1])
        y_next = y[t + 1]
        key = (past, y_next)
        states[key] = states.get(key, 0) + 1

    total = sum(states.values())
    if total == 0:
        return 0.0

    # Count marginal occurrences of each past state and y_next
    past_counts = {}
    ynext_counts = {}
    joint_counts = {}
    for (past, y_next), count in states.items():
        past_counts[past] = past_counts.get(past, 0) + count
        ynext_counts[y_next] = ynext_counts.get(y_next, 0) + count
        joint_counts[(past, y_next)] = count

    # TE = sum p(past, y_next) * log2(p(y_next|past) / p(y_next))
    # where p(y_next|past) = p(past, y_next) / p(past)
    # and p(y_next) = marginal

    te = 0.0
    for (past, y_next), count in joint_counts.items():
        p_joint = count / total
        p_past = past_counts[past] / total
        p_ynext = ynext_counts[y_next] / total
        if p_past > 0 and p_ynext > 0 and p_joint > 0:
            te += p_joint * np.log2(p_joint / (p_past * p_ynext) if p_past * p_ynext > 0 else 0)

    # TE can't be negative, clamp to 0
    return max(0.0, te)


def transfer_entropy_matrix(spike_matrix: np.ndarray, k: int = 1) -> np.ndarray:
    """Compute full N×N TE matrix.

    Args:
        spike_matrix: (n_neurons, n_bins) array of binned spike trains
        k: history length

    Returns:
        (n_neurons, n_neurons) TE matrix. Entry (i,j) = TE(i -> j)
    """
    n = spike_matrix.shape[0]
    te = np.zeros((n, n))
    for i in range(n):
        xi = spike_matrix[i]
        for j in range(n):
            if i == j:
                continue
            te[i, j] = transfer_entropy(xi, spike_matrix[j], k=k)
    return te
