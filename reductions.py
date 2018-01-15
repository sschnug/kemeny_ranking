import numpy as np
import scipy.sparse as sp
from utils import combs


def extended_condorcet_simple(rankings):
    # assumes: cands -> 0,N-1
    n = rankings.shape[1]
    cands = np.arange(n)
    pairs = combs(range(n), 2)

    condorcet_rows, condorcet_cols = [], []

    for cand, other_cand in pairs:
        cand_pos = np.where(rankings == cand)[1]
        other_pos = np.where(rankings == other_cand)[1]

        if np.all(cand_pos < other_pos):
            condorcet_rows.append(cand)
            condorcet_cols.append(other_cand)
        elif np.all(other_pos < cand_pos):
            condorcet_rows.append(other_cand)
            condorcet_cols.append(cand)

    mat = sp.coo_matrix((np.ones(len(condorcet_rows)), (condorcet_rows, condorcet_cols)))
    return mat
