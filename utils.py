import itertools
import numpy as np


def combs(a, r):
    """
    Return successive r-length combinations of elements in the array a.
    Should produce the same output as array(list(combinations(a, r))), but
    faster.
    """
    a = np.asarray(a)
    dt = np.dtype([('', a.dtype)]*r)
    b = np.fromiter(itertools.combinations(a, r), dt)
    b_ = b.view(a.dtype).reshape(-1, r)
    return b_

def perms(a, r):
    """
    Same as above with permutations
    """
    a = np.asarray(a)
    dt = np.dtype([('', a.dtype)]*r)
    b = np.fromiter(itertools.permutations(a, r), dt)
    b_ = b.view(a.dtype).reshape(-1, r)
    return b_
