"""Standard utilities for multiprocessing and multithreading.
"""
from multiprocessing.shared_memory import SharedMemory

import numpy as np


def split(a, n):
    k, m = divmod(len(a), n)
    return [a[i * k + min(i, m) : (i + 1) * k + min(i + 1, m)] for i in range(n)]


def created_shared_memory_equivalent(array):
    header = SharedMemory(create=True, size=array.nbytes)
    copy = np.ndarray(array.shape, dtype=array.dtype, buffer=header.buf)
    copy[:] = array[:]
    return header
