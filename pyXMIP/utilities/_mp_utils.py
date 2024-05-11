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


def map_to_threads(mappable, *args, threading_kw=None):
    """
    Perform a certain process ``map`` under certain threading behaviors.
    Parameters
    ----------
    map
    args

    Returns
    -------

    """
    from concurrent.futures import ThreadPoolExecutor

    # ---------------------------------------- #
    # Setup the threading environment          #
    # ---------------------------------------- #
    # Create the kwargs if they don't exist.
    if threading_kw is None:
        threading_kw = {}

    _max_workers = threading_kw.pop("max_workers", 1)
    if _max_workers == 1:
        _threading = False
    else:
        _threading = True

    # ---------------------------------------- #
    # Run the threads                          #
    # ---------------------------------------- #
    if _threading:
        with ThreadPoolExecutor(
            max_workers=_max_workers,
            thread_name_prefix=threading_kw.pop("thread_name_prefix", ""),
        ) as executor:
            results = executor.map(mappable, *args)
    else:
        results = map(mappable, *args)

    return results
