"""
Utility module for optimized operations on datasets or via parallelism
"""
from multiprocessing.shared_memory import SharedMemory
from typing import Callable, Collection, Iterator, TypeVar

import numpy as np

_T = TypeVar("_T")


def split(a: Collection[_T], n: int) -> Collection[Collection[_T]]:
    """
    Split the collection ``a`` into ``n`` chunks of near-equally sized data.

    Parameters
    ----------
    a: list
        The list to split into chunks.
    n: int
        The number of chunks.

    Returns
    -------
    list of list
        A nested list of each chunk.

    """
    k, m = divmod(len(a), n)
    return [a[i * k + min(i, m) : (i + 1) * k + min(i + 1, m)] for i in range(n)]


def created_shared_memory_equivalent(array: np.ndarray) -> SharedMemory:
    """
    Create a shared memory equivalent of an array. Effectively produces a pointer-like variable directed
    to the position in memory where this array is stored.

    Parameters
    ----------
    array: array
        The array to generate a shared memory equivalent of.

    Returns
    -------
    SharedMemory
        The buffer specification for the data in memory.

    """
    header = SharedMemory(create=True, size=array.nbytes)
    copy = np.ndarray(array.shape, dtype=array.dtype, buffer=header.buf)
    copy[:] = array[:]
    return header


def map_to_threads(
    mappable: Callable[[...], ...], *args, threading_kw: dict = None
) -> Iterator | map:
    """
    Map a function (``mappable``) with arguments ``*args`` to threads or run without threads if threading is
    not enabled.

    Parameters
    ----------
    mappable: Callable
        The mapping function which acts on ``args`` and returns a result.
    *args:
        Arguments to pass through ``mappable``. These would be any of the standard arguments.
    threading_kw: dict
        Dictionary containing the threading parameters. TODO: threading kwargs.

    Returns
    -------
    Iterator
        The output of the mapping operation. Equivalent to ``map(mappable,*args)``.

    Notes
    -----

    This function utilizes the ``concurrent.futures`` ``ThreadPoolExecutor``.
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
