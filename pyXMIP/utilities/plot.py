"""
Plotting utilities for use in the pyXMIP backend.
"""
import functools

from matplotlib import pyplot as plt

from pyXMIP.utilities.core import xsparams


# ======================================================================================================================#
# STYLE FUNCTIONS                                                                                                      #
# ======================================================================================================================#
def _enforce_style(func):
    """Enforces the mpl style."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        _rcp_copy = plt.rcParams.copy()

        for _k, _v in xsparams["plotting"]["defaults"].items():
            plt.rcParams[_k] = _v

        out = func(*args, **kwargs)

        plt.rcParams = _rcp_copy
        del _rcp_copy

        return out

    return wrapper


def set_style():
    for _k, _v in xsparams["plotting"]["defaults"].items():
        plt.rcParams[_k] = _v


# ======================================================================================================================#
# Plotting Systems                                                                                                      #
# ======================================================================================================================#
