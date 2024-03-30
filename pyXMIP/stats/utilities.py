"""
Statistical Utility functions
"""
import numpy as np


def uniform_sample_spherical(n_points):
    r"""
    Return a uniform sample from a spherical surface.

    Parameters
    ----------
    n_points: int
        The number of samples to draw.

    Returns
    -------
    phi
        The azimuthal positions (in rad)
    theta
        The elevation angles (in rad)

    Notes
    -----

    This function relies on functional transform to produce the random sample. A coordinate segment of a spherical surface :math:`d\Omega` has
    area

    .. math::

        d\Omega = \sin(\theta) d\phi d\theta

    Thus, if the probability of finding :math:`k` points in a region of size :math:`d\Omega`, then

    .. math::

        P(k|d\Omega) \sim U(d\Omega) \sim \sin(\theta) U(d\phi d\theta) \sim \alpha d\psi d\omega.

    To accomplish such a transformation, we require only that

    .. math::

        \psi(\phi,\theta) = \phi,\; \omega(\phi,\theta) = \sin^{-1} \theta

    Thus, for the distribution to be independent of position we simply let

    .. math::

        \cos \phi \sim U(N).
    """
    return 2 * np.pi * np.random.uniform(0, 1, size=n_points), np.arccos(
        2 * np.random.uniform(0, 1, size=n_points) - 1
    )
