.. _object_type_reduction:
===========================================
Object Type Reduction
===========================================

One of the most obvious reduction processes that might come to mind is to eliminate potential matches which, given the parameters
of the instruments involved, are unlikely to actually be detected in the band-pass of the observing telescope. In other words, eliminating
unlikely candidates for a given pattern of emission.

This is the intention of the :py:class:`reduction.ObjectTypeReduction` process, which takes into account characteristic behaviors of
various classes of astronomical objects along with the parameters of the observing mission to determine the likelihood of a given match. In this
reference, we will discuss in-depth the core concepts behind this procedure.

Modeling Astronomical Objects
-----------------------------

A seemingly obvious first question: how does one go about modeling the emission behavior of a given class of astronomical objects?

To do this in ``pyXMIP``, we provide so-called continuum models of each of our recognized object types. The continuum model is a probability density function
:math:`f(L,\lambda;\Theta)`, where :math:`\lambda` is the emission wavelength, :math:`L` is the object luminosity in that band, and :math:`f(L,\lambda;\Theta)` is the
probability of such a source existing subject to our model parameters :math:`\Theta`.



Modeling Instruments
--------------------
