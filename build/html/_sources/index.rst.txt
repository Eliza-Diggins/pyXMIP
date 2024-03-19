.. image:: images/logos_icons/logo_main.png
    :scale: 50%

|precom| |docs| |isort Status| |black| |astropy| |astroQuery|

pyXMIP (the python X-Matching and Identification Package) is an astronomical software package for cross referencing source catalogs
against known databases. The package provides a variety of statistical tools for quantifying the confidence of a given match and automates
queries against many of the most common astronomical databases.



.. raw:: html

   <hr style="color:black">

Features
========

- Cross match catalogs of sources from survey missions against a wide array of known source databases.
- Use statistical methodologies to model the distribution of sources and produce match probabilities.
- Create mock all sky surveys.

.. grid:: 3

    .. grid-item::

        .. dropdown:: Astronomical Databases

            - |NED| IPAC/NED
            - |SIMBAD|

    .. grid-item::

        .. dropdown:: Tools

            Coming Soon

    .. grid-item::

        .. dropdown:: Statistic Methods

            Coming Soon

Resources
=========

.. grid:: 2
    :padding: 3
    :gutter: 5

    .. grid-item-card::
        :img-top: images/logos_icons/stopwatch_icon.png

        Quickstart Guide
        ^^^^^^^^^^^^^^^^
        New to pyXMIP? The quickstart guide page includes a variety of examples for utilizing our software for easy use cases.

        +++

        .. button-ref:: getting_started
            :expand:
            :color: secondary
            :click-parent:

            To The Quickstart Page

    .. grid-item-card::
        :img-top: images/logos_icons/lightbulb.png

        Examples
        ^^^^^^^^
        Interested in seeing a variety of use cases for the software? This example catalog is worth checking out!

        +++

        .. button-ref:: examples
            :expand:
            :color: secondary
            :click-parent:

            To the Examples Page

    .. grid-item-card::
        :img-top: images/logos_icons/book.svg

        User References
        ^^^^^^^^^^^^^^^^

        The user reference page contains a brief overview of the critical topics relating to the development and use of this
        software.
        +++

        .. button-ref:: user_guide
            :expand:
            :color: secondary
            :click-parent:

            To the User Guide

    .. grid-item-card::
        :img-top: images/logos_icons/api_icon.png

        API Reference
        ^^^^^^^^^^^^^

        Doing a deep dive into our code? Looking to contribute to development? The API reference is a comprehensive resource
        complete with source code and type hinting so that you can find every detail you might need.

        +++

        .. button-ref:: api
            :expand:
            :color: secondary
            :click-parent:

            API Reference



.. raw:: html

   <hr style="color:black">

Pages
-----

.. toctree::
   :maxdepth: 1

   api
   user_guide
   getting_started
   examples


.. |docs| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg
   :target: https://jzuhone.github.io/cluster_generator/build/html/index.html
.. |precom| image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit
   :target: https://github.com/pre-commit/pre-commit
.. |Github Page| image:: https://github.com/eliza-diggins/pyXs/actions/workflows/build_docs.yml/badge.svg
.. |astropy| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg
.. |isort Status| image:: https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336
    :target: https://pycqa.github.io/isort/
    :alt: isort Status
.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
.. |astropy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/
    :target: https://www.astropy.org
.. |astroquery| image:: http://img.shields.io/badge/powered%20by-AstroQuery-orange.svg?style=flat)](http://www.astropy.org/
    :target: https://www.astropy.org
.. |NED| image:: images/logos_icons/NED.png
    :scale: 20%
    :target: https://ned.ipac.caltech.edu
.. |SIMBAD| image:: images/logos_icons/SIMBAD.jpg
    :scale: 20%
    :target: https://simbad.cds.unistra.fr/simbad/
