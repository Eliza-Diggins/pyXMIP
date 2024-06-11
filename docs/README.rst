.. image:: docs/source/images/logos_icons/logo_main.png
    :width: 400

|precom| |linting| |docs| |isort Status| |black| |astropy| |pydantic| |sklearn|

pyXMIP (the **p**\ ython **X**\ -**M**\ atching and **I**\ dentification **P**\ ackage) is an astronomical software package for cross referencing source catalogs
against known databases. The package provides a variety of statistical tools for quantifying the confidence of a given match and automates
queries against many of the most common astronomical databases.


Features
========

- Cross match catalogs of sources from survey missions against a wide array of known source databases.
- Use statistical methodologies to model the distribution of sources and produce match probabilities.

Documentation
=============

|docs|

The documentation for ``pyXMIP`` is hosted `here <https://eliza-diggins.github.io/pyXMIP>`_. It includes various example
notebooks, a complete API documentation and further information. It also includes guides on custom code implementations to
meet your science goals.

Installation
============

``pyXMIP`` may be installed from PyPI (stable versions) as follows:

.. code-block:: shell

    >>> pip install pyxmip

To install ``pyXMIP`` from source, you need only clone the git repository

.. code-block:: shell

    >>> git clone https://www.github.com/eliza-diggins/pyXMIP

and then install using

.. code-block:: shell

    >>> cd pyXMIP
    >>> pip install .

You can also achieve this in one line using

.. code-block:: shell

    >>> pip install git+https://www.github.com/eliza-diggins/pyXMIP.git



.. |docs| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg
   :target: https://eliza-diggins.github.io/pyXMIP
.. |precom| image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit
   :target: https://github.com/pre-commit/pre-commit
.. |linting| image:: https://img.shields.io/badge/linting-Flake8-brightgreen.svg?style=flat
.. |Github Page| image:: https://github.com/eliza-diggins/pyXMIP/actions/workflows/build_docs.yml/badge.svg
.. |isort Status| image:: https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336
    :target: https://pycqa.github.io/isort/
    :alt: isort Status
.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
.. |astropy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: https://www.astropy.org
.. |pydantic| image:: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/pydantic/pydantic/main/docs/badge/v2.json
    :target: https://docs.pydantic.dev/latest/
.. |sklearn| image:: http://img.shields.io/badge/powered%20by-sklearn-cyan.svg?style=flat
    :target: https://scikit-learn.org/stable/index.html
