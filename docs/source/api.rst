API
===

API documentation for each of the modules in the ``pyXMIP`` library can be found below. Relevant user guides are linked
in the API documentation for which they are relevant.

Core Modules
------------

Below are the central datatype modules of the ``pyXMIP`` module. These modules include the core types that users should
be familiar with and will interact with most often.

.. autosummary::
    :toctree: _as_gen
    :recursive:
    :template: module.rst
    :nosignatures:

    structures.databases
    structures.table
    structures.map
    structures.reduction
    schema
    cross_reference

Statistics
----------

.. autosummary::
    :toctree: _as_gen
    :recursive:
    :template: module.rst
    :nosignatures:

    stats.utilities
    stats.gaussian_process
    stats.map_regression


Other
-----

These sub-modules are largely collections of convenience functions or utility functions for use elsewhere. Nonetheless,
users interested in the logging system, the config system, or plotting should look at the relevant documentation here.

.. autosummary::
    :toctree: _as_gen
    :recursive:
    :template: module.rst
    :nosignatures:

    utilities.core
    utilities.types
    utilities.logging
    utilities.optimize
    utilities.plot
    utilities.text
    utilities.sql
