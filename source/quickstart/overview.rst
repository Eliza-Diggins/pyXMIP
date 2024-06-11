.. image:: ../images/logos_icons/logo_main.png
    :scale: 50%

.. _overview:
=====================
Overview
=====================


``pyXMIP`` is an all-purpose astronomical software package for cross-referencing, source identification, and catalog production. You
can use ``pyXMIP`` to catalog novel telescoping sources, check if a source is an unknown object, and even identify different types of
astronomical objects. In addition, ``pyXMIP`` uses a comprehensive statistical framework to quantify the confidence in identifications and
catalog generation.

What Can I Do With ``pyXMIP``?
-------------------------------

- Rapidly and uniformly perform position-cone queries to a variety of common astronomical databases.
- Cross match tables of sources against any number of astronomical databases.
- Perform match reduction processes to identify unique (best-match) candidates for a given source object.
- Account for instrumentation parameters, Poissonian source rates, dust mapping, and object type constraints at all stages of catalog generation.

How Does it Work?
-----------------

``pyXMIP`` is an all-purpose package and is thus composed of a variety of sub-modules for accomplishing different tasks. Here, we'll briefly summarize
some of the core functionality and provide links for further reading.

In a general sense, ``pyXMIP`` provides a toolkit for working with unidentified astronomical sources. We provide tools for identifying sources from
all types of astronomical study. This is predominantly done by referencing known databases against the provided sources and using a mix of statistics, physics, and
instrumental knowledge to identify likely match candidates. ``pyXMIP`` will continue to grow to include other tasks like characterizing spectra of objects,
performing cluster identification, and more.

.. image:: ../images/diagrams/purpose_flow.png
