.. databases_
==========
Databases
==========

pyXMIP's core functionality is the cross-mapping / cross-referencing of source data against known databases of astronomical objects.
The :py:class:`pyXMIP.structures.databases.SourceDatabase` class is the corner-stone of this functionality. On this page, we'll describe
how these objects work and how you can utilize them to meet your scientific needs.

Basic Information
-----------------

Every database object (:py:class:`pyXMIP.structures.databases.RemoteDatabase` or :py:class:`pyXMIP.structures.databases.LocalDatabase`) is effectively
a simplified structure for interacting with :py:mod:`astroquery`, and from there with the query URLs for each database. Databases have 1 key method; the
:py:meth:`pyXMIP.structures.databases.SourceDatabase.query_radius`, which allows the user to search for any matches in that database around the specified
position and within a given radius. Running a query returns an astropy :py:class:`astropy.table.Table` instance.

.. code-block:: python

    >>> from pyXMIP.structures.databases import NED
    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units
    >>> # Create the necessary positions
    >>> position = SkyCoord(ra=83.63,dec=22.0144,unit='deg')
    >>> search_radius = 1*units.arcmin
    >>> # Conduct the query.
    >>> output_table = NED.query_radius(position,search_radius)
