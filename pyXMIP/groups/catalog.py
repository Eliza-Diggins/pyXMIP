"""
Module containing the catalog classes for `pyXs`.
"""
import pathlib as pt

import numpy as np
from astropy import units as units
from astropy.coordinates import SkyCoord
from astropy.table import Table

from pyXMIP.schema import SkyCollectionSchema
from pyXMIP.utilities.core import _enforce_style, mainlog


class SourceCatalog:
    """
    Class representation of the core catalog in pyXs. Catalogs contain observational data that can then be cross-matched
    to databases and other catalogs.
    """

    _default_file_format = "fits"

    def __init__(self, table, image=None, header=None, schema=None):
        """
        Load the :py:class:`catalogs.catalog.SourceCatalog` object.

        Parameters
        ----------
        table: Astropy Table
            The table object containing the sources.
        image: array-like, optional
            The density image if it is provided.
        header: dict, optional
            The object header, if it is available.
        schema: :py:class:`catalogs.catalog.CatalogSchema`, optional
            The schema associated with this source catalog.

        """
        # -- Setting the base parameters -- #
        self.table = table
        self.header = header
        self._schema = schema
        self._detection_density_map = image

        # -- Schema management -- #
        if isinstance(self._schema, dict):
            self._schema = SkyCollectionSchema(self._schema)
        elif isinstance(self._schema, str):
            self._schema = SkyCollectionSchema.from_file(self._schema)
        elif isinstance(self._schema, SkyCollectionSchema):
            pass
        elif self._schema is None:
            self._schema = SkyCollectionSchema.construct(self.table)
        else:
            raise TypeError(f"The schema {self._schema} is invalid.")

        assert self._schema.isvalid()

        # -- additional attribute management -- #
        if not self.header:
            self.header = {}

    def __len__(self):
        return len(self.table)

    def __getitem__(self, item):
        return self.table[self.names == item]

    def __str__(self):
        if self.name is not None:
            return f"<SourceCatalog {self.name}>"
        else:
            return type(self).__name__

    def __repr__(self):
        return self.__str__()

    def __add__(self, other):
        assert isinstance(
            other, SourceCatalog
        ), f"Only SourceCatalog types can be added to SourceCatalog, not {type(other)}."
        from astropy.table import vstack

        return SourceCatalog(
            vstack([self.table, other.table]), schema=self.schema, header=self.header
        )

    def __radd__(self, other):
        return self.__add__(other)

    @property
    def name(self):
        """The name of the catalog object."""
        if "NAME" in self.header:
            return self.header["NAME"]
        else:
            return None

    @name.setter
    def name(self, value):
        self.header["NAME"] = value

    @property
    def schema(self):
        """The schema for this catalog."""
        return self._schema

    @property
    def columns(self):
        """The available columns in the table."""
        return self.table.columns

    def get_coordinates(self):
        """
        Provides the source coordinates from the catalog.

        Returns
        -------
        SkyCoord
            The coordinates on the sky of the objects in the catalog.
        """
        _coord_frame, _coord_cols = self.schema.coordinate_system

        return SkyCoord(
            self.table[_coord_cols[0]], self.table[_coord_cols[1]], _coord_frame
        )

    def names(self):
        """The names of the objects in the catalog."""
        name_schema = self.schema.names
        return self.table[name_schema]

    def redshift(self):
        """The redshifts of the objects in the catalog"""
        _rs_col = self.schema.redshift_column

        if _rs_col is None:
            raise ValueError(f"The schema for {self} indicates no redshift column.")
        else:
            return self.table[self.schema.redshift_column]

    @classmethod
    def from_file(cls, path, format=None):
        """
        Load a catalog from a file.

        .. note::

            Files must be either ``fits`` or ``csv``.

        Parameters
        ----------
        path: str
            The path to the file.
        format: str, optional
            Specified to force a specific formatting choice. Options are ``fits`` and ``csv``. If not provided, the
            format will be deduced from the file extension.

        Returns
        -------
        :py:class:`pyXs.catalogs.catalog.SourceCatalog`

        """
        path = pt.Path(path)
        mainlog.info(f"Loading cluster catalog from {str(path)}...")

        # sanity checks on the path
        assert path.exists(), f"The file {path} could not be found."

        # -- Managing the format type -- #
        if format is None:
            if path.suffix in [".csv"]:
                format = "csv"
            elif path.suffix in [".fits", ".FITS"]:
                format = "fits"
            else:
                raise ValueError(
                    f"The {path} does not have a recognized format and none was specified."
                )

            mainlog.debug(f"Determined {str(path)} was formatted in {format}.")
        else:
            pass

        # -- Loading based on the format -- #
        return getattr(cls, f"_load_from_{format}")(path)

    @classmethod
    def from_url(cls, url):
        raise NotImplementedError()

    @classmethod
    def _load_from_fits(cls, path, extension_index=None):
        from astropy.io import fits

        # -- opening the fits file -- #
        hudl = fits.open(path)

        primary_header = hudl[0].header
        header = {k: primary_header[k] for k in primary_header.keys()}

        if "ISCAT" in primary_header and primary_header["ISCAT"]:
            # This is a standardized catalog file (from pyXs) and we can get extra information out of it.
            mainlog.debug("FITS file recognized as CATALOG.")

            table = hudl["CATALOG"]

            try:
                schema = SkyCollectionSchema.from_file(header["SCMPTH"])
            except Exception:
                mainlog.warning(f"Failed to load schema from {header['SCMPTH']}.")
                schema = None

            if header["HAS_DET"]:
                image = hudl["DETIMG"]
            else:
                image = None

        else:
            # This is a non-standardized format.
            tables = [tbl for tbl in hudl if isinstance(tbl, fits.BinTableHDU)]
            images = [tbl for tbl in hudl if isinstance(tbl, fits.ImageHDU)]

            if extension_index is None:
                table = tables[0]
                if len(tables) > 1:
                    mainlog.warning(
                        f"Defaulting to first available table: {table.name}."
                    )
            else:
                table = tables[extension_index]

            if "DENIMG" in [img.name for img in images]:
                image = hudl["DENIMG"]
            else:
                image = None

            schema = None

        # ---- Loading ----#
        mainlog.info(f"TABLE  ={table.name} from {path}.")
        if schema is not None:
            mainlog.info(f"SCHEMA ={schema.name} from {path}")
        if image is not None:
            mainlog.info(f"DENIMG ={image.name} from {path}")

        ret = cls(
            Table.read(table, format="fits"),
            image=(image.data if image is not None else None),
            header=header,
            schema=schema,
        )

        hudl.close()
        return ret

    def _load_from_csv(self, path):
        raise NotImplementedError()

    def to_file(self, path, overwrite=False):
        """
        Write the catalog object to a fits file.

        Parameters
        ----------
        path: str
            The path to write to.
        overwrite: bool
            If ``True``, the function will be able to overwrite an existing fits file if one is found.

        Returns
        -------
        None

        """
        from astropy.io import fits

        mainlog.info(f"Writing {self} to {path}.")
        # Path management
        path = pt.Path(path)
        if not path.parents[0].exists():
            path.parents[0].mkdir(parents=True, exist_ok=True)

        # -- Write the schema -- #
        self.schema.to_file(str(path) + ".yaml")
        mainlog.info(f"Wrote schema to {path}.yaml.")

        # -- Building the header -- #
        header = self.header
        header["SCMPTH"] = str(path) + ".yaml"
        header["ISCAT"] = True
        header["HAS_DET"] = False

        # -- Collecting fits tables -- #
        table_hduls = []
        _primary_table = fits.table_to_hdu(self.table)
        _primary_table.header.set("name", "CATALOG")

        table_hduls.append(_primary_table)

        # -- Collecting images -- #
        image_hduls = []
        if self._detection_density_map is not None:
            image_hduls.append(
                fits.ImageHDU(self._detection_density_map, name="DETIMG")
            )
            header["HAS_DET"] = True

        # -- Building the header -- #
        _header = fits.Header()

        for k, v in header.items():
            _header[k] = v

        # -- building the primary -- #
        _primary = fits.PrimaryHDU(header=_header)

        hdul = fits.HDUList([_primary] + table_hduls + image_hduls)

        hdul.writeto(path, overwrite=overwrite)
        mainlog.info(f"Wrote {self} to {path}.")

    @_enforce_style
    def plot(
        self,
        color_column=None,
        size_column=None,
        size_range=None,
        norm=None,
        fig=None,
        ax=None,
        projection="aitoff",
        coordinate_system="galactic",
        **kwargs,
    ):
        import matplotlib.pyplot as plt

        # -- Setup the figure / axes -- #
        if fig is None:
            fig = plt.figure(figsize=(8, 8))

        if ax is None:
            ax = fig.add_subplot(111, projection=projection)

        # -- manage the coordinates -- #
        coords = self.get_coordinates()
        if coordinate_system == "galactic":
            phi, theta = (
                coords.galactic.l.wrap_at(180 * units.deg).radian,
                coords.galactic.b.radian,
            )
        elif coordinate_system == "icrs":
            phi, theta = coords.ra.wrap_at(180 * units.deg).radian, coords.dec.radian
        else:
            raise ValueError(
                f"{coordinate_system} is not a known coordinate system for plotting."
            )

        # -- managing the colors -- #
        if color_column is not None:
            if norm is None:
                from matplotlib.colors import Normalize

                norm = Normalize(
                    vmin=np.amin(self.table[color_column]),
                    vmax=np.amax(self.table[color_column]),
                )

        kwargs["c"] = norm(self.table[color_column])
        ax.scatter(phi, theta, **kwargs)

        return fig, ax


if __name__ == "__main__":
    u = SourceCatalog.from_file(
        "/home/ediggins/pyROSITA_test/erass1_hard_eromapper_v1.fits"
    )
    print(u.schema)
    u.schema.redshift_column = "BEST_Z"

    print(u.redshift())
