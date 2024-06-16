"""
Module with core functionality for managing different source types in ``pyXMIP``.
"""
import os
import pathlib as pt
from typing import Mapping, TypeVar

import ruamel.yaml as yaml
from pydantic import BaseModel, model_validator
from tqdm.auto import tqdm

from pyXMIP.utilities.core import bin_directory
from pyXMIP.utilities.types import PydanticRegistry, Registry

# -- Importing SELF -- #
try:
    from typing import Self  # noqa
except ImportError:
    from typing_extensions import Self as Self  # noqa

# ------------------------------------------------------------------------ #
# Core Variables                                                           #
# ------------------------------------------------------------------------ #
yaml = yaml.YAML()
ASS = TypeVar("ASS")
_default_objects_path = pt.Path(os.path.join(bin_directory, "objects.yaml"))


class SourceTypeRegistry(Registry):
    """
    Registry class for collections of astronomical source systems.
    """

    def register(self, source_object: ASS, overwrite: bool = False) -> None:
        if source_object.name in self.keys():
            # The source object is at least in the key set. We still need to check for equality.
            if not source_object == self[source_object.name]:
                assert (
                    overwrite
                ), f"Found a matching registry item {source_object.name}, but overwrite = False and the matching item was not that passed."

                self[source_object.name] = source_object
        else:
            self[source_object.name] = source_object

    @classmethod
    def from_yaml(cls, path: str | pt.Path) -> Self:
        """
        Load the :py:class:`SourceTypeRegistry` from a YAML file.

        Parameters
        ----------
        path: str
            The path to the file.

        Returns
        -------
        :py:class:`SourceTypeRegistry`
            The :py:class:`SourceTypeRegistry` loaded from the YAML file.
        """
        # -------------------------------------- #
        # Validating the file and arg management #
        # -------------------------------------- #
        path = pt.Path(path)  # type enforcement.
        assert path.exists(), f"{path} does not exist"

        with open(path, "r") as stream:
            data = yaml.load(stream)

        # -------------------------------------- #
        # Loading the registry                   #
        # -------------------------------------- #
        return cls.load_from_dict(data)

    @classmethod
    def load_from_dict(cls, mapping: Mapping, instance: Self = None):
        if instance is None:
            instance = cls({})

        for k, v in tqdm(mapping.items()):
            instance[k] = AstronomicalSourceSystem(name=k, **v, registry=instance)

        return instance

    def __contains__(self, item: str | Self):
        if isinstance(item, AstronomicalSourceSystem):
            return item.name in self.keys()
        else:
            return item in self.keys()


class AstronomicalSourceSystem(BaseModel):
    """
    Base class representing a generic astronomical source type.
    """

    name: str
    """ str: The name of this astronomical source type."""
    description: str = None
    """ str: The description of this AstronomicalSourceSystem."""
    aliases: list[str] = None
    """ list of str: The recognized aliases for this astronomical source type."""
    registry: PydanticRegistry | None = None
    """ :py:class:`SourceTypeRegistry`: The registry in which to register this object."""
    parents: list[Self | str] = None
    """ list of :py:class:`AstronomicalSourceSystem`: The parents of this astronomical source system."""
    children: list[Self | str] = None
    """ list of :py:class:`AstronomicalSourceSystem`: The descendants of this astronomical source system."""

    @model_validator(mode="after")
    def validator(self):
        if self.aliases is None:
            self.aliases = []
        if self.parents is None:
            self.parents = []
        if self.children is None:
            self.children = []

        # Registry management
        if self.registry is not None:
            # We are going to check for contains and then add.
            self.registry.register(self, overwrite=True)

            # -- register parents and children -- #
            # This is necessary to enforce self-consistent typing norms.
            # Otherwise the types would potentially have missing instances on either side of the instance.
            if isinstance(self.registry, SourceTypeRegistry):
                self.parents = [
                    self.registry[p] for p in self.parents if p in self.registry
                ]
                self.children = [
                    self.registry[p] for p in self.children if p in self.registry
                ]

                for child in self.children:
                    if self not in child.parents:
                        child.parents.append(self)

                for parent in self.parents:
                    if self not in parent.children:
                        parent.children.append(self)

        return self

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.__str__()


DEFAULT_OBJECT_REGISTRY = SourceTypeRegistry.from_yaml(_default_objects_path)

if __name__ == "__main__":
    print(DEFAULT_OBJECT_REGISTRY["G"].children)
