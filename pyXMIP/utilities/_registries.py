"""
Registry classes for storage and tracking of class instances.
"""


class _Registry(dict):
    """
    Bare-bones registry class.
    """

    def __init__(self, mapping):
        """
        Initialize the generalized registry class.

        Parameters
        ----------
        mapping: dict
            The dictionary underlying this registry.
        """
        super().__init__(mapping)

    def as_list(self):
        return list(self.values())

    def add(self, key, value):
        assert key not in self, f"Key {key} already in registry."
        self[key] = value
