.. _configuration:
========================
Configuring ``pyXMIP``
========================

For most use cases, ``pyXMIP`` doesn't need configuration; however, the code is equipped with a relatively sophisticated set
of configurable options should need / desire drive the user to want something to change. You can access the configuration settings
either through the python module (code-side) or through the terminal commands (terminal-side). In all of the following examples, we
describe both approaches.

Accessing the Configuration
---------------------------

The configuration file for ``pyXMIP`` is located in the ``/bin`` directory of your installation. If you are unsure of your installation
directory, you can find it using the following:

.. tab-set::

    .. tab-item:: Python

        To access the configuration path from ``python``, simply import it from the utilities module:

        .. code-block:: python

            >>> from pyXMIP.utilities.core import config_directory
            >>> print(config_directory)
                PosixPath('.../pyXMIP/bin/config.yaml')

        In your case, you'll see the installation path of your package where the ``...`` are.

    .. tab-item:: CLI

        The configuration path can be obtained using the following CLI command.

        .. code-block:: shell

            USER:~$ pyxmip config path
            .../pyXMIP/bin/config.yaml

The actual configuration appears as an object in ``pyXMIP.utilities.core``. From this object, you can make alterations, print
values, etc.

.. tab-set::

    .. tab-item:: Python

        To access the configuration object, simply import :py:attr:`utilities.core.pxconfig`

        .. code-block:: python

            >>> from pyXMIP.utilities.core import pxconfig

        The actual configuration dictionary is located in ``pxconfig.config``:

        .. code-block:: python

            >>> from rich.pretty import pprint # Just for nicer printing
            >>> pprint(pxconfig.config)
            {
            │   'system': {'preferences': {'disable_progress_bars': False}},
            │   'logging': {
            │   │   'mainlog': {
            │   │   │   'format': '%(name)-3s : [%(levelname)-9s] %(asctime)s %(message)s',
            │   │   │   'level': 'INFO',
            │   │   │   'stream': 'stderr'
            │   │   },
            │   │   'devlog': {
            │   │   │   'enabled': False,
            │   │   │   'format': '%(name)-3s : [%(levelname)-9s] %(asctime)s %(message)s',
            │   │   │   'level': 'DEBUG',
            │   │   │   'stream': 'stderr'
            │   │   }
            │   },
            │   'plotting': {
            │   │   'defaults': {
            │   │   │   'text.usetex': True,
            │   │   │   'xtick.major.size': 8,
            │   │   │   'ytick.major.size': 8,
            │   │   │   'xtick.minor.size': 5,
            │   │   │   'ytick.minor.size': 5,
            │   │   │   'xtick.top': True,
            │   │   │   'xtick.bottom': True,
            │   │   │   'xtick.labeltop': False,
            │   │   │   'xtick.labelbottom': True,
            │   │   │   'xtick.minor.visible': True,
            │   │   │   'ytick.left': True,
            │   │   │   'ytick.right': True,
            │   │   │   'ytick.labelleft': True,
            │   │   │   'ytick.labelright': False,
            │   │   │   'ytick.minor.visible': True,
            │   │   │   'xtick.direction': 'in',
            │   │   │   'ytick.direction': 'in',
            │   │   │   'font.style': 'normal',
            │   │   │   'font.size': 12,
            │   │   │   'image.cmap': 'inferno'
            │   │   }
            │   }
            }

        Configuration settings can be set in the current context simply by making alterations to the dictionary.

    .. tab-item:: CLI

        From the CLI, you can view the current configuration; however, (obviously) you cannot alter the current local context values.
        In the next section, we'll discuss how to make permanent changes to the configuration file.

        To view the configuration or some part of it, you simple need to use the following commands

        .. code-block:: shell

            USER:~$ pyxmip config view
            {
            │   'system': {'preferences': {'disable_progress_bars': False}},
            │   'logging': {
            │   │   'mainlog': {
            │   │   │   'format': '%(name)-3s : [%(levelname)-9s] %(asctime)s %(message)s',
            │   │   │   'level': 'INFO',
            │   │   │   'stream': 'stderr'
            │   │   },
            │   │   'devlog': {
            │   │   │   'enabled': False,
            │   │   │   'format': '%(name)-3s : [%(levelname)-9s] %(asctime)s %(message)s',
            │   │   │   'level': 'DEBUG',
            │   │   │   'stream': 'stderr'
            │   │   }
            │   },
            │   'plotting': {
            │   │   'defaults': {
            │   │   │   'text.usetex': True,
            │   │   │   'xtick.major.size': 8,
            │   │   │   'ytick.major.size': 8,
            │   │   │   'xtick.minor.size': 5,
            │   │   │   'ytick.minor.size': 5,
            │   │   │   'xtick.top': True,
            │   │   │   'xtick.bottom': True,
            │   │   │   'xtick.labeltop': False,
            │   │   │   'xtick.labelbottom': True,
            │   │   │   'xtick.minor.visible': True,
            │   │   │   'ytick.left': True,
            │   │   │   'ytick.right': True,
            │   │   │   'ytick.labelleft': True,
            │   │   │   'ytick.labelright': False,
            │   │   │   'ytick.minor.visible': True,
            │   │   │   'xtick.direction': 'in',
            │   │   │   'ytick.direction': 'in',
            │   │   │   'font.style': 'normal',
            │   │   │   'font.size': 12,
            │   │   │   'image.cmap': 'inferno'
            │   │   }
            │   }
            }

            USER:~$ pyxmip config view system.preferences
            {'disable_progress_bars': False}


Altering The Configuration
--------------------------

When using the ``pxconfig`` variable, changes made to the ``.config`` attribute will not be saved at the end of the runtime context. If you
want to change something permanently, you will need to use a command, alter the file directly, or use the :py:meth:`utilities.Core.YAMLConfiguration.set_on_disk` method.

.. tab-set::

    .. tab-item:: Python

        To make a change on disk using ``python``, simply execute the ``set_on_disk`` function. For example, here
        we change the value of the ``disable_progress_bars`` setting to ``True``.

        .. code-block:: python

            >>> from pyXMIP.utilities.core import pxconfig
            >>> pxconfig.set_on_disk(pxconfig.path,'system.preferences.disable_progress_bars',True)


    .. tab-item:: CLI

        To make the same change using the CLI, simply use the ``config set`` command:

        .. code-block:: shell

            USER:~$ pyxmip config set system.preferences.disable_progress_bars false
