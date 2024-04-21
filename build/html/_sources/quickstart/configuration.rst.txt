.. _configuration:
========================
Configuring ``pyXMIP``
========================

``pyXMIP`` has a lot of configurable settings that you can interact with. In this mini-guide, we'll show you the basics
of messing around with the settings for your installation.

To start, let's take a look at our configuration settings:

.. code-block:: shell
    
    >>> pyxmip config show

        {
        'system': {
            'logging': {
                'main': {
                    'enabled': True,
                    'format': '%(name)-3s : [%(levelname)-9s] %(asctime)s %(message)s',
                    'level': 'DEBUG',
                    'stream': 'STDERR'
                },
                'developer': {
                    'enabled': False,
                    'output_directory': None
                }
            },
            'add_ons': {
                'rich': False
            }
        },
        'plotting': {
            'defaults': {
                'text.usetex': True,
                'xtick.major.size': 8,
                'ytick.major.size': 8,
                'xtick.minor.size': 5,
                'ytick.minor.size': 5,
                'xtick.top': True,
                'xtick.bottom': True,
                'xtick.labeltop': False,
                'xtick.labelbottom': True,
                'xtick.minor.visible': True,
                'ytick.left': True,
                'ytick.right': True,
                'ytick.labelleft': True,
                'ytick.labelright': False,
                'ytick.minor.visible': True,
                'xtick.direction': 'in',
                'ytick.direction': 'in',
                'font.style': 'normal',
                'font.size': 12,
                'image.cmap': 'inferno'
            }
        }
    }

There's a lot going on there!

.. adonition:: TODO

    TODO: describe in better detail how to interact with the configuration information.
