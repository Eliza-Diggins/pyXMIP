"""
Logging configuration for ``pyXMIP``.

Notes
-----

``pyXMIP`` has the following loggers:

- ``mainlog``: the main logger for the module. This is used for any standard display needs within any
  of the processes.
- ``devlog``: The development logger. Produces logging statements which are more useful for debugging.


"""
import logging
import sys

from pyXMIP.utilities.core import pxconfig

# --------------------------------------------------------- #
# Establishing relevant streams                             #
# --------------------------------------------------------- #
streams = dict(
    mainlog=getattr(sys, pxconfig.config.logging.mainlog.stream),
    devlog=getattr(sys, pxconfig.config.logging.devlog.stream),
)

_loggers = dict(mainlog=logging.Logger("pyXMIP"), devlog=logging.Logger("pyXMIP-DEV"))

_handlers = {}

for k, v in _loggers.items():
    # Construct the formatter string.
    _handlers[k] = logging.StreamHandler(streams[k])
    _handlers[k].setFormatter(
        logging.Formatter(getattr(pxconfig.config.logging, k).format)
    )
    v.addHandler(_handlers[k])
    v.setLevel(getattr(pxconfig.config.logging, k).level)
    v.propagate = False

    if k != "mainlog":
        v.disabled = not getattr(pxconfig.config.logging, k).enabled

mainlog: logging.Logger = _loggers["mainlog"]
# :py:class:`logging.Logger`: The main logger for ``pyXMIP``.
devlog: logging.Logger = _loggers["devlog"]
# :py:class:`logging.Logger`: The development logger for ``pyXMIP``.
