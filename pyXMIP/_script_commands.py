"""
Command's specific to the CLI interface.
"""
from argparse import ArgumentParser
from typing import Mapping

from pyXMIP.utilities.text import get_package_version, print_version

CLI_FUNCTIONS = {f.__name__: f for f in [get_package_version, print_version]}


class CLIParser(ArgumentParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def build_recursive(cls, parser_config: Mapping, base_parser=None):
        if base_parser is None:
            base_parser = cls(prog="pyXMIP")

            _w = base_parser.add_subparsers(title="commands", help="Available commands")

            cls.build_recursive(parser_config, base_parser=_w)

            return base_parser

        # -- Iterate through the layers -- #
        for k, v in parser_config.items():
            _is_command = (
                "subparsers" not in v.keys()
            )  # determines if its a command or not.

            if _is_command:
                # This entry is a command.
                _s = base_parser.add_parser(name=k, help=v.get("help", None))

                for arg_name, arg_dict in v.get("args", {}).items():
                    _s.add_argument(arg_name, **arg_dict)
                for kwarg_name, kwarg_dict in v.get("kwargs", {}).items():
                    _s.add_argument(
                        f"--{kwarg_name}",
                        kwarg_dict.pop("shortcut", f"-{kwarg_name}"),
                        **kwarg_dict,
                    )

                _s.set_defaults(operation=v.get("function", "none"))

            else:
                _s = base_parser.add_parser(name=k, help=v.get("help", None))
                _q = _s.add_subparsers(title=v.get("title"), help=v.get("help", None))

                cls.build_recursive(v["subparsers"], base_parser=_q)

        return base_parser


def print_config_path():
    from pyXMIP.utilities.core import config_directory

    print(str(config_directory))


def view_config(dictionary_position="all"):
    from pyXMIP.utilities.core import pxconfig, rgetattr

    if dictionary_position == "all":
        print(pxconfig.config)
    else:
        print(rgetattr(pxconfig.config, dictionary_position))


def set_config(dictionary_position=None, value=None):
    from yaml import safe_load

    from pyXMIP.utilities.core import pxconfig, rgetattr
    from pyXMIP.utilities.logging import mainlog

    if isinstance(rgetattr(pxconfig.config, dictionary_position), dict):
        mainlog.error(
            "Cannot alter option dictionary from CLI, please provide a valid option path."
        )
    else:
        pxconfig.set_on_disk(pxconfig.path, dictionary_position, safe_load(value))


CLI_FUNCTIONS["view_config"] = view_config
CLI_FUNCTIONS["set_config"] = set_config
CLI_FUNCTIONS["print_config_path"] = print_config_path
