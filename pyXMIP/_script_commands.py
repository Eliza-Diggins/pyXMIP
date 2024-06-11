"""
CLI specific commands
"""
from pyXMIP.utilities._text import (
    get_citation,
    get_package_version,
    print_cli_header,
    show_config,
)
from pyXMIP.utilities.core import config_directory, getFromDict, setInDict, xsparams

try:
    from rich.pretty import pprint as print
except ImportError:
    pass


def build_cli_parsers_recursive(yaml_dict, parent_parser=None):
    """Build CLI parsers in a recursive manner."""
    _output = {}

    for k, v in yaml_dict.items():
        if not isinstance(v, dict):
            # We want all of the leaves to be dictionaries.
            raise ValueError(
                "Encountered non-dictionary leaf ending on CLI configuration."
            )

        _is_leaf = v.get("is_leaf", False)  # --> tells us if this is a leaf.

        if not _is_leaf:
            # setup the command parser then proceed to next level
            _output[k] = {
                "parser": parent_parser.add_parser(name=k, help=v.get("help", None))
            }
            _output[k]["command_parser"] = _output[k]["parser"].add_subparsers(
                title="commands", help="Available sub-commands."
            )

            # -- now pass to the next level -- #
            _output[k]["subparsers"] = build_cli_parsers_recursive(
                v["leaves"], _output[k]["command_parser"]
            )
        else:
            _output[k] = {
                "parser": parent_parser.add_parser(name=k, help=v.get("help", None))
            }

            for arg in v["args"]:
                args, kwargs = arg
                _output[k]["parser"].add_argument(*args, **kwargs)

            _output[k]["parser"].set_defaults(
                operation_function=v.get("function", "unimplemented")
            )

    return _output


# ----------------------------------------------------------- #
# Basic Functions                                             #
# ----------------------------------------------------------- #
def print_version(**kwargs):
    """print the current pyXMIP version"""
    print(get_package_version())


def cite(**kwargs):
    print(print_cli_header())
    print(get_citation())


# ----------------------------------------------------------- #
# CONFIG Functions                                            #
# ----------------------------------------------------------- #
def print_config(**kwargs):
    show_config()


def set_config(option_name=None, option_value=None):
    from ruamel.yaml import YAML

    yaml = YAML()
    with open(config_directory, "r+") as f:
        config = yaml.load(f)

        try:
            setInDict(config, option_name, option_value)
            print(f"Set config {option_name} to {option_value}.")
        except KeyError:
            print("Failed to change configuration.")

    with open(config_directory, "w") as f:
        yaml.dump(config, f)


def get_config(option_name=None):
    if not isinstance(option_name, list):
        option_name = [option_name]
    try:
        print(getFromDict(xsparams, option_name))
    except KeyError:
        print("No configuration setting item.")


# ----------------------------------------------------------- #
# Database Commands                                           #
# ----------------------------------------------------------- #
def get_available_databases(**kwargs):
    from pyXMIP.structures.databases import DEFAULT_DATABASE_REGISTRY

    print(DEFAULT_DATABASE_REGISTRY)


def get_poisson_path(database=None):
    from pyXMIP.structures.databases import DEFAULT_DATABASE_REGISTRY

    print(DEFAULT_DATABASE_REGISTRY.get(database).poisson_map)


# ----------------------------------------------------------- #
# X-match Commands                                            #
# ----------------------------------------------------------- #
def xmatch_run(input_path=None, output_path=None, databases=None, force=None):
    from pyXMIP.cross_reference import cross_match

    cross_match(input_path, output_path, databases=databases, overwrite=force)
