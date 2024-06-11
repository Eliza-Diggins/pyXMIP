"""
Pytest configuration for the ``pyXMIP`` testing suite.

Notes
-----

There are a few command line arguments that can be provided when running the testing suite:

- ``answer_dir`` (Required) The directory in which the answers are stored.
- ``answer_store`` If this flag is active, we create new answers in ``answer_dir`` and don't check answers.
- ``runslow`` If active, slow tests will also be run.
"""
import os

import pytest


def pytest_addoption(parser):
    parser.addoption("--answer_dir", help="Directory where answers are stored.")
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )
    parser.addoption(
        "--answer_store",
        action="store_true",
        help="Generate new answers, but don't test.",
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    # ------------------------------------------ #
    # Managing slow skipping                     #
    # ------------------------------------------ #
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


@pytest.fixture()
def answer_store(request):
    """fetches the ``--answer_store`` option."""
    return request.config.getoption("--answer_store")


@pytest.fixture()
def answer_dir(request):
    """fetches the ``--answer_dir`` option."""
    ad = os.path.abspath(request.config.getoption("--answer_dir"))
    if not os.path.exists(ad):
        os.makedirs(ad)
    return ad
