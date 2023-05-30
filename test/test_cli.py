"""
Tests for command line interface (CLI)
"""
from importlib import import_module
from importlib_metadata import version
from os import linesep
from unittest.mock import patch

from cli_test_helpers import shell



def test_runas_module():
    """
    Can this package be run as a Python module?
    """
    result = shell('chromTools --help')
    assert result.exit_code == 0



def test_main_module():
    """
    Exercise (most of) the code in the ``__main__`` module.
    """
    import_module('chromTools')



def test_version():
    """
    Does --version display information as expected?
    """
    expected_version = version('chromTools')
    result = shell('chromTools --version')

    assert result.stdout == 'chromTools 'f"{expected_version}{linesep}"
    assert result.exit_code == 0


