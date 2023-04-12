"""
Tests for command line interface (CLI)
"""
from importlib import import_module
from importlib_metadata import version
from os import linesep
from unittest.mock import patch

import entirety
import pytest

from cli_test_helpers import ArgvContext, shell



def test_runas_module():
    """
    Can this package be run as a Python module?
    """
    result = shell('entirety --help')
    assert result.exit_code == 0



def test_main_module():
    """
    Exercise (most of) the code in the ``__main__`` module.
    """
    import_module('entirety')





def test_version():
    """
    Does --version display information as expected?
    """
    expected_version = version('entirety')
    result = shell('entirety --version')

    assert result.stdout == 'entirety 'f"{expected_version}{linesep}"
    assert result.exit_code == 0

'''
@patch('entirety')
def test_usage(mock_dispatch):
    """
    Does CLI abort w/o arguments, displaying usage instructions?
    """
    with ArgvContext('entirety'), pytest.raises(SystemExit):
        entirety()

    assert not mock_dispatch.called, 'CLI should stop execution'

    result = shell('entirety')

    assert 'usage:' in result.stderr
'''
