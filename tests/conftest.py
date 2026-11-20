import pytest


@pytest.fixture(autouse=True)
def _isolated_interface_memory(tmp_path_factory, monkeypatch):
    """The interface remembers the last mode and folders in the user's home;
    tests use a private location instead."""
    monkeypatch.setenv("MSPARSER_HOME", str(tmp_path_factory.mktemp("msparser_home")))
