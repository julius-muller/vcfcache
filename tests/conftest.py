# tests/conftest.py
import os
import pytest
from pathlib import Path


@pytest.fixture(scope="session", autouse=True)
def set_vepstash_root():
    """Set the VEPSTASH_ROOT environment variable for all tests."""
    if 'VEPSTASH_ROOT' not in os.environ:
        package_root = Path(__file__).parent.parent.absolute()
        os.environ['VEPSTASH_ROOT'] = str(package_root)
        print(f"Set VEPSTASH_ROOT to {package_root}")
    return os.environ['VEPSTASH_ROOT']