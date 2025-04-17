import os
import tempfile
import pytest
from pathlib import Path
from src.utils.validation import compute_md5

def test_compute_md5():
    """Test that compute_md5 correctly calculates MD5 hash of a file."""
    # Create a temporary file with known content
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file.write(b"test content for MD5 calculation")
        temp_file_path = temp_file.name
    
    try:
        # Calculate MD5 hash
        calculated_md5 = compute_md5(Path(temp_file_path))
        
        # Expected MD5 hash for "test content for MD5 calculation"
        # This was pre-calculated using: echo -n "test content for MD5 calculation" | md5sum
        expected_md5 = "11c1ea4414f9b160b0b9f98a3e53f3a2"
        
        # Assert that the calculated hash matches the expected hash
        assert calculated_md5 == expected_md5, f"Expected {expected_md5}, got {calculated_md5}"
    finally:
        # Clean up the temporary file
        os.unlink(temp_file_path)