#!/usr/bin/env python3
"""
VCFStash - A VCF Annotation Cache Manager

Author: Julius Müller, PhD
Organization: GHGA - German Human Genome-Phenome Archive
Date: 16-03-2025
"""

import sys
from pathlib import Path
from src.utils.paths import get_vcfstash_root  # This will set VCFSTASH_ROOT
from src.cli import main

if __name__ == "__main__":
    # Add package root to Python path if needed
    package_root = get_vcfstash_root()
    if package_root != Path.cwd():
        sys.path.insert(0, str(package_root))

    main()
