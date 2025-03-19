#!/usr/bin/env python3
"""
VEPStash - A VEP Annotation Cache Manager

Author: Julius Müller, PhD
Organization: GHGA - German Human Genome-Phenome Archive
Date: 16-03-2025
"""

import sys
from pathlib import Path
from src.cli import main

if __name__ == "__main__":
    if not Path(sys.argv[0]).parent.resolve() == Path.cwd():
        sys.path.insert(0, str(Path(sys.argv[0]).parent.resolve()))
    main()