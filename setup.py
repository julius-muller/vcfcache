from setuptools import setup, find_packages

setup(
    name="vepstash",
    version="1.0.0",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "vepstash=src.cli:main",
        ],
    },
    install_requires=[
        "pandas",
        "pyarrow",
        "pysam",
    ],
)