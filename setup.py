"""Setup script used for package installation"""

import os
from setuptools import setup, find_packages


def read(filename):
    return open(os.path.join(os.path.dirname(__file__), filename)).read()


setup(
    name="nds",
    version="1.0",
    description="Path integral simulations for nonlinear dissipative systems.",
    packages=find_packages(include=["nds", "nds.*"]),
    author="Oli Bridge",
    author_email="ob344@cam.ac.uk",
)
