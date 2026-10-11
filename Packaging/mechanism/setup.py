"""Metadata for the staged, already-compiled fork; not a C++ build backend."""

import json
from pathlib import Path

from setuptools import Distribution, find_packages, setup


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True


root = Path(__file__).resolve().parent
info = json.loads((root / "rdkit_mechanism" / "build-info.json").read_text())
packages = find_packages()
data = {}
for package in packages:
    directory = root.joinpath(*package.split("."))
    data[package] = [
        str(path.relative_to(directory))
        for path in directory.rglob("*")
        if path.is_file() and "__pycache__" not in path.parts
    ]

setup(
    name="rdkit",
    version=info["version"],
    description="RDKit source fork with native endpoint-only mechanism drawing",
    url="https://github.com/jayshewale025/rdkit",
    license="BSD-3-Clause",
    license_files=["license.txt"],
    packages=packages,
    package_data=data,
    include_package_data=False,
    python_requires=">=3.12",
    install_requires=["numpy>=2,<3", "Pillow>=10"],
    entry_points={"console_scripts": ["rdkit-mechanism=rdkit_mechanism.__main__:main"]},
    distclass=BinaryDistribution,
    zip_safe=False,
)
