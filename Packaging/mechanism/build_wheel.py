"""Assemble and repair a wheel from a complete, matching CMake Python build."""

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REQUIRED = (
    "rdBase", "rdchem", "rdmolfiles", "rdmolops", "rdDepictor", "rdMolDraw2D",
    "rdMolDescriptors", "rdDistGeom", "rdForceFieldHelpers", "rdMolStandardize",
    "rdChemReactions", "rdinchi",
)


def run(*args, **kwargs):
    subprocess.run([str(arg) for arg in args], check=True, **kwargs)


def conda_licenses(prefix, destination):
    records = []
    for metadata in sorted((prefix / "conda-meta").glob("*.json")):
        record = json.loads(metadata.read_text(encoding="utf-8"))
        source = Path(record.get("link", {}).get("source", "")) / "info" / "licenses"
        if source.is_dir():
            shutil.copytree(source, destination / record["name"])
            records.append({
                "name": record["name"], "version": record["version"],
                "license": record.get("license", ""),
                "package_url": record.get("url", ""),
            })
    if not records:
        raise RuntimeError("No dependency licenses found in the conda package cache")
    (destination / "conda-packages.json").write_text(
        json.dumps(records, indent=2), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--runtime-prefix", type=Path)
    parser.add_argument("--licenses-dir", type=Path)
    parser.add_argument("--version", default="2026.9.1+mechanism.1")
    parser.add_argument("--platform", default="manylinux_2_28_aarch64")
    args = parser.parse_args()
    build, output = args.build_dir.resolve(), args.output.resolve()
    if sys.platform == "win32" and args.runtime_prefix is None:
        parser.error("--runtime-prefix is required on Windows")
    if sys.platform == "linux" and args.licenses_dir is None:
        parser.error("--licenses-dir is required on Linux")
    output.mkdir(parents=True, exist_ok=True)
    suffix = ".pyd" if sys.platform == "win32" else ".so"
    extensions = sorted((build / "rdkit").rglob("*" + suffix))
    names = {path.name.split(".")[0] for path in extensions}
    missing = set(REQUIRED) - names
    if missing:
        raise RuntimeError(f"Incomplete CMake build; missing: {sorted(missing)}")
    cache = (build / "CMakeCache.txt").read_text(encoding="utf-8")
    if "RDK_BUILD_BOOST_PYTHON_WRAPPERS:BOOL=ON" not in cache:
        raise RuntimeError("This packaging recipe currently requires Boost.Python")
    options = dict(re.findall(r"^(RDK_[A-Z0-9_]+):BOOL=(ON|OFF)$", cache, re.M))
    revision = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
    dirty = bool(subprocess.check_output(
        ["git", "-C", str(ROOT), "status", "--porcelain"], text=True).strip())

    with tempfile.TemporaryDirectory(prefix="rdkit-mechanism-wheel-") as temporary:
        stage = Path(temporary)
        package = stage / "rdkit"
        ignored = shutil.ignore_patterns(
            "__pycache__", "CMakeFiles", "*.pyc", "*.pyd", "*.so", "*.cmake",
            "CMakeLists.txt", "CTestTestfile.cmake", "test_data", "testData",
            "UnitTest*", "test*.py", "conftest.py", "pytest.ini",
        )
        shutil.copytree(ROOT / "rdkit", package, ignore=ignored)
        shutil.copytree(ROOT / "Data", package / "Data",
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        for extension in extensions:
            destination = package / extension.relative_to(build / "rdkit")
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(extension, destination)
        shutil.copy2(ROOT / "External" / "INCHI-API" / "python" / "inchi.py",
                     package / "Chem" / "inchi.py")
        shutil.copy2(HERE / "RDPaths.py", package / "RDPaths.py")
        for name in ("setup.py", "pyproject.toml"):
            shutil.copy2(HERE / name, stage / name)
        shutil.copy2(ROOT / "license.txt", stage / "license.txt")
        helper = stage / "rdkit_mechanism"
        shutil.copytree(HERE / "rdkit_mechanism", helper, ignore=ignored)
        examples = helper / "examples"
        examples.mkdir()
        for example in (ROOT / "Docs" / "Book" / "data").glob("mechanism_*.py"):
            shutil.copy2(example, examples / example.name)
        notices = helper / "licenses"
        notices.mkdir()
        shutil.copy2(ROOT / "license.txt", notices / "RDKit.txt")
        if args.runtime_prefix is not None:
            conda_licenses(args.runtime_prefix.resolve(), notices)
        if args.licenses_dir is not None:
            shutil.copytree(args.licenses_dir, notices / "native-dependencies")
        for directory in (ROOT / "External", ROOT / "Code"):
            for notice in directory.rglob("*"):
                if notice.is_file() and (
                    notice.name.lower().startswith(("license", "copying"))
                    or notice.name.lower() == "copyright"
                ):
                    relative = notice.relative_to(ROOT)
                    destination = notices / "source-components" / relative
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(notice, destination)
        (helper / "build-info.json").write_text(json.dumps({
            "version": args.version,
            "source_revision": revision, "uncommitted_changes": dirty,
            "python": sys.version.split()[0], "native_modules": sorted(names),
            "cmake_options": options,
        }, indent=2), encoding="utf-8")
        raw = stage / "unrepaired"
        run(sys.executable, "-m", "pip", "wheel", "--no-deps", "--no-build-isolation",
            "--wheel-dir", raw, stage)
        wheels = list(raw.glob("*.whl"))
        if len(wheels) != 1:
            raise RuntimeError(f"Expected one wheel, found {len(wheels)}")
        if sys.platform == "win32":
            run(sys.executable, "-m", "delvewheel", "repair", wheels[0],
                "--add-path", args.runtime_prefix.resolve() / "Library" / "bin",
                "--wheel-dir", output)
        elif sys.platform == "linux":
            run(sys.executable, "-m", "auditwheel", "repair", wheels[0],
                "--plat", args.platform, "--wheel-dir", output)
        else:
            parser.error("Only Windows and Linux packaging are currently supported")
    repaired = list(output.glob(f"rdkit-{args.version}-*.whl"))
    if not repaired:
        raise RuntimeError("Wheel repair did not produce the expected artifact")
    for wheel in repaired:
        with zipfile.ZipFile(wheel) as archive:
            paths = archive.namelist()
            assert "rdkit_mechanism/examples/mechanism_question6.py" in paths
            assert "rdkit/Data/Fonts/Telex-Regular.ttf" in paths
        print(f"Ready: {wheel} ({wheel.stat().st_size / 1024**2:.1f} MiB)")


if __name__ == "__main__":
    main()
