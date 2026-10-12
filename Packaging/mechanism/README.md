# Install the mechanism-drawing fork

These wheels contain the native RDKit fork, its Python API, runtime libraries,
data/fonts, and runnable examples. They are **not** plugins for stock RDKit.
The distribution is named `rdkit`; install it in a new virtual environment and
do not subsequently install stock `rdkit` or `rdkit-pypi` into that environment.

## Supported packages

| Wheel suffix | Target |
| --- | --- |
| `cp312-cp312-win_amd64.whl` | 64-bit Windows, standard CPython 3.12 |
| `cp313-cp313-manylinux_2_28_aarch64.whl` | 64-bit ARM Linux, standard CPython 3.13, glibc 2.28+ |

Python patch versions share the ABI: the ARM package is appropriate for
Python **3.13.5**. It is not for 32-bit Raspberry Pi OS, x86 Linux, musl/Alpine,
or free-threaded Python. On a Pi, check `uname -m` (must say `aarch64`) and
`python3 --version` before installing. Raspberry Pi OS Bookworm/Trixie meet
the glibc baseline, but the Python version must match the wheel.

The ARM workflow builds and exercises the package on native ARM64 Linux, not
on physical Pi hardware. A successful CI artifact is the release gate; merely
adding the workflow is not evidence that an ARM wheel was built.

The fork currently has InChI, SVG/Cairo drawing, core descriptors, conformer
generation and force fields. Optional CoordGen, MaeParser, ChemDraw,
PubChemShape, Eigen-based 3D descriptors, and Boost serialization/iostreams
are disabled in this build recipe. Exact options and module names are in
`rdkit_mechanism/build-info.json` inside each wheel.

## Windows: install and draw

Download/copy the Windows wheel into the current directory. These commands use
the Python launcher; substitute the full path to your Python 3.12 executable
if `py` is not installed.

```powershell
py -3.12 -m venv .venv
$python = ".\.venv\Scripts\python.exe"
$wheel = Get-ChildItem -File "rdkit-2026.9.1+mechanism.1-cp312-cp312-win_amd64.whl"
& $python -m pip install --upgrade pip
& $python -m pip install $wheel.FullName
& $python -m rdkit_mechanism doctor
& $python -m rdkit_mechanism demo --output drawings
& $python -m rdkit_mechanism question6 --compact --output question6
```

Open `drawings\carbonyl-addition.svg` or `question6\textbook.html`.
No C++ compiler, conda activation, `RDBASE`, `PYTHONPATH`, or DLL-path changes
are needed for an installed wheel.

## Raspberry Pi OS: install and draw

In the fork's **Actions** tab, open **Mechanism ARM64 wheel**. Relevant native
code, Python, data, packaging, and example changes trigger it for pull requests
targeting `master` and pushes to `master`. You can also use **Run workflow**
to build `master` manually. After the job succeeds, download and unzip
`rdkit-mechanism-cp313-aarch64`. This is a workflow artifact, not a PyPI release.

Run these commands **on the Pi**, in the directory containing the ARM wheel:

```bash
python3.13 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install ./rdkit-2026.9.1+mechanism.1-cp313-cp313-manylinux_2_28_aarch64.whl
.venv/bin/python -m rdkit_mechanism doctor
.venv/bin/python -m rdkit_mechanism demo --output drawings
.venv/bin/python -m rdkit_mechanism question6 --compact --output question6
```

If venv creation reports missing `ensurepip`, install the venv component for
the same Python interpreter. Do not use `sudo pip` or
`--break-system-packages`. Python 3.11, commonly installed with Bookworm, needs
a separate CPython 3.11 wheel or a separately installed Python 3.13.

## Use it in your own Python code

```python
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D as draw

A, Arrow = draw.MechanismAnchor, draw.MechanismArrow
aldehyde = Chem.MolFromSmiles("CC=O")
hydroxide = Chem.MolFromSmiles("[OH-]")
carbonyl = aldehyde.GetBondBetweenAtoms(1, 2).GetIdx()
arrows = [
    Arrow(A.LonePair(1, 0), A.Atom(0, 1)),
    Arrow(A.Bond(0, carbonyl), A.Atom(0, 2)),
]

drawer = draw.MolDraw2DSVG(640, 320)
drawer.DrawMechanism([aldehyde, hydroxide], arrows)
drawer.FinishDrawing()
Path("addition.svg").write_text(drawer.GetDrawingText(), encoding="utf-8")
```

`fragment` numbers are positions in the molecule list; atom/bond indices belong
to the corresponding original molecule. No coordinates or controls are needed.
Use `MolDraw2DCairo` and `write_bytes()` for PNG. The seven-scene gallery is:

```console
python -m rdkit_mechanism gallery --png --output gallery
```

The installed examples are ordinary editable Python files:

```console
python -c "import pathlib, rdkit_mechanism; print(pathlib.Path(rdkit_mechanism.__file__).parent / 'examples')"
```

Copy examples into your own project rather than editing site-packages.
The source-controlled versions are in `Docs/Book/data/mechanism_*.py`.
`mechanism_question6.py --compact` makes the one-page scheme; omit `--compact`
for the detailed version. Pillow composes the page and straight reaction arrows;
RDKit natively renders molecules and curved electron arrows.

## Limits and troubleshooting

- This is a drawing API, not a reaction predictor. Supply the intermediate
  structures and chemically meaningful electron movements yourself.
- Fragment placement is rigid. Internal bonds are not flexed to accommodate a
  mechanism, and crossing-free routing is not guaranteed.
- Explicit offsets and manual control points remain available. Multiple-bond
  automatic tails attach to the last visible stroke on their departure side.
- If `MechanismArrow` is missing, run `doctor` and check `rdkit_location`:
  stock RDKit or a source checkout may be shadowing the wheel. Remove stale
  `PYTHONPATH`/`RDBASE` settings and run outside the source checkout.
- "Not a supported wheel": verify architecture, Python minor version, and
  whether the interpreter is free-threaded. Updating pip cannot make a Windows
  wheel run on Linux.
- Dependencies are installed from PyPI by pip. Internet-free installation
  additionally requires matching NumPy and Pillow wheels downloaded beforehand.
- Licenses and dependency inventories ship in `rdkit_mechanism/licenses`.

## Reproduce the package

Packaging consumes a complete CMake build; it does not compile missing modules
or quietly borrow binaries from stock RDKit.

```powershell
python -m pip install -r Packaging\mechanism\requirements-build.txt
python Packaging\mechanism\build_wheel.py `
  --build-dir build-mechanism --runtime-prefix C:\path\to\build-environment `
  --output build-mechanism\wheelhouse
```

Use the same Python minor version and architecture that compiled the native
extensions. Build every enabled Python target first, not just `rdMolDraw2D`.
Windows DLLs are bundled and load paths repaired with delvewheel. The runtime
prefix must retain its conda package-cache license records.

The ARM CI uses `build_linux_arm.sh` inside the manylinux ARM64 image, verifies
the Boost source checksum, bundles shared libraries with auditwheel, and tests
an installed wheel after clearing build-library paths. Neither workflow nor
builder uploads packages to PyPI or creates releases.
