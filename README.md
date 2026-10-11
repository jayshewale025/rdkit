# RDKit
[![Azure build Status](https://dev.azure.com/rdkit-builds/RDKit/_apis/build/status/rdkit.rdkit?branchName=master)](https://dev.azure.com/rdkit-builds/RDKit/_build/latest?definitionId=1&branchName=master)
[![DOI](https://zenodo.org/badge/10009991.svg)](https://zenodo.org/badge/latestdoi/10009991)


## What is it?

The [RDKit](https://www.rdkit.org) is a collection of cheminformatics and machine-learning software written in C++ and Python.

  * [BSD license](https://github.com/rdkit/rdkit/blob/master/license.txt) - a business friendly license for open source
  * Core data structures and algorithms in C++
  * [Python 3.x wrapper](https://www.rdkit.org/docs/GettingStartedInPython.html) generated using Boost.Python
  * Java and C# wrappers generated with SWIG
  * JavaScript (generated with emscripten) and CFFI wrappers around important functionality
  * 2D and 3D molecular operations
  * [Descriptor](https://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-descriptors) and [Fingerprint](http://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-fingerprints) generation for machine learning
  * Molecular database [cartridge](https://www.rdkit.org/docs/Cartridge.html) for PostgreSQL supporting substructure and similarity searches as well as many descriptor calculators
  * Cheminformatics nodes for [KNIME](https://www.knime.com/rdkit)
  * [Contrib](https://github.com/rdkit/rdkit/tree/master/Contrib) folder with useful community-contributed software harnessing the power of the RDKit


## Automatic electron-flow diagrams (this source fork)

This fork adds native **endpoint-only mechanism drawing** to `MolDraw2D`.
Supply molecules and the electron-flow endpoints; RDKit generates missing 2D
coordinates, chooses fragment placement and curved arrows heuristically, and
uses a deterministic arrow-colour palette. Molecules, lone pairs, and arrows
are rendered together in SVG or, with Cairo support, PNG—not as SVG overlays.

**Installing stock RDKit with pip or conda does not provide this API.**
Use this fork's wheel or build its `feature/mechanism-drawing` branch.
The [wheel installation and quickstart guide](Packaging/mechanism/README.md)
covers Windows x64 / Python 3.12, ARM64 Linux / Python 3.13 CI artifacts,
isolated environments, PNG/SVG output, and the one-page textbook example.
For native development, follow the existing
[RDKit build instructions](Docs/Book/Install.md) using this checkout.

### Start with endpoints, not coordinates

This complete example draws electron-pair donation from hydroxide to a proton.
The isolated `H+` is a simplified protonation illustration, not a representation
of proton solvation. No conformers, offsets, rotations, or control points are
supplied:

```python
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

Anchor = rdMolDraw2D.MechanismAnchor
Arrow = rdMolDraw2D.MechanismArrow
hydroxide = Chem.MolFromSmiles("[OH-]")
proton = Chem.MolFromSmiles("[H+]")
assert hydroxide.GetAtomWithIdx(0).GetAtomicNum() == 8
assert proton.GetAtomWithIdx(0).GetAtomicNum() == 1

arrow = Arrow(Anchor.LonePair(0, 0), Anchor.Atom(1, 0))
drawer = rdMolDraw2D.MolDraw2DSVG(480, 240)
atom_pixels = drawer.DrawMechanism([hydroxide, proton], [arrow])
drawer.FinishDrawing()
Path("protonation.svg").write_text(drawer.GetDrawingText(), encoding="utf-8")
```

Carbonyl addition needs **two** simultaneous electron-pair arrows: donation
to carbonyl carbon and movement of the carbonyl pi pair to oxygen. Bond indices
are obtained from the actual molecule:

```python
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

Anchor = rdMolDraw2D.MechanismAnchor
Arrow = rdMolDraw2D.MechanismArrow
aldehyde = Chem.MolFromSmiles("CC=O")
hydroxide = Chem.MolFromSmiles("[OH-]")
assert [a.GetSymbol() for a in aldehyde.GetAtoms()] == ["C", "C", "O"]
assert hydroxide.GetAtomWithIdx(0).GetAtomicNum() == 8
carbonyl = aldehyde.GetBondBetweenAtoms(1, 2)
assert carbonyl.GetBondType() == Chem.BondType.DOUBLE
arrows = [
    Arrow(Anchor.LonePair(1, 0), Anchor.Atom(0, 1)),
    Arrow(Anchor.Bond(0, carbonyl.GetIdx()), Anchor.Atom(0, 2)),
]
drawer = rdMolDraw2D.MolDraw2DSVG(640, 320)
drawer.DrawMechanism([aldehyde, hydroxide], arrows)
drawer.FinishDrawing()
Path("carbonyl-addition.svg").write_text(drawer.GetDrawingText(), encoding="utf-8")
```

### Automatic defaults and manual escape hatches

* `MechanismArrow(start, end, head=ElectronArrowHead.Full)` sets `autoRoute=True`
  and `autoColour=True`. `FishhookLeft` and `FishhookRight` are available for
  one-electron arrows. With automatic routing, stored `control1` and `control2`
  are ignored.
* `MechanismAnchor.LonePair(fragment, atomIndex, offset=None)` sets
  `autoPosition=True`: RDKit chooses where to draw the requested pair of dots.
  It does not infer whether an atom has a lone pair. An explicit offset selects
  manual dot placement.
* `DrawMechanism(..., autoLayout=True)` arranges fragments only for a nonempty
  set of all-automatic arrows, with no absolute `Point` anchors and no supplied
  fragment `offsets` or `rotations`. A manual arrow, absolute point, or explicit
  fragment placement opts out. `autoLayout=False` retains fragment placement
  while still choosing curves for arrows whose `autoRoute` is true.
* The original `MechanismArrow(start, end, control1, control2, head=...)`
  constructor remains available for manual cubic curves; it defaults both
  automatic arrow flags to false. Assigning an arrow's `.colour` in Python
  disables `autoColour` without disabling automatic routing.

Automatic fragment layout uses **rigid rotations and translations**, preserving
intramolecular geometry and stereochemistry. Missing conformers are generated
on drawing copies. The API does not change input graphs, infer mechanisms, or
generate products. Its layout and routing are heuristic: **there is no
crossing-free or chemically correct placement guarantee**. Review every figure;
use explicit geometry for a particular stereochemical view, such as an
anti-periplanar E2 depiction.

### Runnable galleries

From the repository root, with the fork available in your Python environment:

```console
python Docs/Book/data/mechanism_automatic.py --output mechanism-automatic
python Docs/Book/data/mechanism_automatic.py --output mechanism-automatic --png
```

Open `mechanism-automatic/index.html`. The
[automatic gallery source](Docs/Book/data/mechanism_automatic.py) supplies only
molecules and endpoints for seven scenes, including resonance, fishhooks, and
a clearly labelled multi-donor **geometry test**, not a simultaneous reaction.
`--png` additionally requires Cairo. See the
[Cookbook's mechanism section](Docs/Book/Cookbook.rst#drawing-user-specified-electron-flow-arrows)
for API details and the
[existing manual gallery](Docs/Book/data/mechanism_gallery.py) for advanced
conformer, placement, and curve controls.


## Installation and getting started

The package commands below install standard RDKit, not this fork's mechanism API.

If you are working in Python and using conda (our recommendation), installation is super easy:

```shell-session
$ conda install -c conda-forge rdkit
```

You can then take a look at our [Getting Started in Python](https://rdkit.org/docs/GettingStartedInPython.html) guide.

More detailed installation instructions are available in [Docs/Book/Install.md](https://github.com/rdkit/rdkit/blob/master/Docs/Book/Install.md).

## Documentation
Available on the [RDKit page](https://www.rdkit.org/docs/index.html)
and in the [Docs](https://github.com/rdkit/rdkit/tree/master/Docs) folder on GitHub

The [RDKit blog](https://greglandrum.github.io/rdkit-blog/) often has useful tips and tricks.

## Support and Community

If you have questions, comments, or suggestions, the best places for those are:

  * [GitHub discussions](https://github.com/rdkit/rdkit/discussions)
  * The [mailing list](https://sourceforge.net/p/rdkit/mailman/)

If you've found a bug or would like to request a feature, please [create an issue](https://github.com/rdkit/rdkit/issues)

We also have a [LinkedIn group](https://www.linkedin.com/groups/RDKit-8192558/about)

We have a yearly user group meeting (the UGM) where members of the community do presentations and lightning talks on things they've done with the RDKit. Materials from past UGMs, which can quite useful, are also online:
  * [2012 UGM, London](http://www.rdkit.org/UGM/2012/)
  * [2013 UGM, Hinxton](https://github.com/rdkit/UGM_2013)
  * [2014 UGM, Darmstadt](https://github.com/rdkit/UGM_2014)
  * [2015 UGM, Zurich](https://github.com/rdkit/UGM_2015)
  * [2016 UGM, Basel](https://github.com/rdkit/UGM_2016)
  * [2017 UGM, Berlin](https://github.com/rdkit/UGM_2017)
  * [2018 UGM, Cambridge](https://github.com/rdkit/UGM_2018)
  * [2019 UGM, Hamburg](https://github.com/rdkit/UGM_2019)
  * [2020 UGM, virtual](https://github.com/rdkit/UGM_2020)
  * [2021 UGM, virtual](https://github.com/rdkit/UGM_2021)
  * [2022 UGM, Berlin](https://github.com/rdkit/UGM_2022)
  * [2023 UGM, Mainz](https://github.com/rdkit/UGM_2023)
  * [2024 UGM, Zurich](https://github.com/rdkit/UGM_2024)

## License

Code released under the [BSD license](https://github.com/rdkit/rdkit/blob/master/license.txt).
