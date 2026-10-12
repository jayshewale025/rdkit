#
# Copyright (C) 2026 RDKit contributors
# This file is part of the RDKit, covered by the BSD license in license.txt.
#
"""Manually specified resonance and reversible benzaldehyde hydration scenes."""

from rdkit import Chem
from rdkit.Geometry import Point2D as P
from rdkit.Chem.Draw import rdMolDraw2D

from mechanism_common import Scene, positioned_molecule, electron_arrow, lone_pair

Anchor = rdMolDraw2D.MechanismAnchor

_PHENYL = [
  (-1.6, 0.0), (-2.4, 1.386), (-4.0, 1.386),
  (-4.8, 0.0), (-4.0, -1.386), (-2.4, -1.386),
]


def _bond(mol, begin, end, bond_type):
  bond = mol.GetBondBetweenAtoms(begin, end)
  assert bond is not None and bond.GetBondType() == bond_type
  return bond.GetIdx()


def _composition(mols):
  """Count implicit ring hydrogens as well as explicitly retained H atoms."""
  atoms = {}
  for mol in mols:
    for atom in mol.GetAtoms():
      symbol = atom.GetSymbol()
      atoms[symbol] = atoms.get(symbol, 0) + 1
      atoms["H"] = atoms.get("H", 0) + atom.GetTotalNumHs()
  return atoms, sum(Chem.GetFormalCharge(mol) for mol in mols)


def _resonance():
  benzene = positioned_molecule(
    "c1ccccc1",
    [(0.0, 1.6), (1.386, 0.8), (1.386, -0.8),
     (0.0, -1.6), (-1.386, -0.8), (-1.386, 0.8)],
    kekulize=True,
  )
  assert benzene.GetNumAtoms() == benzene.GetNumBonds() == 6
  assert _composition([benzene]) == ({"C": 6, "H": 6}, 0)
  double = Chem.BondType.DOUBLE
  single = Chem.BondType.SINGLE
  arrows = [
    electron_arrow(
      Anchor.Bond(0, _bond(benzene, 1, 2, double), offset=P(-0.25, 0.0)),
      Anchor.Bond(0, _bond(benzene, 2, 3, single), offset=P(-0.06, 0.104)),
      (1.10, -0.43), (1.00, -0.84),
    ),
    electron_arrow(
      Anchor.Bond(0, _bond(benzene, 3, 4, double), offset=P(0.125, 0.217)),
      Anchor.Bond(0, _bond(benzene, 4, 5, single), offset=P(0.12, 0.0)),
      (-0.922, -0.738), (-1.227, -0.446),
    ),
    electron_arrow(
      Anchor.Bond(0, _bond(benzene, 5, 0, double), offset=P(0.125, -0.217)),
      Anchor.Bond(0, _bond(benzene, 0, 1, single), offset=P(-0.06, -0.104)),
      (-0.178, 1.168), (0.227, 1.286),
    ),
  ]
  return Scene(
    slug="aromatic-resonance",
    title="Aromatic resonance: three simultaneous pi shifts",
    description=(
      "Three simultaneous electron-pair arrows on one Kekule contributor of "
      "benzene redistribute the pi bonds to the alternate resonance contributor. "
      "This is resonance-contributor electron bookkeeping, not a reaction or "
      "molecules interconverting in time; benzene has delocalized pi electrons."
    ),
    mols=[benzene], arrows=arrows,
  )


def _addition():
  # Original substrate indices: H 0, carbonyl C 1, O 2, phenyl C 3--8.
  substrate = positioned_molecule(
    "[H]C(=O)c1ccccc1",
    [(0.8, -1.386), (0.0, 0.0), (0.8, 1.386)] + _PHENYL,
    kekulize=True,
  )
  hydroxide = positioned_molecule("[O-][H]", [(0.0, 0.0), (1.3, 0.0)])
  assert [mol.GetNumAtoms() for mol in (substrate, hydroxide)] == [9, 2]
  assert [Chem.GetFormalCharge(mol) for mol in (substrate, hydroxide)] == [0, -1]
  carbonyl = _bond(substrate, 1, 2, Chem.BondType.DOUBLE)
  _bond(hydroxide, 0, 1, Chem.BondType.SINGLE)
  arrows = [
    electron_arrow(
      lone_pair(1, 0, (-0.36, 0.10)),
      Anchor.Atom(0, 1, P(0.14, 0.0)),
      (3.0, 0.35), (1.25, 0.22),
    ),
    electron_arrow(
      Anchor.Bond(0, carbonyl, offset=P(0.15, -0.087)),
      Anchor.Atom(0, 2),
      (1.35, 0.55), (1.60, 1.30),
    ),
  ]
  return Scene(
    slug="carbonyl-addition",
    title="Hydration, step 1: hydroxide addition",
    description=(
      "Manually supplied benzaldehyde and hydroxide: the hydroxide O lone pair "
      "forms a bond to carbonyl C while the C=O pi pair moves to O, giving the "
      "tetrahedral alkoxide shown in step 2. The phenyl ring and explicit "
      "aldehydic H crowd the center without abbreviating its connectivity. "
      "This illustrates one direction of reversible base-catalyzed hydration, "
      "not a prediction that the hydrate is favored."
    ),
    mols=[substrate, hydroxide], arrows=arrows,
    offsets=[P(0.0, 0.0), P(3.85, -0.35)],
  )


def _protonation():
  # Intermediate: H 0, C 1, former carbonyl O- 2, added O 3/H 4, ring 5--10.
  intermediate = positioned_molecule(
    "[H]C([O-])(O[H])c1ccccc1",
    [(0.8, -1.386), (0.0, 0.0), (0.8, 1.386),
     (1.6, 0.0), (2.7, -0.55)] + _PHENYL,
    kekulize=True,
  )
  water = positioned_molecule(
    "[H]O[H]", [(0.0, 0.0), (1.3, 0.0), (1.95, -1.1)])
  assert [mol.GetNumAtoms() for mol in (intermediate, water)] == [11, 3]
  assert intermediate.GetAtomWithIdx(1).GetDegree() == 4
  assert intermediate.GetAtomWithIdx(2).GetFormalCharge() == -1
  _bond(intermediate, 1, 2, Chem.BondType.SINGLE)
  _bond(intermediate, 1, 3, Chem.BondType.SINGLE)
  donor_bond = _bond(water, 0, 1, Chem.BondType.SINGLE)
  arrows = [
    electron_arrow(
      lone_pair(0, 2, (0.05, 0.38), angle=0),
      Anchor.Atom(1, 0),
      (1.2, 2.75), (2.65, 2.90),
    ),
    electron_arrow(
      Anchor.Bond(1, donor_bond, offset=P(0.0, -0.12)),
      Anchor.Atom(1, 1),
      (3.55, 1.10), (4.10, 1.10),
    ),
  ]
  return Scene(
    slug="carbonyl-protonation",
    title="Hydration, step 2: a separate proton-transfer intermediate",
    description=(
      "The manually supplied tetrahedral alkoxide takes an explicit H from "
      "water: the former carbonyl O lone pair points to water H, and that O-H "
      "bond pair returns to water O. This separate proton-transfer step yields "
      "the gem-diol and regenerates hydroxide. No structures are inferred "
      "from the arrows; the reversible process need not favor hydration."
    ),
    mols=[intermediate, water], arrows=arrows,
    offsets=[P(0.0, 0.0), P(3.1, 2.0)],
  )


def _product():
  # Product: H 0, C 1, former carbonyl O 2/new H 3, added O 4/H 5, ring 6--11.
  product = positioned_molecule(
    "[H]C(O[H])(O[H])c1ccccc1",
    [(0.8, -1.386), (0.0, 0.0), (0.8, 1.386), (2.0, 1.8),
     (1.6, 0.0), (2.7, -0.55)] + _PHENYL,
    kekulize=True,
  )
  hydroxide = positioned_molecule("[O-][H]", [(0.0, 0.0), (0.65, -1.1)])
  assert [mol.GetNumAtoms() for mol in (product, hydroxide)] == [12, 2]
  assert product.GetAtomWithIdx(1).GetDegree() == 4
  assert [Chem.GetFormalCharge(mol) for mol in (product, hydroxide)] == [0, -1]
  for begin, end in [(1, 2), (2, 3), (1, 4), (4, 5)]:
    _bond(product, begin, end, Chem.BondType.SINGLE)
  return Scene(
    slug="carbonyl-product",
    title="Hydration: manually supplied product and regenerated base",
    description=(
      "Manually supplied benzaldehyde hydrate (the gem-diol, Ph-CH(OH)2) and "
      "regenerated hydroxide after step 2. The aldehydic H and both O-H bonds "
      "remain explicit. This is a stoichiometrically balanced electron-flow "
      "illustration of reversible hydration, not a claim of an isolated or "
      "thermodynamically favored product. No reaction inference was used."
    ),
    mols=[product, hydroxide], arrows=[],
    offsets=[P(0.0, 0.0), P(4.4, 2.0)],
  )


def scenes() -> list[Scene]:
  """Return four native scenes with explicit graphs, coordinates and anchors."""
  resonance, addition, protonation, product = (
    _resonance(), _addition(), _protonation(), _product())
  assert _composition(addition.mols) == _composition(protonation.mols[:1])
  assert _composition(protonation.mols) == _composition(product.mols)
  assert _composition(addition.mols) == ({"C": 7, "H": 7, "O": 2}, -1)
  assert _composition(product.mols) == ({"C": 7, "H": 9, "O": 3}, -1)
  return [resonance, addition, protonation, product]
