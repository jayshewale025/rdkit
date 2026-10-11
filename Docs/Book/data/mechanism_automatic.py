#
# Copyright (C) 2026 RDKit contributors
# This file is part of the RDKit, covered by the BSD license in license.txt.
#
"""Render endpoint-only mechanisms using automatic native layout and routing."""

import argparse
from pathlib import Path
import re

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from mechanism_common import Scene
from mechanism_gallery import generate_gallery


def _molecule(smiles, symbols, charge=0):
  parameters = Chem.SmilesParserParams()
  parameters.removeHs = False
  mol = Chem.MolFromSmiles(smiles, parameters)
  if mol is None:
    raise ValueError(f"Invalid gallery SMILES: {smiles!r}")
  assert [atom.GetSymbol() for atom in mol.GetAtoms()] == symbols
  assert Chem.GetFormalCharge(mol) == charge
  assert mol.GetNumConformers() == 0
  return mol


def _bond(mol, begin, end, bond_type=Chem.BondType.SINGLE):
  bond = mol.GetBondBetweenAtoms(begin, end)
  assert bond is not None and bond.GetBondType() == bond_type
  return bond.GetIdx()


def scenes() -> list[Scene]:
  """Return fresh, conformer-free scenes with only automatic endpoint arrows."""
  arrow_type = getattr(rdMolDraw2D, "MechanismArrow", None)
  anchor_type = getattr(rdMolDraw2D, "MechanismAnchor", None)
  if not (
      hasattr(rdMolDraw2D.MolDraw2DSVG, "DrawMechanism")
      and hasattr(arrow_type, "autoRoute")
      and hasattr(arrow_type, "autoColour")
      and hasattr(anchor_type, "autoPosition")
  ):
    raise RuntimeError(
      "This gallery requires a build of this RDKit fork with the endpoint-only "
      "MechanismArrow API, not stock RDKit or the older manual-only API.")
  Anchor = rdMolDraw2D.MechanismAnchor
  Arrow = rdMolDraw2D.MechanismArrow
  Head = rdMolDraw2D.ElectronArrowHead

  hydroxide = _molecule("[OH-]", ["O"], -1)
  proton = _molecule("[H+]", ["H"], 1)
  protonation = Scene(
    "protonation", "Protonation: hydroxide donates to H+",
    "One hydroxide lone pair donates to a proton, representing O-H bond "
    "formation in water. Isolated H+ is a simplified endpoint illustration, "
    "not a model of proton solvation. No product is generated.",
    [hydroxide, proton],
    [Arrow(Anchor.LonePair(0, 0), Anchor.Atom(1, 0))],
    width=480, height=260)

  hydroxide = _molecule("[OH-]", ["O"], -1)
  methyl_bromide = _molecule("CBr", ["C", "Br"])
  assert methyl_bromide.GetAtomWithIdx(0).GetTotalNumHs() == 3
  methyl_bromide.GetAtomWithIdx(0).SetProp("_displayLabel", "CH<sub>3</sub>")
  carbon_bromine = _bond(methyl_bromide, 0, 1)
  sn2 = Scene(
    "sn2", "SN2: attack and leaving-group departure",
    "Hydroxide donates to methyl carbon while the C-Br bond pair moves to Br, "
    "describing one concerted substitution giving methanol and bromide. "
    "SN2 requires backside approach; this heuristic 2D layout is not a "
    "transition-state geometry or a guarantee of that approach direction.",
    [hydroxide, methyl_bromide],
    [Arrow(Anchor.LonePair(0, 0), Anchor.Atom(1, 0)),
     Arrow(Anchor.Bond(1, carbon_bromine), Anchor.Atom(1, 1))],
    width=600, height=300)

  aldehyde = _molecule("CC=O", ["C", "C", "O"])
  hydroxide = _molecule("[OH-]", ["O"], -1)
  carbonyl = _bond(aldehyde, 1, 2, Chem.BondType.DOUBLE)
  addition = Scene(
    "carbonyl-addition", "Carbonyl addition: two simultaneous arrows",
    "Hydroxide donates to acetaldehyde's carbonyl carbon while the C=O pi "
    "electron pair moves to oxygen. Together these arrows describe the "
    "tetrahedral alkoxide-forming step of reversible hydration, not a claim "
    "that hydration is favored. The reactant graph remains unchanged.",
    [aldehyde, hydroxide],
    [Arrow(Anchor.LonePair(1, 0), Anchor.Atom(0, 1)),
     Arrow(Anchor.Bond(0, carbonyl), Anchor.Atom(0, 2))],
    width=640, height=320)

  benzene = _molecule("c1ccccc1", ["C"] * 6)
  Chem.Kekulize(benzene, clearAromaticFlags=True)
  assert benzene.GetNumBonds() == 6
  ring, = benzene.GetRingInfo().AtomRings()
  assert len(ring) == 6
  resonance_arrows = []
  # Follow one ring direction, using the actual Kekule bond assignment.
  for index, begin in enumerate(ring):
    middle = ring[(index + 1) % len(ring)]
    end = ring[(index + 2) % len(ring)]
    source = benzene.GetBondBetweenAtoms(begin, middle)
    if source.GetBondType() == Chem.BondType.DOUBLE:
      target = _bond(benzene, middle, end)
      resonance_arrows.append(
        Arrow(Anchor.Bond(0, source.GetIdx()), Anchor.Bond(0, target)))
  assert len(resonance_arrows) == 3
  resonance = Scene(
    "aromatic-resonance", "Benzene resonance: three pi-pair shifts",
    "Three simultaneous arrows move each pi pair to the next single bond "
    "around an actual Kekule form of benzene. This is electron bookkeeping "
    "between resonance contributors, not a reaction or an interconversion "
    "in time. Benzene's electrons are delocalized.",
    [benzene], resonance_arrows, width=560, height=420)

  hydroxide = _molecule("[OH-]", ["O"], -1)
  bromoethane = _molecule("[H]CCBr", ["H", "C", "C", "Br"])
  beta_ch = _bond(bromoethane, 0, 1)
  carbon_carbon = _bond(bromoethane, 1, 2)
  leaving_bond = _bond(bromoethane, 2, 3)
  elimination = Scene(
    "e2-schematic", "E2 electron flow: not an anti-specific depiction",
    "The base accepts an explicit beta-H, its C-H bond pair forms the C=C "
    "pi bond, and the C-Br bond pair moves to Br. This illustrates an E2 "
    "pathway to ethene, water, and bromide, not predicted selectivity over "
    "SN2. RDKit's generated 2D conformer does not guarantee anti-periplanar "
    "H-C-C-Br geometry; use a manually chosen conformer for that depiction.",
    [hydroxide, bromoethane],
    [Arrow(Anchor.LonePair(0, 0), Anchor.Atom(1, 0)),
     Arrow(Anchor.Bond(1, beta_ch), Anchor.Bond(1, carbon_carbon)),
     Arrow(Anchor.Bond(1, leaving_bond), Anchor.Atom(1, 3))])

  hub_mols = [
    _molecule("[H+]", ["H"], 1),
    _molecule("O", ["O"]),
    _molecule("N", ["N"]),
    _molecule("[OH-]", ["O"], -1),
    _molecule("[Cl-]", ["Cl"], -1),
  ]
  hub = Scene(
    "donor-hub-geometry", "Geometry test: alternative donors at a shared hub",
    "Four alternative lone-pair donors target the same H+ to stress fragment "
    "placement, crowding, and routing. The arrows are superimposed alternatives, "
    "NOT a simultaneous chemical step: a proton does not accept four bonds. "
    "This panel is a geometry test only, with no crossing-free guarantee.",
    hub_mols,
    [Arrow(Anchor.LonePair(fragment, 0), Anchor.Atom(0, 0))
     for fragment in range(1, len(hub_mols))],
    width=820, height=480)

  chlorine = _molecule("ClCl", ["Cl", "Cl"])
  chlorine_bond = _bond(chlorine, 0, 1)
  homolysis = Scene(
    "homolysis", "Chlorine homolysis: two single-electron fishhooks",
    "Each fishhook moves one electron from the Cl-Cl bond to a different "
    "chlorine atom, representing homolytic cleavage under irradiation or "
    "heating. The input retains its bond; radical products are not generated.",
    [chlorine],
    [Arrow(Anchor.Bond(0, chlorine_bond), Anchor.Atom(0, 0),
           head=Head.FishhookLeft),
     Arrow(Anchor.Bond(0, chlorine_bond), Anchor.Atom(0, 1),
           head=Head.FishhookRight)],
    width=480, height=280)

  result = [protonation, sn2, addition, resonance, elimination, hub, homolysis]
  for scene in result:
    assert not scene.offsets and not scene.rotations
    assert all(mol.GetNumConformers() == 0 for mol in scene.mols)
    assert scene.arrows and all(
      arrow.autoRoute and arrow.autoColour for arrow in scene.arrows)
    for arrow in scene.arrows:
      for anchor in (arrow.start, arrow.end):
        assert anchor.type != rdMolDraw2D.MechanismAnchorType.Point
        if anchor.type == rdMolDraw2D.MechanismAnchorType.LonePair:
          assert anchor.autoPosition
  return result


_AUTOMATIC_HEADER = """<header>
<h1>Automatic native RDKit mechanism gallery</h1>
<p>Every scene supplies molecules and electron-flow endpoints only: no
conformers, fragment placements, cubic control points, or arrow-colour overrides.
<code>DrawMechanism</code> generates missing coordinates, chooses rigid fragment
placements and curves heuristically, and uses its default arrow palette.
Molecules, lone pairs, and arrows are rendered together, not as SVG overlays.</p>
<p>The API does not predict mechanisms, generate products, or change input
graphs. Intramolecular geometry and stereochemistry are preserved, but layout
is not guaranteed crossing-free or chemically appropriate. Review each figure.
The E2 schematic does not guarantee an anti-periplanar view. The multi-donor
hub is a geometry stress test, NOT a simultaneous chemical step; aromatic
resonance is electron bookkeeping, not a reaction.</p>
</header>"""


def _automatic_index(output):
  """Replace only the shared runner's manual-gallery introduction, not figures."""
  index = output / "index.html"
  page = index.read_text(encoding="utf-8")
  page, replacements = re.subn(
    r"<header>.*?</header>", _AUTOMATIC_HEADER, page, count=1, flags=re.DOTALL)
  if replacements != 1:
    raise RuntimeError("The shared gallery's HTML header could not be found")
  page = page.replace(
    "<title>Native RDKit mechanism gallery</title>",
    "<title>Automatic native RDKit mechanism gallery</title>", 1)
  index.write_text(page, encoding="utf-8")


def main(argv=None):
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("--output", type=Path, default=Path("mechanism-automatic"),
                      help="Output directory (default: mechanism-automatic)")
  parser.add_argument("--png", action="store_true",
                      help="Also render native Cairo PNGs; requires Cairo support")
  args = parser.parse_args(argv)
  try:
    manifest = generate_gallery(scenes(), args.output, png=args.png)
    _automatic_index(args.output)
  except (OSError, ValueError, RuntimeError) as error:
    parser.exit(1, f"Automatic mechanism gallery: {error}\n")
  print(f"Rendered {len(manifest['scenes'])} scenes: {(args.output / 'index.html').resolve()}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
