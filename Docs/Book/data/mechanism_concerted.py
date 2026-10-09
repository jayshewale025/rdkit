#
# Copyright (C) 2026 RDKit contributors
# This file is part of the RDKit, covered by the BSD license in license.txt.
#
"""Explicitly authored, concerted two-electron mechanisms for the gallery."""

from rdkit.Geometry import Point2D as P
from rdkit.Chem.Draw.rdMolDraw2D import MechanismAnchor as Anchor

from mechanism_common import Scene, positioned_molecule, electron_arrow, lone_pair


def _sn2():
  hydroxide = positioned_molecule("[OH-]", [(-2.65, 0)])
  methyl_bromide = positioned_molecule("CBr", [(0, 0), (1.65, 0)])
  # A visible methyl label avoids a bare line end at the site of backside attack.
  methyl_bromide.GetAtomWithIdx(0).SetProp("_displayLabel", "CH<sub>3</sub>")
  arrows = [
    electron_arrow(lone_pair(0, 0, (0.35, 0)), Anchor.Atom(1, 0),
                   (-1.65, 0.5), (-0.75, 0.65)),
    electron_arrow(Anchor.Bond(1, 0, offset=P(0, 0.04)), Anchor.Atom(1, 1),
                   (0.95, 0.95), (2.2, 0.8)),
  ]
  return Scene(
    "sn2", "SN2: backside displacement",
    "Hydroxide donates a lone pair to the methyl carbon from the side opposite "
    "bromine while the C–Br bond pair moves to Br. Both arrows describe one "
    "concerted step, forming methanol and bromide; no intermediate is implied.",
    [hydroxide, methyl_bromide], arrows, width=560, height=280)


def _e2():
  hydroxide = positioned_molecule("[OH-]", [(-3.25, 1.32)])
  # H(0)–Cbeta(1)–Calpha(3)–Br(5) is planar with a 180-degree dihedral.
  # The explicit beta-H and Br are on opposite sides of the central C–C bond.
  bromobutane = positioned_molecule(
    "[H]C(C)C(C)Br",
    [(-0.78, 1.32), (0, 0), (-0.78, -1.32),
     (1.5, 0), (2.28, 1.32), (2.28, -1.32)])
  arrows = [
    electron_arrow(lone_pair(0, 0, (0.35, 0)), Anchor.Atom(1, 0),
                   (-2.2, 1.6), (-1.35, 2.1)),
    electron_arrow(Anchor.Bond(1, 0), Anchor.Bond(1, 2, offset=P(0, 0.12)),
                   (-0.03, 1.05), (0.8, 0.9)),
    electron_arrow(Anchor.Bond(1, 4), Anchor.Atom(1, 5),
                   (2.8, -0.28), (3.1, -1.62)),
  ]
  return Scene(
    "e2", "E2: anti-periplanar elimination",
    "A planar side view of H–Cβ–Cα–Br with a 180° dihedral, not a syn "
    "arrangement. Hydroxide abstracts the explicit β-H, its C–H bond pair "
    "forms the Cβ=Cα π bond, and the Cα–Br bond pair moves to Br. The three "
    "arrows are simultaneous: but-2-ene, water, and bromide form in one step.",
    [hydroxide, bromobutane], arrows)


def _proton_transfer():
  ammonia = positioned_molecule("N", [(-3.05, 0)])
  # The transferable H is original atom 4, not an inferred implicit hydrogen.
  acetic_acid = positioned_molecule(
    "CC(=O)O[H]", [(3.6, 0), (2.3, 0.75), (2.3, 2.25), (1, 0), (-0.5, 0)])
  arrows = [
    electron_arrow(lone_pair(0, 0, (0, 0.4), angle=0), Anchor.Atom(1, 4),
                   (-2.9, 1.25), (-1.2, 1.25)),
    electron_arrow(Anchor.Bond(1, 3, offset=P(0, -0.04)), Anchor.Atom(1, 3),
                   (0.1, -0.8), (1.7, -0.75)),
  ]
  return Scene(
    "proton-transfer", "Proton transfer: follow the explicit H",
    "Ammonia's explicitly drawn nitrogen lone pair attacks the explicit acidic "
    "hydrogen of acetic acid. At the same time, the O–H bond pair returns to "
    "that oxygen, forming ammonium and acetate without a free-proton "
    "intermediate.",
    [ammonia, acetic_acid], arrows)


def scenes() -> list[Scene]:
  """Return fresh molecules and arrows, preserving original atom/bond indices."""
  return [_sn2(), _e2(), _proton_transfer()]
