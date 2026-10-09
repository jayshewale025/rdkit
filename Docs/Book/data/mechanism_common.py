#
# Copyright (C) 2026 RDKit contributors
# This file is part of the RDKit, covered by the BSD license in license.txt.
#
"""Shared scene descriptions for the manually authored mechanism gallery."""

from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Geometry import Point2D, Point3D
from rdkit.Chem.Draw import rdMolDraw2D


@dataclass
class Scene:
  slug: str
  title: str
  description: str
  mols: list
  arrows: list
  offsets: list = field(default_factory=list)
  rotations: list = field(default_factory=list)
  width: int = 720
  height: int = 360


def positioned_molecule(smiles, positions, *, kekulize=False):
  """Preserve explicit H atoms and give each original atom fixed coordinates."""
  parameters = Chem.SmilesParserParams()
  parameters.removeHs = False
  mol = Chem.MolFromSmiles(smiles, parameters)
  if mol is None:
    raise ValueError(f"Invalid gallery SMILES: {smiles!r}")
  if mol.GetNumAtoms() != len(positions):
    raise ValueError(f"{smiles!r} has {mol.GetNumAtoms()} atoms, not {len(positions)}")
  if kekulize:
    Chem.Kekulize(mol, clearAromaticFlags=True)
  conf = Chem.Conformer(mol.GetNumAtoms())
  conf.Set3D(False)
  for index, (x, y) in enumerate(positions):
    conf.SetAtomPosition(index, Point3D(x, y, 0))
  mol.AddConformer(conf, assignId=True)
  return mol


def electron_arrow(start, end, control1, control2, head=None):
  """Use scene-scaled, restrained linework consistently across the examples."""
  if head is None:
    head = rdMolDraw2D.ElectronArrowHead.Full
  arrow = rdMolDraw2D.MechanismArrow(
    start, end, Point2D(*control1), Point2D(*control2), head)
  arrow.colour = (0.12, 0.22, 0.55, 1.0)
  arrow.lineWidth = 0.022
  arrow.headLength = 0.17
  arrow.headWidth = 0.065
  return arrow


def lone_pair(fragment, atom, offset, angle=90):
  anchor = rdMolDraw2D.MechanismAnchor.LonePair(
    fragment, atom, Point2D(*offset), angle=angle)
  anchor.lonePairSeparation = 0.13
  anchor.lonePairRadius = 0.022
  return anchor


def render_scene(scene, backend="svg", *, font_size=22, rotation=0, size=None):
  """Return native SVG text or PNG bytes, plus original-atom pixel positions."""
  if not hasattr(rdMolDraw2D.MolDraw2DSVG, "DrawMechanism"):
    raise RuntimeError("This gallery requires the RDKit fork with DrawMechanism.")
  dimensions = size if size is not None else (scene.width, scene.height)
  if backend == "svg":
    drawer = rdMolDraw2D.MolDraw2DSVG(*dimensions)
  elif backend == "png":
    if not hasattr(rdMolDraw2D, "MolDraw2DCairo"):
      raise RuntimeError("PNG rendering requires RDKit built with Cairo support.")
    drawer = rdMolDraw2D.MolDraw2DCairo(*dimensions)
  else:
    raise ValueError(f"Unknown rendering backend: {backend!r}")
  options = drawer.drawOptions()
  options.fixedFontSize = font_size
  options.fixedBondLength = 65
  options.bondLineWidth = 2
  options.padding = 0.10
  options.rotate = rotation
  options.useBWAtomPalette()
  positions = drawer.DrawMechanism(
    scene.mols, scene.arrows, offsets=scene.offsets, rotations=scene.rotations)
  drawer.FinishDrawing()
  return drawer.GetDrawingText(), positions
