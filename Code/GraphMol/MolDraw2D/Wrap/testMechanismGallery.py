#
# Copyright (C) 2026 RDKit contributors
# This file is part of the RDKit, covered by the BSD license in license.txt.
#
import contextlib
from collections import Counter
from dataclasses import replace
import io
import json
import math
from pathlib import Path
import struct
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

import numpy as np
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from testMechanismDrawing import SVG, pathPoints, snapshot, withClass

DATA = Path(__file__).resolve().parents[4] / 'Docs' / 'Book' / 'data'
sys.path.insert(0, str(DATA))
from mechanism_common import positioned_molecule, render_scene
from mechanism_gallery import build_scenes, generate_gallery, main


def curve_extrema(points):
  result = [points[0], points[-1]]
  parameters = set()
  for axis in ('x', 'y'):
    p0, p1, p2, p3 = [getattr(point, axis) for point in points]
    for root in np.roots([-p0 + 3*p1 - 3*p2 + p3,
                          2*(p0 - 2*p1 + p2), p1 - p0]):
      if abs(root.imag) < 1e-10 and 0 < root.real < 1:
        parameters.add(float(root.real))
  for t in parameters:
    u = 1 - t
    result.append(points[0] * (u**3) + points[1] * (3*u*u*t) +
                  points[2] * (3*u*t*t) + points[3] * (t**3))
  return result


class TestMechanismGallery(unittest.TestCase):

  def setUp(self):
    self.scenes = build_scenes()

  def testGalleryCoverageAndExplicitHydrogens(self):
    expected = {
      'sn2': 2, 'e2': 3, 'proton-transfer': 2, 'aromatic-resonance': 3,
      'carbonyl-addition': 2, 'carbonyl-protonation': 2, 'carbonyl-product': 0,
      'homolysis': 2,
    }
    self.assertEqual({scene.slug: len(scene.arrows) for scene in self.scenes}, expected)
    for scene in self.scenes:
      with self.subTest(scene=scene.slug):
        self.assertTrue(scene.title and scene.description)
        self.assertTrue(all(mol.GetNumConformers() == 1 for mol in scene.mols))
        self.assertTrue(all(not mol.GetConformer().Is3D() for mol in scene.mols))
        if scene.slug in ('e2', 'proton-transfer', 'carbonyl-addition',
                          'carbonyl-protonation'):
          self.assertTrue(any(atom.GetAtomicNum() == 1
                              for mol in scene.mols for atom in mol.GetAtoms()))
    resonance = next(scene for scene in self.scenes if scene.slug == 'aromatic-resonance')
    for arrow in resonance.arrows:
      self.assertEqual(arrow.start.type, rdMolDraw2D.MechanismAnchorType.Bond)
      self.assertEqual(arrow.end.type, rdMolDraw2D.MechanismAnchorType.Bond)
      self.assertEqual(resonance.mols[arrow.start.fragment].GetBondWithIdx(
        arrow.start.index).GetBondType(), Chem.BondType.DOUBLE)
      self.assertEqual(resonance.mols[arrow.end.fragment].GetBondWithIdx(
        arrow.end.index).GetBondType(), Chem.BondType.SINGLE)

  def testNativeGeometryAndNonmutationAcrossLayouts(self):
    for scene in self.scenes:
      original = [snapshot(mol) for mol in scene.mols]
      for font, rotation, size in ((22, 0, None), (28, 0, None),
                                   (18, 35, None), (-1, 0, (-1, -1))):
        with self.subTest(scene=scene.slug, font=font, rotation=rotation, size=size):
          svg, positions = render_scene(scene, font_size=font, rotation=rotation, size=size)
          root = ET.fromstring(svg)
          width, height = float(root.get('width').removesuffix('px')), float(
            root.get('height').removesuffix('px'))
          self.assertEqual([len(fragment) for fragment in positions],
                           [mol.GetNumAtoms() for mol in scene.mols])
          self.assertEqual([snapshot(mol) for mol in scene.mols], original)
          shafts = withClass(root, 'mechanism-shaft')
          heads = withClass(root, 'mechanism-head')
          self.assertEqual(len(shafts), len(scene.arrows))
          self.assertEqual(len(heads), len(scene.arrows))
          for index, (shaft, head) in enumerate(zip(shafts, heads)):
            self.assertIn(f'mechanism-arrow-{index}', shaft.get('class'))
            self.assertIn(f'mechanism-arrow-{index}', head.get('class'))
            self.assertEqual(shaft.get('d').count('C'), 1)
            self.assertNotIn('L', shaft.get('d'))
            points = pathPoints(shaft)
            self.assertEqual(len(points), 4)
            self.assertIn('fill:none', shaft.get('style'))
            for point in curve_extrema(points) + pathPoints(head):
              self.assertTrue(math.isfinite(point.x) and math.isfinite(point.y))
              self.assertGreaterEqual(point.x, 0.5)
              self.assertLessEqual(point.x, width - 0.5)
              self.assertGreaterEqual(point.y, 0.5)
              self.assertLessEqual(point.y, height - 0.5)
          for dot in withClass(root, 'mechanism-lone-pair'):
            for center, radius, limit in (('cx', 'rx', width), ('cy', 'ry', height)):
              c, r = float(dot.get(center)), float(dot.get(radius))
              self.assertGreaterEqual(c - r, 0.5)
              self.assertLessEqual(c + r, limit - 0.5)

  def testConcertedArrowChemistryAndAntiGeometry(self):
    scenes = {scene.slug: scene for scene in self.scenes}
    anchor_type = rdMolDraw2D.MechanismAnchorType
    sn2 = scenes['sn2']
    self.assertEqual(Chem.MolToSmiles(sn2.mols[1]), 'CBr')
    self.assertEqual((sn2.arrows[0].start.fragment, sn2.arrows[0].start.index), (0, 0))
    self.assertEqual((sn2.arrows[0].end.fragment, sn2.arrows[0].end.index), (1, 0))
    self.assertEqual(sn2.arrows[0].start.type, anchor_type.LonePair)
    self.assertEqual(sn2.arrows[1].start.type, anchor_type.Bond)
    self.assertEqual((sn2.arrows[1].start.index, sn2.arrows[1].end.index), (0, 1))

    e2 = scenes['e2']
    substrate = e2.mols[1]
    self.assertEqual([substrate.GetAtomWithIdx(i).GetAtomicNum()
                      for i in (0, 1, 3, 5)], [1, 6, 6, 35])
    conf = substrate.GetConformer()
    h, beta, alpha, bromine = [
      np.array(tuple(conf.GetAtomPosition(i))) for i in (0, 1, 3, 5)]
    axis = alpha - beta
    axis /= np.linalg.norm(axis)
    h_vector, br_vector = h - beta, bromine - alpha
    h_projection = h_vector - np.dot(h_vector, axis) * axis
    br_projection = br_vector - np.dot(br_vector, axis) * axis
    cosine = np.dot(h_projection, br_projection) / (
      np.linalg.norm(h_projection) * np.linalg.norm(br_projection))
    self.assertAlmostEqual(cosine, -1, places=7)
    self.assertEqual((e2.arrows[0].end.fragment, e2.arrows[0].end.index), (1, 0))
    for arrow, source, target in ((e2.arrows[1], 0, 2), (e2.arrows[2], 4, 5)):
      self.assertEqual(arrow.start.type, anchor_type.Bond)
      self.assertEqual((arrow.start.fragment, arrow.start.index), (1, source))
      self.assertEqual((arrow.end.fragment, arrow.end.index), (1, target))
    self.assertEqual(e2.arrows[1].end.type, anchor_type.Bond)
    self.assertEqual(e2.arrows[2].end.type, anchor_type.Atom)

    transfer = scenes['proton-transfer']
    nitrogen = transfer.mols[0].GetAtomWithIdx(0)
    self.assertEqual((nitrogen.GetAtomicNum(), nitrogen.GetFormalCharge()), (7, 0))
    self.assertEqual((transfer.arrows[0].end.fragment, transfer.arrows[0].end.index), (1, 4))
    self.assertEqual(transfer.mols[1].GetAtomWithIdx(4).GetAtomicNum(), 1)
    self.assertEqual((transfer.arrows[1].start.index, transfer.arrows[1].end.index), (3, 3))

  def testCarbonylStepsHaveSuppliedBalancedStructures(self):
    scenes = {scene.slug: scene for scene in self.scenes}
    addition, protonation, product = [scenes[key] for key in (
      'carbonyl-addition', 'carbonyl-protonation', 'carbonyl-product')]

    def formula(mols):
      elements = Counter()
      for mol in mols:
        for atom in mol.GetAtoms():
          elements[atom.GetAtomicNum()] += 1
          elements[1] += atom.GetNumImplicitHs() + atom.GetNumExplicitHs()
      return +elements

    self.assertEqual(formula(addition.mols), formula(protonation.mols[:1]))
    self.assertEqual(formula(protonation.mols), formula(product.mols))
    for scene in (addition, protonation, product):
      self.assertEqual(sum(Chem.GetFormalCharge(mol) for mol in scene.mols), -1)
    self.assertEqual(addition.mols[0].GetBondWithIdx(1).GetBondType(), Chem.BondType.DOUBLE)
    self.assertEqual(protonation.mols[0].GetAtomWithIdx(2).GetFormalCharge(), -1)
    self.assertEqual(Chem.GetFormalCharge(product.mols[0]), 0)
    self.assertEqual(product.arrows, [])
    self.assertEqual((addition.arrows[0].start.fragment, addition.arrows[0].start.index), (1, 0))
    self.assertEqual((addition.arrows[0].end.fragment, addition.arrows[0].end.index), (0, 1))
    self.assertEqual((addition.arrows[1].start.index, addition.arrows[1].end.index), (1, 2))
    self.assertEqual((protonation.arrows[0].start.fragment,
                      protonation.arrows[0].start.index), (0, 2))
    self.assertEqual((protonation.arrows[0].end.fragment,
                      protonation.arrows[0].end.index), (1, 0))

  def testGalleryArtifactsAndManifest(self):
    with tempfile.TemporaryDirectory() as directory:
      output = Path(directory)
      manifest = generate_gallery(self.scenes, output)
      self.assertEqual(manifest, json.loads((output / 'manifest.json').read_text()))
      self.assertEqual(manifest['schemaVersion'], 1)
      self.assertEqual(len(manifest['scenes']), 8)
      document = (output / 'index.html').read_text(encoding='utf-8')
      self.assertEqual(document.count('<article '), 8)
      self.assertEqual(document.count('<img '), 8)
      self.assertNotIn('<svg', document)
      for entry, scene in zip(manifest['scenes'], self.scenes):
        self.assertEqual(entry['files'], {'svg': scene.slug + '.svg'})
        self.assertEqual(entry['originalAtoms'], [mol.GetNumAtoms() for mol in scene.mols])
        self.assertEqual(entry['arrows'], len(scene.arrows))
        self.assertEqual(entry['fragments'], len(scene.mols))
        self.assertIn(entry['files']['svg'], document)
        root = ET.parse(output / entry['files']['svg']).getroot()
        self.assertEqual(root.tag, SVG + 'svg')
        self.assertEqual(float(root.get('width').removesuffix('px')), entry['width'])
        self.assertEqual(float(root.get('height').removesuffix('px')), entry['height'])

  @unittest.skipUnless(hasattr(rdMolDraw2D, 'MolDraw2DCairo'), 'Cairo support not enabled')
  def testNativeCairoGallery(self):
    with tempfile.TemporaryDirectory() as directory:
      manifest = generate_gallery(self.scenes, directory, png=True)
      for entry in manifest['scenes']:
        with self.subTest(scene=entry['slug']):
          png = (Path(directory) / entry['files']['png']).read_bytes()
          self.assertEqual(png[:8], b'\x89PNG\r\n\x1a\n')
          self.assertEqual(png[12:16], b'IHDR')
          self.assertEqual(struct.unpack('>II', png[16:24]), (entry['width'], entry['height']))
          self.assertIn(b'IDAT', png)
          self.assertEqual(png[-8:-4], b'IEND')

  def testOutputValidationAndHtmlEscaping(self):
    with tempfile.TemporaryDirectory() as directory:
      scene = self.scenes[0]
      for scenes in ([], [scene, scene], [replace(scene, slug='../outside')]):
        with self.subTest(scenes=scenes), self.assertRaises(ValueError):
          generate_gallery(scenes, directory)
      self.assertEqual(list(Path(directory).iterdir()), [])
      special = replace(scene, title='Attack <A> & B', description='"Donor" < acceptor')
      generate_gallery([special], directory)
      document = (Path(directory) / 'index.html').read_text(encoding='utf-8')
      self.assertIn('Attack &lt;A&gt; &amp; B', document)
      self.assertIn('&quot;Donor&quot; &lt; acceptor', document)
      self.assertNotIn('Attack <A>', document)

  def testCliSelectionAndErrors(self):
    with tempfile.TemporaryDirectory() as directory:
      output = Path(directory) / 'gallery'
      with contextlib.redirect_stdout(io.StringIO()) as listing:
        self.assertEqual(main(['--list', '--output', str(output)]), 0)
      self.assertEqual(len(listing.getvalue().splitlines()), 8)
      self.assertFalse(output.exists())
      with contextlib.redirect_stdout(io.StringIO()):
        self.assertEqual(main(['--only', 'e2', '--only', 'homolysis',
                               '--output', str(output)]), 0)
      manifest = json.loads((output / 'manifest.json').read_text())
      self.assertEqual([scene['slug'] for scene in manifest['scenes']], ['e2', 'homolysis'])
      with contextlib.redirect_stderr(io.StringIO()) as errors:
        with self.assertRaises(SystemExit) as raised:
          main(['--only', 'no-such-scene', '--output', str(output)])
      self.assertEqual(raised.exception.code, 2)
      self.assertIn('Unknown scenes', errors.getvalue())

  def testExampleHelperRejectsIncorrectCoordinateCounts(self):
    with self.assertRaises(ValueError):
      positioned_molecule('CC', [(0, 0)])
    water = positioned_molecule('[H]O[H]', [(-1, 0), (0, 0), (1, 0)])
    self.assertEqual([atom.GetAtomicNum() for atom in water.GetAtoms()], [1, 8, 1])
    with self.assertRaises(ValueError):
      render_scene(self.scenes[0], 'pdf')


if __name__ == '__main__':
  unittest.main()
