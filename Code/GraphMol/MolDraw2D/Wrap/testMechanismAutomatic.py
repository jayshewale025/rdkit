#
# Copyright (C) 2026 RDKit contributors
# This file is part of the RDKit, covered by the BSD license in license.txt.
#
import itertools
import math
from pathlib import Path
import re
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

from rdkit import Chem, Geometry
from rdkit.Chem.Draw import rdMolDraw2D

from testMechanismDrawing import (Anchor, Arrow, Head, P, affineMap, molecule,
                                  pathPoints, snapshot, withClass)

DATA = Path(__file__).resolve().parents[4] / 'Docs' / 'Book' / 'data'
sys.path.insert(0, str(DATA))
from mechanism_automatic import scenes
from mechanism_gallery import generate_gallery


def smilesMol(smiles):
  params = Chem.SmilesParserParams()
  params.removeHs = False
  mol = Chem.MolFromSmiles(smiles, params)
  assert mol is not None
  assert mol.GetNumConformers() == 0
  return mol


def anchorState(anchor):
  return (anchor.type, anchor.fragment, anchor.index,
          (anchor.offset.x, anchor.offset.y), anchor.fraction,
          anchor.clipToLabel, anchor.autoPosition, anchor.lonePairAngle,
          anchor.lonePairSeparation, anchor.lonePairRadius)


def arrowState(arrow):
  return (anchorState(arrow.start), anchorState(arrow.end),
          (arrow.control1.x, arrow.control1.y),
          (arrow.control2.x, arrow.control2.y), arrow.head,
          arrow.headLength, arrow.headWidth, arrow.lineWidth, arrow.colour,
          arrow.autoRoute, arrow.autoColour)


def moleculeState(mol):
  # Include stereochemistry, computed/private properties and every conformer,
  # not just the graph and coordinates exposed by the older snapshot helper.
  flags = Chem.PropertyPickleOptions.AllProps | Chem.PropertyPickleOptions.CoordsAsDouble
  return snapshot(mol), mol.ToBinary(flags)


def distance(first, second):
  return math.hypot(first.x - second.x, first.y - second.y)


def stroke(node):
  return re.search(r'(?:^|;)stroke:([^;]+)', node.get('style', '')).group(1)


class TestMechanismAutomatic(unittest.TestCase):

  def render(self, mols, arrows, size=(600, 360), font=20, **kwargs):
    drawer = rdMolDraw2D.MolDraw2DSVG(*size, -1, -1, True)
    drawer.drawOptions().fixedFontSize = font
    positions = drawer.DrawMechanism(mols, arrows, **kwargs)
    drawer.FinishDrawing()
    return drawer, positions, ET.fromstring(drawer.GetDrawingText())

  def assertPoint(self, actual, expected, tolerance=0.11):
    self.assertAlmostEqual(actual.x, expected.x, delta=tolerance)
    self.assertAlmostEqual(actual.y, expected.y, delta=tolerance)

  def assertGeometry(self, drawer, positions, root, mols, arrows):
    self.assertEqual([len(group) for group in positions],
                     [mol.GetNumAtoms() for mol in mols])
    width, height = drawer.Width(), drawer.Height()
    self.assertGreater(width, 0)
    self.assertGreater(height, 0)

    def contained(point):
      self.assertTrue(math.isfinite(point.x) and math.isfinite(point.y))
      self.assertGreaterEqual(point.x, -0.11)
      self.assertLessEqual(point.x, width + 0.11)
      self.assertGreaterEqual(point.y, -0.11)
      self.assertLessEqual(point.y, height + 0.11)

    for point in itertools.chain.from_iterable(positions):
      contained(point)
    shafts, heads = withClass(root, 'mechanism-shaft'), withClass(root, 'mechanism-head')
    self.assertEqual(len(shafts), len(arrows))
    self.assertEqual(len(heads), len(arrows))
    for index, (shaft, head) in enumerate(zip(shafts, heads)):
      self.assertIn(f'mechanism-arrow-{index}', shaft.get('class', '').split())
      self.assertIn(f'mechanism-arrow-{index}', head.get('class', '').split())
      self.assertEqual(shaft.get('d').count('C'), 1)
      self.assertNotIn('L', shaft.get('d'))
      self.assertIn('fill:none', shaft.get('style'))
      curve = pathPoints(shaft)
      self.assertEqual(len(curve), 4)
      self.assertTrue(all(math.isfinite(p.x) and math.isfinite(p.y) for p in curve))
      self.assertGreater(distance(curve[0], curve[1]), 1e-6)
      self.assertGreater(distance(curve[2], curve[3]), 1e-6)
      # Control points need not fit inside the canvas; the actual curve must.
      for step in range(101):
        t, u = step / 100, 1 - step / 100
        contained(curve[0] * (u**3) + curve[1] * (3*u*u*t) +
                  curve[2] * (3*u*t*t) + curve[3] * (t**3))
      vertices = pathPoints(head)
      self.assertTrue(vertices)
      self.assertLess(min(distance(p, curve[3]) for p in vertices), 0.2)
      self.assertEqual(stroke(head), stroke(shaft))
      for point in vertices:
        contained(point)
    for dot in withClass(root, 'mechanism-lone-pair'):
      center = P(float(dot.get('cx')), float(dot.get('cy')))
      rx, ry = float(dot.get('rx', dot.get('r'))), float(dot.get('ry', dot.get('r')))
      self.assertGreater(rx, 0)
      self.assertGreater(ry, 0)
      contained(P(center.x - rx, center.y - ry))
      contained(P(center.x + rx, center.y + ry))

  def assertSeparated(self, positions, font):
    # A font-sized envelope adds label clearance, unlike checking atom boxes
    # alone (which are points for hydroxide and a proton).
    margin = 0.75 * font
    boxes = [(min(p.x for p in group) - margin, min(p.y for p in group) - margin,
              max(p.x for p in group) + margin, max(p.y for p in group) + margin)
             for group in positions]
    for first, second in itertools.combinations(boxes, 2):
      self.assertTrue(first[2] < second[0] or second[2] < first[0] or
                      first[3] < second[1] or second[3] < first[1],
                      f'Fragment envelopes overlap: {first}, {second}')

  def assertGlobalPlacement(self, positions, scene):
    transform, scale = affineMap(positions[0], scene[0])
    self.assertGreater(scale, 0)
    for actual, expected in zip(itertools.chain.from_iterable(positions),
                                itertools.chain.from_iterable(scene)):
      self.assertPoint(actual, transform(expected), 1e-7)
    return transform

  def testAutomaticAndManualFactories(self):
    start, end = Anchor.LonePair(fragment=0, atomIndex=1), Anchor.Atom(1, 0)
    self.assertTrue(start.autoPosition)
    self.assertTrue(Anchor.LonePair(0, 1, offset=None).autoPosition)
    self.assertFalse(Anchor.LonePair(0, 1, P(0, 0)).autoPosition)
    self.assertFalse(Anchor.LonePair(0, 1, P(0.4, 0.2), angle=30).autoPosition)
    self.assertFalse(Anchor.Atom(0, 1).autoPosition)
    for head in (Head.Full, Head.FishhookLeft, Head.FishhookRight):
      arrow = Arrow(start=start, end=end, head=head)
      self.assertTrue(arrow.autoRoute)
      self.assertTrue(arrow.autoColour)
      self.assertEqual(arrow.head, head)
    self.assertEqual(Arrow(start, end).head, Head.Full)
    manual = Arrow(start, end, P(0, 2), P(2, 2))
    self.assertFalse(manual.autoRoute)
    self.assertFalse(manual.autoColour)
    arrow.colour = (0.25, 0.5, 0.75)
    self.assertFalse(arrow.autoColour)
    self.assertTrue(arrow.autoRoute)
    arrow.autoColour = True
    self.assertTrue(arrow.autoColour)
    arrow.autoRoute = False
    self.assertFalse(arrow.autoRoute)
    start.autoPosition = False
    self.assertFalse(start.autoPosition)

  def testHydroxideToProtonWithoutCoordinates(self):
    mols = [smilesMol('[OH-]'), smilesMol('[H+]')]
    self.assertEqual([mol.GetNumAtoms() for mol in mols], [1, 1])
    before = [moleculeState(mol) for mol in mols]
    arrow = Arrow(Anchor.LonePair(0, 0), Anchor.Atom(1, 0))
    for size, font in (((600, 360), 20), ((600, 360), 28), ((-1, -1), 20)):
      with self.subTest(size=size, font=font):
        drawer, positions, root = self.render(mols, [arrow], size=size, font=font)
        self.assertGeometry(drawer, positions, root, mols, [arrow])
        self.assertSeparated(positions, font)
        dots = withClass(root, 'mechanism-lone-pair')
        self.assertEqual(len(dots), 2)
        centers = [P(float(dot.get('cx')), float(dot.get('cy'))) for dot in dots]
        pair = (centers[0] + centers[1]) / 2
        curve = pathPoints(withClass(root, 'mechanism-shaft')[0])
        self.assertPoint(curve[0], pair)
        chord = curve[3] - curve[0]
        midpoint = (curve[0] + curve[1] * 3 + curve[2] * 3 + curve[3]) / 8
        bow = abs(chord.x * (midpoint.y - curve[0].y) -
                  chord.y * (midpoint.x - curve[0].x)) / distance(curve[0], curve[3])
        self.assertGreater(bow, 1)
        # The simple two-endpoint route must not double back along its chord.
        samples = [curve[0] * ((1-t)**3) + curve[1] * (3*(1-t)**2*t) +
                   curve[2] * (3*(1-t)*t*t) + curve[3] * (t**3)
                   for t in (step / 100 for step in range(101))]
        projections = [((p.x - curve[0].x) * chord.x +
                        (p.y - curve[0].y) * chord.y) /
                       (chord.x**2 + chord.y**2) for p in samples]
        self.assertTrue(all(-0.01 <= t <= 1.01 for t in projections))
        self.assertTrue(all(a <= b + 0.001 for a, b in zip(projections, projections[1:])))
        self.assertGreater(distance(pair, positions[0][0]), 0.25 * font)
        self.assertLess(distance(pair, positions[0][0]),
                        0.5 * distance(positions[0][0], positions[1][0]))
        self.assertGreater(distance(curve[-1], positions[1][0]), 1)
        self.assertGreater(distance(centers[0], centers[1]),
                           2 * float(dots[0].get('rx')))
    self.assertEqual([moleculeState(mol) for mol in mols], before)
    self.assertTrue(all(mol.GetNumConformers() == 0 for mol in mols))

  def testCarbonylAdditionPaletteDeterminismAndOverrides(self):
    mols = [smilesMol('CC=O'), smilesMol('[OH-]')]
    arrows = [Arrow(Anchor.LonePair(1, 0), Anchor.Atom(0, 1)),
              Arrow(Anchor.Bond(0, 1), Anchor.Atom(0, 2))]
    before = [arrowState(arrow) for arrow in arrows]
    drawer, positions, root = self.render(mols, arrows)
    self.assertGeometry(drawer, positions, root, mols, arrows)
    self.assertSeparated(positions, 20)
    self.assertEqual(len({stroke(p) for p in withClass(root, 'mechanism-shaft')}), 2)
    repeated, repeatedPositions, _ = self.render(mols, arrows)
    self.assertEqual(drawer.GetDrawingText(), repeated.GetDrawingText())
    for actual, expected in zip(itertools.chain.from_iterable(repeatedPositions),
                                itertools.chain.from_iterable(positions)):
      self.assertPoint(actual, expected, 1e-9)
    self.assertEqual([arrowState(arrow) for arrow in arrows], before)
    arrows[0].control1, arrows[0].control2 = P(100, -200), P(-300, 400)
    sentinel = [arrowState(arrow) for arrow in arrows]
    ignored, _, _ = self.render(mols, arrows)
    self.assertEqual(drawer.GetDrawingText(), ignored.GetDrawingText())
    self.assertEqual([arrowState(arrow) for arrow in arrows], sentinel)
    arrows[0].colour = (0.25, 0.5, 0.75, 0.5)
    self.assertFalse(arrows[0].autoColour)
    _, _, overridden = self.render(mols, arrows)
    self.assertEqual(stroke(withClass(overridden, 'mechanism-shaft')[0]), '#3F7FBF7F')
    self.assertEqual(stroke(withClass(overridden, 'mechanism-shaft')[1]),
                     stroke(withClass(root, 'mechanism-shaft')[1]))

  def testAutomaticMultipleBondTailsUseOutermostStroke(self):
    for smiles, bondTarget in (
        ('CC=O', False), ('C=O', False), ('C#N', False), ('c1ccccc1', True)):
      for rotation, spacing, split in ((0, 0.15, False), (37, 0.24, True)):
        with self.subTest(smiles=smiles, rotation=rotation, spacing=spacing,
                          split=split):
          mol = smilesMol(smiles)
          Chem.Kekulize(mol, clearAromaticFlags=True)
          source = next(bond for bond in mol.GetBonds() if bond.GetBondType() in
                        (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE))
          bondIndex = source.GetIdx()
          target = Anchor.Atom(0, source.GetEndAtomIdx())
          if bondTarget:
            neighbor = next(bond for bond in source.GetEndAtom().GetBonds()
                            if bond.GetBondType() == Chem.BondType.SINGLE)
            target = Anchor.Bond(0, neighbor.GetIdx())
          arrow = Arrow(Anchor.Bond(0, bondIndex), target)
          drawer = rdMolDraw2D.MolDraw2DSVG(500, 320, -1, -1, True)
          options = drawer.drawOptions()
          options.rotate = rotation
          options.multipleBondOffset = spacing
          options.splitBonds = split
          positions = drawer.DrawMechanism([mol], [arrow])
          drawer.FinishDrawing()
          root = ET.fromstring(drawer.GetDrawingText())
          self.assertGeometry(drawer, positions, root, [mol], [arrow])
          curve = pathPoints(withClass(root, 'mechanism-shaft')[0])
          tail, direction = curve[0], curve[1] - curve[0]
          direction /= math.hypot(direction.x, direction.y)
          strokes = [pathPoints(node) for node in withClass(root, f'bond-{bondIndex}')]
          self.assertGreaterEqual(len(strokes), 2)
          distances = []
          for start, end in strokes:
            segment = end - start
            offset = tail - start
            lengthSq = segment.x**2 + segment.y**2
            self.assertGreater(lengthSq, 0)
            fraction = max(0, min(1, (offset.x * segment.x +
                                      offset.y * segment.y) / lengthSq))
            distances.append(distance(tail, start + segment * fraction))
            denominator = direction.x * segment.y - direction.y * segment.x
            if abs(denominator) > 1e-8:
              relative = start - tail
              alongRay = (relative.x * segment.y -
                          relative.y * segment.x) / denominator
              alongStroke = (relative.x * direction.y -
                             relative.y * direction.x) / denominator
              if 0 <= alongStroke <= 1:
                self.assertLessEqual(alongRay, 0.2)
          self.assertLessEqual(min(distances), 0.15)

  def testCoincidentPlacementCandidateDoesNotAbortAutomaticLayout(self):
    # Geometry stress case, not a chemical mechanism: one explored placement
    # coincides with the explicit source offset, but other placements are valid.
    mols = [Chem.MolFromSmiles('[Cl-]') for _ in range(2)]
    before = [moleculeState(mol) for mol in mols]
    arrows = [Arrow(Anchor.Atom(0, 0), Anchor.Atom(1, 0)),
              Arrow(Anchor.Atom(0, 0, P(1.6, 0)), Anchor.Atom(1, 0))]
    drawer = rdMolDraw2D.MolDraw2DSVG(500, 300)
    positions = drawer.DrawMechanism(mols, arrows)
    drawer.FinishDrawing()
    root = ET.fromstring(drawer.GetDrawingText())
    self.assertGeometry(drawer, positions, root, mols, arrows)
    self.assertEqual(len(withClass(root, 'mechanism-shaft')), 2)
    self.assertEqual([moleculeState(mol) for mol in mols], before)
    self.assertTrue(all(mol.GetNumConformers() == 0 for mol in mols))

  def testIntramolecularFlowsAndMultipleFragmentHub(self):
    scenarios = [
      ([smilesMol('C=CC=O')],
       [Arrow(Anchor.Bond(0, 0), Anchor.Bond(0, 1)),
        Arrow(Anchor.Bond(0, 2), Anchor.Atom(0, 3))]),
      ([smilesMol('O=C=O'), smilesMol('[OH-]'), smilesMol('N'), smilesMol('[Cl-]')],
       [Arrow(Anchor.LonePair(1, 0), Anchor.Atom(0, 1)),
        Arrow(Anchor.LonePair(2, 0), Anchor.Atom(0, 1), Head.FishhookLeft),
        Arrow(Anchor.LonePair(3, 0), Anchor.Atom(0, 1), Head.FishhookRight)]),
    ]
    for mols, arrows in scenarios:
      with self.subTest(fragments=len(mols)):
        before = [moleculeState(mol) for mol in mols]
        drawer, positions, root = self.render(mols, arrows, size=(800, 480))
        self.assertGeometry(drawer, positions, root, mols, arrows)
        self.assertSeparated(positions, 20)
        self.assertEqual(len({stroke(p) for p in withClass(root, 'mechanism-shaft')}),
                         len(arrows))
        self.assertEqual([moleculeState(mol) for mol in mols], before)

  def testExistingConformersRemainRigidAndInputsRemainExact(self):
    first = molecule('[13CH3][C@H](O)/C=C/F',
                     [(0, 0), (1.5, 0), (1.5, 1.5), (3, -0.5), (4, 0), (5, 1)])
    first.SetProp('mechanism-test', 'retained')
    first.GetAtomWithIdx(2).SetAtomMapNum(17)
    first.GetAtomWithIdx(1).SetProp('_source', 'chiral center')
    first.GetBondWithIdx(3).SetProp('source', 'stereo bond')
    selected = Chem.Conformer(first.GetConformer())
    selected.SetId(7)
    selected.SetProp('source', 'selected conformer')
    selected.SetAtomPosition(5, Geometry.Point3D(5, 2, 0))
    first.AddConformer(selected, assignId=False)
    second = molecule('CO', [(0, 0), (0, 2)])
    mols = [first, second]
    before = [moleculeState(mol) for mol in mols]
    arrows = [Arrow(Anchor.LonePair(0, 2), Anchor.Atom(1, 0))]
    drawer, positions, root = self.render(mols, arrows, confIds=[7, -1])
    self.assertGeometry(drawer, positions, root, mols, arrows)
    self.assertEqual([moleculeState(mol) for mol in mols], before)
    scale = distance(positions[0][0], positions[0][1]) / 1.5
    for mol, confId, group in zip(mols, (7, -1), positions):
      conf = mol.GetConformer(confId)
      for i, j in itertools.combinations(range(mol.GetNumAtoms()), 2):
        expected = distance(conf.GetAtomPosition(i), conf.GetAtomPosition(j)) * scale
        self.assertAlmostEqual(distance(group[i], group[j]), expected, delta=1e-7)
    a, b, c = positions[0][:3]
    # A canvas Y flip reverses the signed area; a rigid rotation cannot mirror it again.
    self.assertLess((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x), 0)

  def testLayoutSuppressionPreservesOneGlobalPlacement(self):
    mols = [molecule('CC', [(0, 0), (1.5, 0)]),
            molecule('CC', [(5, 1), (6.5, 1)])]
    scene = [[P(0, 0), P(1.5, 0)], [P(5, 1), P(6.5, 1)]]
    automatic = Arrow(Anchor.Atom(0, 1), Anchor.Atom(1, 0))
    manual = Arrow(Anchor.Atom(0, 0), Anchor.Atom(1, 1), P(0, 3), P(7, 3))
    cases = [
      ('disabled', [automatic], {'autoLayout': False}),
      ('offsets', [automatic], {'offsets': [P(0, 0), P(0, 0)]}),
      ('rotations', [automatic], {'rotations': [0, 0]}),
      ('mixed', [automatic, manual], {}),
      ('absolute source', [automatic, Arrow(Anchor.Point(P(-2, 2)), Anchor.Atom(0, 0))], {}),
      ('absolute target', [automatic, Arrow(Anchor.Atom(0, 0), Anchor.Point(P(-2, 2)))], {}),
      ('no arrows', [], {}),
    ]
    for name, arrows, kwargs in cases:
      with self.subTest(case=name):
        drawer, positions, root = self.render(mols, arrows, **kwargs)
        self.assertGeometry(drawer, positions, root, mols, arrows)
        transform = self.assertGlobalPlacement(positions, scene)
        if name == 'mixed':
          curve = pathPoints(withClass(root, 'mechanism-shaft')[1])
          self.assertPoint(curve[1], transform(manual.control1))
          self.assertPoint(curve[2], transform(manual.control2))
    drawer, positions, root = self.render(
      mols, [automatic], offsets=[P(1, -2), P(-2, 3)], rotations=[90, 0])
    self.assertGlobalPlacement(positions, [[P(1, -2), P(1, -0.5)], [P(3, 4), P(4.5, 4)]])
    self.assertGeometry(drawer, positions, root, mols, [automatic])

  def testAutomaticLonePairsAndSameAtomLoop(self):
    mol = smilesMol('NCC')
    pair = Anchor.LonePair(0, 0)
    arrows = [Arrow(pair, Anchor.Atom(0, 1)), Arrow(pair, Anchor.Atom(0, 2))]
    before = [arrowState(arrow) for arrow in arrows]
    drawer, positions, root = self.render([mol], arrows)
    self.assertGeometry(drawer, positions, root, [mol], arrows)
    self.assertIn(len(withClass(root, 'mechanism-lone-pair')), (2, 4))
    self.assertEqual([arrowState(arrow) for arrow in arrows], before)
    loop = Arrow(Anchor.Atom(0, 1), Anchor.Atom(0, 1))
    drawer, positions, root = self.render([mol], [loop])
    self.assertGeometry(drawer, positions, root, [mol], [loop])
    curve = pathPoints(withClass(root, 'mechanism-shaft')[0])
    atom = positions[0][1]
    for endpoint in (curve[0], curve[3]):
      self.assertLess(distance(endpoint, atom),
                      0.5 * distance(positions[0][0], atom))
    midpoint = (curve[0] + curve[1] * 3 + curve[2] * 3 + curve[3]) / 8
    self.assertGreater(distance(midpoint, atom),
                       max(distance(curve[0], atom), distance(curve[3], atom)) + 1)

  def testCoincidentFixedPointsProduceLoopOrClearError(self):
    mol = smilesMol('CC')
    arrow = Arrow(Anchor.Point(P(3, 2)), Anchor.Point(P(3, 2)))
    try:
      drawer, positions, root = self.render([mol], [arrow], autoLayout=False)
    except ValueError as error:
      self.assertRegex(str(error).lower(), r'coincident|degenerate|zero|tangent|rout|distinct')
    else:
      self.assertGeometry(drawer, positions, root, [mol], [arrow])
      curve = pathPoints(withClass(root, 'mechanism-shaft')[0])
      self.assertPoint(curve[0], curve[3])
      midpoint = (curve[0] + curve[1] * 3 + curve[2] * 3 + curve[3]) / 8
      self.assertGreater(distance(midpoint, curve[0]), 1)

  def testAutomaticValidationBeforeDrawing(self):
    mol = smilesMol('CC=O')
    for factory in (Anchor.Atom, Anchor.Bond, Anchor.LonePair):
      for fragment, index in ((-1, 0), (0, -1)):
        with self.subTest(factory=factory, fragment=fragment, index=index):
          with self.assertRaises(ValueError):
            factory(fragment, index)
    invalid = [Anchor.Atom(1, 0), Anchor.Atom(0, 3), Anchor.Bond(0, 2),
               Anchor.LonePair(1, 0), Anchor.LonePair(0, 3), Anchor.Bond(0, 0, 1.1)]
    for field, value in (('lonePairRadius', 0), ('lonePairRadius', float('nan')),
                         ('lonePairSeparation', -0.1), ('lonePairSeparation', float('inf')),
                         ('lonePairSeparation', 0.04)):
      pair = Anchor.LonePair(0, 2)
      setattr(pair, field, value)
      invalid.append(pair)
    for anchor in invalid:
      for atStart in (False, True):
        with self.subTest(anchor=anchorState(anchor), atStart=atStart):
          drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
          before = drawer.GetDrawingText()
          other = Anchor.Atom(0, 0)
          arrow = Arrow(anchor, other) if atStart else Arrow(other, anchor)
          with self.assertRaises(ValueError):
            drawer.DrawMechanism([mol], [arrow])
          self.assertEqual(drawer.GetDrawingText(), before)

  def testAutomaticGalleryUsesEndpointOnlyScenesAndAccurateIntroduction(self):
    gallery = scenes()
    self.assertEqual(len(gallery), 7)
    mols = [mol for scene in gallery for mol in scene.mols]
    before = [moleculeState(mol) for mol in mols]
    for scene in gallery:
      with self.subTest(scene=scene.slug):
        self.assertFalse(scene.offsets)
        self.assertFalse(scene.rotations)
        self.assertTrue(scene.arrows)
        self.assertTrue(all(arrow.autoRoute for arrow in scene.arrows))
        self.assertTrue(all(mol.GetNumConformers() == 0 for mol in scene.mols))
    with tempfile.TemporaryDirectory(prefix='mechanism-automatic-gallery-', dir='.') as directory:
      output = Path(directory)
      manifest = generate_gallery(gallery, output, png=False)
      self.assertEqual(len(manifest['scenes']), 7)
      entries = {entry['slug']: entry for entry in manifest['scenes']}
      for scene in gallery:
        with self.subTest(scene=scene.slug):
          entry = entries[scene.slug]
          self.assertEqual(set(entry['files']), {'svg'})
          self.assertEqual(entry['arrows'], len(scene.arrows))
          root = ET.parse(output / entry['files']['svg']).getroot()
          self.assertEqual(len(withClass(root, 'mechanism-shaft')), len(scene.arrows))
          self.assertEqual(len(withClass(root, 'mechanism-head')), len(scene.arrows))
      page = (output / 'index.html').read_text(encoding='utf-8')
      self.assertIn('automatic', page.lower())
      self.assertIn('heuristic', page.lower())
      self.assertNotIn(
        'Every structure, electron-flow endpoint, and curve is explicitly authored', page)
      self.assertFalse(list(output.glob('*.png')))
    self.assertTrue(all(mol.GetNumConformers() == 0 for mol in mols))
    self.assertEqual([moleculeState(mol) for mol in mols], before)

  @unittest.skipUnless(hasattr(rdMolDraw2D, 'MolDraw2DCairo'), 'Cairo not enabled')
  def testAutomaticCairoAndSvg(self):
    mols = [smilesMol('CCC'), smilesMol('CC')]
    arrows = [Arrow(Anchor.Bond(0, 0), Anchor.Atom(1, 0))]
    positions = []
    for backend in (rdMolDraw2D.MolDraw2DSVG, rdMolDraw2D.MolDraw2DCairo):
      drawer = backend(600, 360, -1, -1, True)
      positions.append(drawer.DrawMechanism(mols, arrows))
      drawer.FinishDrawing()
      if backend is rdMolDraw2D.MolDraw2DCairo:
        self.assertTrue(drawer.GetDrawingText().startswith(b'\x89PNG\r\n\x1a\n'))
        self.assertGreater(len(drawer.GetDrawingText()), 100)
    for actual, expected in zip(itertools.chain.from_iterable(positions[0]),
                                itertools.chain.from_iterable(positions[1])):
      self.assertPoint(actual, expected, 0.5)
    charged = [smilesMol('[OH-]'), smilesMol('[H+]')]
    arrow = Arrow(Anchor.LonePair(0, 0), Anchor.Atom(1, 0))
    drawer = rdMolDraw2D.MolDraw2DCairo(600, 360, -1, -1, True)
    chargedPositions = drawer.DrawMechanism(charged, [arrow])
    self.assertSeparated(chargedPositions, 20)
    for point in itertools.chain.from_iterable(chargedPositions):
      self.assertTrue(math.isfinite(point.x) and math.isfinite(point.y))
      self.assertTrue(0 <= point.x <= 600 and 0 <= point.y <= 360)
    drawer.FinishDrawing()
    self.assertTrue(drawer.GetDrawingText().startswith(b'\x89PNG\r\n\x1a\n'))


if __name__ == '__main__':
  unittest.main()
