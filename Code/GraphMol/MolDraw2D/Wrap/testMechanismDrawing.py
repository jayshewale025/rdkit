#
# Copyright (C) 2026 RDKit contributors
# This file is part of the RDKit, covered by the BSD license in license.txt.
#
import math
import re
import unittest
import xml.etree.ElementTree as ET

from rdkit import Chem, Geometry
from rdkit.Chem.Draw import rdMolDraw2D


P = Geometry.Point2D
Anchor = rdMolDraw2D.MechanismAnchor
Arrow = rdMolDraw2D.MechanismArrow
Head = rdMolDraw2D.ElectronArrowHead
SVG = '{http://www.w3.org/2000/svg}'
NUMBER = r'[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?'


def molecule(smiles, coordinates):
  mol = Chem.MolFromSmiles(smiles)
  conf = Chem.Conformer(mol.GetNumAtoms())
  conf.Set3D(False)
  for idx, (x, y) in enumerate(coordinates):
    conf.SetAtomPosition(idx, Geometry.Point3D(x, y, 0))
  mol.AddConformer(conf, assignId=True)
  return mol


def snapshot(mol):
  return (
    mol.GetNumAtoms(), mol.GetNumBonds(),
    mol.GetPropsAsDict(includePrivate=True, includeComputed=False),
    [(atom.GetAtomicNum(), atom.GetFormalCharge(), atom.GetAtomMapNum(),
      atom.GetPropsAsDict(includePrivate=True, includeComputed=False))
     for atom in mol.GetAtoms()],
    [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType(),
      bond.GetBondDir(), bond.GetStereo(),
      bond.GetPropsAsDict(includePrivate=True, includeComputed=False))
     for bond in mol.GetBonds()],
    [(conf.GetId(), conf.Is3D(),
      [(pos.x, pos.y, pos.z) for pos in
       (conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms()))])
     for conf in mol.GetConformers()],
  )


def pathPoints(path):
  values = [float(value) for value in re.findall(NUMBER, path.attrib['d'])]
  return [P(x, y) for x, y in zip(values[::2], values[1::2])]


def cubics(root):
  return [path for path in root.iter(SVG + 'path')
          if 'C' in path.attrib.get('d', '') and
          'fill:none' in path.attrib.get('style', '')]


def withClass(root, name):
  return [node for node in root.iter() if name in node.get('class', '').split()]


def affineMap(positions, scene):
  dx, dy = scene[1].x - scene[0].x, scene[1].y - scene[0].y
  scale = ((positions[1].x - positions[0].x) * dx -
           (positions[1].y - positions[0].y) * dy) / (dx * dx + dy * dy)
  tx = positions[0].x - scale * scene[0].x
  ty = positions[0].y + scale * scene[0].y
  return lambda point: P(tx + scale * point.x, ty - scale * point.y), scale


class TestMechanismDrawing(unittest.TestCase):

  def setUp(self):
    self.carbonyl = molecule('C=O', [(0, 0), (1.5, 0)])
    self.nucleophile = molecule('[OH-]', [(0, 0)])

  def assertPoint(self, actual, expected, tolerance=0.11):
    self.assertAlmostEqual(actual.x, expected.x, delta=tolerance)
    self.assertAlmostEqual(actual.y, expected.y, delta=tolerance)

  def render(self, mols, arrows, size=(480, 320), **kwargs):
    drawer = rdMolDraw2D.MolDraw2DSVG(*size, -1, -1, True)
    positions = drawer.DrawMechanism(mols, arrows, **kwargs)
    drawer.FinishDrawing()
    return drawer, positions, ET.fromstring(drawer.GetDrawingText())

  def pointArrow(self, head=Head.Full):
    return Arrow(Anchor.Point(P(-1, 0)), Anchor.Point(P(2, 0)),
                 P(-0.5, 1), P(1.5, 1), head)

  def testFactoriesAndWritableFields(self):
    default = Anchor()
    self.assertEqual(default.type, rdMolDraw2D.MechanismAnchorType.Point)
    self.assertTrue(default.clipToLabel)
    atom = Anchor.Atom(fragment=1, atomIndex=2)
    self.assertEqual(atom.type, rdMolDraw2D.MechanismAnchorType.Atom)
    self.assertEqual((atom.fragment, atom.index), (1, 2))
    self.assertPoint(atom.offset, P(0, 0))
    bond = Anchor.Bond(2, 3, fraction=0.25, offset=P(0.1, 0.2))
    self.assertEqual(bond.type, rdMolDraw2D.MechanismAnchorType.Bond)
    self.assertEqual((bond.fragment, bond.index, bond.fraction), (2, 3, 0.25))
    pair = Anchor.LonePair(0, 1, P(0.4, 0.1), angle=30)
    self.assertEqual(pair.type, rdMolDraw2D.MechanismAnchorType.LonePair)
    self.assertEqual(pair.lonePairAngle, 30)
    self.assertAlmostEqual(pair.lonePairSeparation, 0.12)
    self.assertAlmostEqual(pair.lonePairRadius, 0.02)
    pair.lonePairAngle = 90
    pair.lonePairSeparation = 0.2
    pair.lonePairRadius = 0.05
    pair.clipToLabel = False
    pair.offset = P(0.5, 0.6)
    pair.fragment = 2
    pair.index = 3
    pair.fraction = 0.75
    pair.type = rdMolDraw2D.MechanismAnchorType.Atom
    self.assertEqual((pair.fragment, pair.index, pair.fraction), (2, 3, 0.75))
    self.assertEqual((pair.lonePairAngle, pair.lonePairSeparation,
                      pair.lonePairRadius), (90, 0.2, 0.05))
    self.assertFalse(pair.clipToLabel)
    self.assertPoint(pair.offset, P(0.5, 0.6))
    arrow = self.pointArrow()
    self.assertEqual(arrow.head, Head.Full)
    self.assertAlmostEqual(arrow.headLength, 0.16)
    self.assertAlmostEqual(arrow.headWidth, 0.06)
    self.assertAlmostEqual(arrow.lineWidth, 0.02)
    self.assertEqual(arrow.colour, (0, 0, 0, 1))
    arrow.start = atom
    arrow.end = bond
    arrow.control1 = P(1, 2)
    arrow.control2 = P(3, 4)
    arrow.head = Head.FishhookLeft
    arrow.headLength, arrow.headWidth, arrow.lineWidth = 0.4, 0.2, 0.08
    arrow.colour = (0.25, 0.5, 0.75, 0.5)
    self.assertEqual(arrow.colour, (0.25, 0.5, 0.75, 0.5))
    self.assertEqual(arrow.start.index, 2)
    self.assertEqual(arrow.end.index, 3)
    self.assertPoint(arrow.control1, P(1, 2))
    self.assertPoint(arrow.control2, P(3, 4))
    self.assertEqual(arrow.head, Head.FishhookLeft)
    self.assertEqual((arrow.headLength, arrow.headWidth, arrow.lineWidth),
                     (0.4, 0.2, 0.08))
    arrow.colour = (0.25, 0.5, 0.75)
    self.assertEqual(arrow.colour, (0.25, 0.5, 0.75, 1))
    for colour in [(1, 0), (0, 0, 0, 1, 1), (float('nan'), 0, 0),
                   (0, float('inf'), 0), (-0.1, 0, 0), (0, 0, 0, 1.1)]:
      with self.subTest(colour=colour), self.assertRaises(ValueError):
        arrow.colour = colour

  def testNativeCubicPrimitive(self):
    for raw in (False, True):
      with self.subTest(rawCoords=raw):
        drawer = rdMolDraw2D.MolDraw2DSVG(320, 240)
        drawer.SetScale(320, 240, P(-1, -1), P(3, 2))
        drawer.SetFillPolys(True)
        drawer.SetLineWidth(2.5)
        drawer.SetColour((0.25, 0.5, 0.75, 0.5))
        points = [P(10, 100), P(40, 20), P(140, 20), P(200, 100)] if raw else [
          P(0, 0), P(0.5, 1), P(1.5, 1), P(2, 0)]
        expected = points if raw else [drawer.GetDrawCoords(p) for p in points]
        drawer.DrawCubicBezier(*points, rawCoords=raw)
        drawer.FinishDrawing()
        paths = cubics(ET.fromstring(drawer.GetDrawingText()))
        self.assertEqual(len(paths), 1)
        self.assertEqual(paths[0].attrib['d'].count('C'), 1)
        self.assertNotIn('L', paths[0].attrib['d'])
        self.assertNotIn('Z', paths[0].attrib['d'])
        self.assertIn('fill:none', paths[0].attrib['style'])
        self.assertIn('stroke:#3F7FBF7F', paths[0].attrib['style'])
        self.assertIn('stroke-width:2.5px', paths[0].attrib['style'])
        for actual, target in zip(pathPoints(paths[0]), expected):
          self.assertPoint(actual, target)
        self.assertTrue(drawer.FillPolys())

  def testCubicRejectsNonfiniteCoordinates(self):
    backends = [rdMolDraw2D.MolDraw2DSVG]
    if hasattr(rdMolDraw2D, 'MolDraw2DCairo'):
      backends.append(rdMolDraw2D.MolDraw2DCairo)
    for backend in backends:
      for index in range(4):
        for value in (float('nan'), float('inf'), -float('inf')):
          for axis in ('x', 'y'):
            with self.subTest(backend=backend, index=index, value=value, axis=axis):
              points = [P(10, 10), P(20, 30), P(40, 30), P(50, 10)]
              setattr(points[index], axis, value)
              with self.assertRaises(ValueError):
                backend(200, 150).DrawCubicBezier(*points, rawCoords=True)

  def testCrossFragmentRotationOffsetsAndOriginalIndices(self):
    fragment = molecule('CCO', [(0, 0), (1.5, 0), (1.5, 1.5)])
    start = Anchor.Bond(0, 0, 0.25, P(0.3, 0.2))
    end = Anchor.Atom(1, 1, P(-0.2, 0.4))
    end.clipToLabel = False
    control1, control2 = P(0, 3), P(-2, 3)
    arrow = Arrow(start, end, control1, control2)
    _, positions, root = self.render(
      (fragment, self.carbonyl), (arrow,), offsets=(P(2, -1), P(-3, 2)),
      rotations=(90, -90), confIds=(-1, -1))
    self.assertIsInstance(positions, list)
    self.assertTrue(all(isinstance(group, list) for group in positions))
    self.assertEqual([len(group) for group in positions], [3, 2])
    self.assertIsInstance(positions[0][0], Geometry.Point2D)
    scene0 = [P(2, -1), P(2, 0.5), P(0.5, 0.5)]
    transform, scale = affineMap(positions[0], scene0)
    self.assertGreater(scale, 0)
    for actual, expected in zip(positions[0] + positions[1],
                                scene0 + [P(-3, 2), P(-3, 0.5)]):
      self.assertPoint(actual, transform(expected))
    curves = cubics(root)
    self.assertEqual(len(curves), 1)
    self.assertEqual(withClass(root, 'mechanism-shaft'), curves)
    self.assertIn(curves[0], withClass(root, 'mechanism-arrow-0'))
    expected = [P(1.8, -0.325), control1, control2, P(-2.6, 0.7)]
    for actual, target in zip(pathPoints(curves[0]), expected):
      self.assertPoint(actual, transform(target))
    atoms = [node for node in root.iter() if 'atom-1' in node.get('class', '').split()]
    bonds = [node for node in root.iter() if 'bond-0' in node.get('class', '').split()]
    self.assertTrue(atoms)
    self.assertTrue(bonds)

  def testMultipleArrowsAndFishhookSides(self):
    arrows = []
    colours = [(0.25, 0.5, 0.75), (0.75, 0.25, 0.5), (0.5, 0.75, 0.25)]
    hexColours = ['#3F7FBF', '#BF3F7F', '#7FBF3F']
    heads = [Head.Full, Head.FishhookLeft, Head.FishhookRight]
    for y, head, colour in zip((1, 2, 3), heads, colours):
      arrow = Arrow(Anchor.Point(P(-1, y)), Anchor.Point(P(2, y)),
                    P(-0.5, y + 0.6), P(1.5, y), head)
      arrow.colour = colour
      arrows.append(arrow)
    _, positions, root = self.render([self.carbonyl], arrows)
    transform, scale = affineMap(positions[0], [P(0, 0), P(1.5, 0)])
    curves = cubics(root)
    self.assertEqual(len(curves), 3)
    for index, (y, head, colour, curve) in enumerate(zip(
        (1, 2, 3), heads, hexColours, curves)):
      with self.subTest(head=head):
        tip = transform(P(2, y))
        self.assertPoint(pathPoints(curve)[-1], tip)
        self.assertIn(curve, withClass(root, 'mechanism-arrow-' + str(index)))
        headPaths = [path for path in withClass(root, 'mechanism-head')
                     if 'mechanism-arrow-' + str(index) in path.get('class', '').split()]
        self.assertTrue(all(colour in path.get('style', '') for path in headPaths))
        vertices = [p for path in headPaths for p in pathPoints(path)]
        self.assertTrue(vertices)
        low, high = min(p.y for p in vertices), max(p.y for p in vertices)
        if head != Head.FishhookRight:
          self.assertAlmostEqual(low, tip.y - arrows[0].headWidth * scale, delta=0.2)
        else:
          self.assertGreaterEqual(low, tip.y - 0.11)
        if head != Head.FishhookLeft:
          self.assertAlmostEqual(high, tip.y + arrows[0].headWidth * scale, delta=0.2)
        else:
          self.assertLessEqual(high, tip.y + 0.11)

  def testInputsAreNotModified(self):
    self.carbonyl.SetProp('mechanism_test', 'untouched')
    self.carbonyl.GetAtomWithIdx(0).SetAtomMapNum(17)
    self.carbonyl.GetAtomWithIdx(1).SetProp('source', 'oxygen')
    self.carbonyl.GetBondWithIdx(0).SetProp('source', 'pi')
    undepicted = Chem.MolFromSmiles('C1CC[C@H]2NCCCC2C1')
    inputs = [self.carbonyl, self.nucleophile, undepicted]
    before = [snapshot(mol) for mol in inputs]
    _, positions, _ = self.render(inputs, [self.pointArrow()],
                                   offsets=[P(0, 0), P(-3, 0), P(6, 0)])
    self.assertEqual([len(group) for group in positions],
                     [mol.GetNumAtoms() for mol in inputs])
    self.assertEqual([snapshot(mol) for mol in inputs], before)
    self.assertEqual(undepicted.GetNumConformers(), 0)

  def testConformerSelectionAndSequenceLifetimes(self):
    mol = molecule('CCC', [(0, 0), (1.5, 0), (1.5, 1.5)])
    conf = Chem.Conformer(mol.GetConformer())
    conf.SetId(7)
    conf.SetAtomPosition(2, Geometry.Point3D(1.5, -1.5, 0))
    mol.AddConformer(conf, assignId=False)
    _, positions, _ = self.render((mol,), (), confIds=(7,))
    self.assertGreater(positions[0][2].y, positions[0][1].y)
    self.assertEqual(mol.GetConformer(0).GetAtomPosition(2).y, 1.5)

    class MoleculeSequence:
      def __len__(self):
        return 2

      def __getitem__(self, index):
        if index >= 2:
          raise IndexError(index)
        return molecule('CC', [(0, 0), (1.5, 0)])

    _, positions, _ = self.render(MoleculeSequence(), (),
                                   offsets=(P(0, 0), P(3, 0)))
    self.assertEqual([len(group) for group in positions], [2, 2])

  def testAtomLabelClippingAndManualOverride(self):
    for clip in (True, False):
      with self.subTest(clipToLabel=clip):
        end = Anchor.Atom(0, 1)
        end.clipToLabel = clip
        arrow = Arrow(Anchor.Point(P(3, 0)), end, P(3, 1), P(2, 1))
        _, positions, root = self.render([self.carbonyl], [arrow])
        tip = pathPoints(cubics(root)[0])[-1]
        oxygen = positions[0][1]
        if clip:
          self.assertGreater(math.hypot(tip.x - oxygen.x, tip.y - oxygen.y), 1)
          self.assertGreater(tip.x, oxygen.x)
          self.assertLess(tip.y, oxygen.y)
        else:
          self.assertPoint(tip, oxygen)
    end.offset = P(0.25, 0.5)
    for clip in (False, True):
      end.clipToLabel = clip
      _, positions, root = self.render([self.carbonyl], [Arrow(
        Anchor.Point(P(3, 0)), end, P(3, 1), P(2, 1))])
      transform, _ = affineMap(positions[0], [P(0, 0), P(1.5, 0)])
      self.assertPoint(pathPoints(cubics(root)[0])[-1], transform(P(1.75, 0.5)))

  def testExplicitLonePairDotsRotateWithFragment(self):
    pair = Anchor.LonePair(0, 0, P(0.5, 0.4), angle=0)
    arrow = Arrow(pair, Anchor.Point(P(2, 2)), P(0, 2), P(1, 2))
    _, positions, root = self.render(
      [self.nucleophile, self.carbonyl], [arrow], rotations=[90, 0],
      offsets=[P(0, 0), P(3, 0)])
    transform, scale = affineMap(positions[1], [P(3, 0), P(4.5, 0)])
    dots = [node for node in root.iter()
            if node.tag in (SVG + 'ellipse', SVG + 'circle')]
    self.assertEqual(len(dots), 2)
    self.assertEqual(withClass(root, 'mechanism-lone-pair'), dots)
    centers = sorted([P(float(dot.get('cx')), float(dot.get('cy'))) for dot in dots],
                     key=lambda point: point.y)
    half_separation = pair.lonePairSeparation / 2
    for actual, scene in zip(centers, [P(-0.4, 0.5 + half_separation),
                                     P(-0.4, 0.5 - half_separation)]):
      self.assertPoint(actual, transform(scene))
    for dot in dots:
      self.assertAlmostEqual(float(dot.get('rx', dot.get('r'))),
                             pair.lonePairRadius * scale, delta=0.11)
    self.assertPoint(pathPoints(cubics(root)[0])[0], transform(P(-0.4, 0.5)))

    second = Arrow(pair, Anchor.Point(P(2, -2)), P(0, -2), P(1, -2))
    _, _, root = self.render([self.nucleophile], [arrow, second])
    self.assertEqual(len(withClass(root, 'mechanism-lone-pair')), 2)

  def testWholeSceneFitsFixedAndAutomaticCanvas(self):
    arrow = Arrow(Anchor.Point(P(-4, -2)), Anchor.Point(P(5, 3)),
                  P(-5, 6), P(6, 6))
    for size in ((480, 320), (-1, -1)):
      with self.subTest(size=size):
        drawer, positions, root = self.render([self.carbonyl], [arrow], size=size)
        width, height = drawer.Width(), drawer.Height()
        self.assertGreater(width, 0)
        self.assertGreater(height, 0)
        curve = pathPoints(cubics(root)[0])
        self.assertEqual(len(curve), 4)
        points = list(positions[0])
        for step in range(101):
          t, u = step / 100, 1 - step / 100
          points.append(P(
            u**3 * curve[0].x + 3*u*u*t * curve[1].x +
            3*u*t*t * curve[2].x + t**3 * curve[3].x,
            u**3 * curve[0].y + 3*u*u*t * curve[1].y +
            3*u*t*t * curve[2].y + t**3 * curve[3].y))
        for point in points:
          self.assertTrue(math.isfinite(point.x) and math.isfinite(point.y))
          self.assertGreaterEqual(point.x, 0)
          self.assertLessEqual(point.x, width)
          self.assertGreaterEqual(point.y, 0)
          self.assertLessEqual(point.y, height)

  def testOptionalSequencesAndDrawingStateRestored(self):
    _, default, _ = self.render([self.carbonyl], [self.pointArrow()])
    for optional in (None, []):
      _, positions, _ = self.render(
        [self.carbonyl], [self.pointArrow()], offsets=optional,
        rotations=optional, confIds=optional)
      for actual, expected in zip(positions[0], default[0]):
        self.assertPoint(actual, expected, 1e-8)
    drawer = rdMolDraw2D.MolDraw2DSVG(480, 320, -1, -1, True)
    options = drawer.drawOptions()
    values = dict(rotate=27, fixedFontSize=17, prepareMolsBeforeDrawing=False,
                  centreMoleculesBeforeDrawing=True, clearBackground=False,
                  scaleBondWidth=False)
    for name, value in values.items():
      setattr(options, name, value)
    drawer.SetLineWidth(3.5)
    drawer.SetFillPolys(True)
    drawer.SetColour((0.25, 0.5, 0.75, 0.5))
    drawer.DrawMechanism([self.carbonyl], [self.pointArrow()])
    invalid = self.pointArrow()
    invalid.control2 = P(2, 0)
    with self.assertRaises(ValueError):
      drawer.DrawMechanism([self.carbonyl], [invalid])
    self.assertEqual(drawer.LineWidth(), 3.5)
    self.assertTrue(drawer.FillPolys())
    for name, value in values.items():
      self.assertEqual(getattr(options, name), value)
    drawer.DrawCubicBezier(P(10, 10), P(20, 20), P(30, 20), P(40, 10),
                           rawCoords=True)
    drawer.FinishDrawing()
    last = cubics(ET.fromstring(drawer.GetDrawingText()))[-1]
    self.assertIn('stroke:#3F7FBF7F', last.attrib['style'])

  def testInvalidAnchorsAndGeometry(self):
    for factory in (Anchor.Atom, Anchor.Bond,
                    lambda fragment, index: Anchor.LonePair(fragment, index, P(0.5, 0))):
      for fragment, index in ((-1, 0), (0, -1)):
        with self.subTest(fragment=fragment, index=index), self.assertRaises(ValueError):
          factory(fragment, index)
    for field in ('fragment', 'index'):
      with self.subTest(field=field), self.assertRaises(ValueError):
        setattr(Anchor(), field, -1)
    for anchor in (Anchor.Atom(2, 0), Anchor.Atom(0, 10),
                   Anchor.Bond(0, 10), Anchor.LonePair(0, 10, P(0.5, 0))):
      with self.subTest(anchor=anchor), self.assertRaises(ValueError):
        self.render([self.carbonyl], [Arrow(
          anchor, Anchor.Point(P(3, 0)), P(0, 2), P(2, 2))])
    for fraction in (-0.1, 1.1, float('nan'), float('inf')):
      with self.subTest(fraction=fraction), self.assertRaises(ValueError):
        start = Anchor.Bond(0, 0, fraction)
        self.render([self.carbonyl], [Arrow(
          start, Anchor.Point(P(3, 0)), P(0, 2), P(2, 2))])
    for field in ('control1', 'control2'):
      for value in (float('nan'), float('inf')):
        with self.subTest(field=field, value=value), self.assertRaises(ValueError):
          arrow = self.pointArrow()
          setattr(arrow, field, P(value, 0))
          self.render([self.carbonyl], [arrow])
    for field in ('headLength', 'headWidth', 'lineWidth'):
      for value in (float('nan'), float('inf'), -0.1):
        with self.subTest(field=field, value=value), self.assertRaises(ValueError):
          arrow = self.pointArrow()
          setattr(arrow, field, value)
          self.render([self.carbonyl], [arrow])
    for factory in (lambda: Anchor.Atom(0, 0, P(float('nan'), 0)),
                    lambda: Anchor.Point(P(0, float('inf'))),
                    lambda: Anchor.LonePair(0, 0, P(0.5, 0), float('nan'))):
      with self.subTest(factory=factory), self.assertRaises(ValueError):
        self.render([self.carbonyl], [Arrow(
          factory(), Anchor.Point(P(3, 0)), P(0, 2), P(2, 2))])
    for field in ('lonePairRadius', 'lonePairSeparation'):
      for value in (0, -0.1, float('nan'), float('inf')):
        with self.subTest(field=field, value=value), self.assertRaises(ValueError):
          pair = Anchor.LonePair(0, 1, P(0.5, 0))
          setattr(pair, field, value)
          self.render([self.carbonyl], [Arrow(
            pair, Anchor.Point(P(3, 0)), P(0, 2), P(2, 2))])
    for field, endpoint in (('control1', P(-1, 0)), ('control2', P(2, 0))):
      arrow = self.pointArrow()
      setattr(arrow, field, endpoint)
      with self.subTest(field=field), self.assertRaises(ValueError):
        self.render([self.carbonyl], [arrow])

  def testInvalidPythonArguments(self):
    for mols in (None, 'CC', 1, [], [Chem.Mol()], [None],
                 [self.carbonyl, None], [P(0, 0)]):
      with self.subTest(mols=mols), self.assertRaises(ValueError):
        self.render(mols, [])
    for arrows in (None, 'arrow', [None], [P(0, 0)]):
      with self.subTest(arrows=arrows), self.assertRaises(ValueError):
        self.render([self.carbonyl], arrows)
    for kwargs in (
      {'offsets': [P(0, 0), P(1, 1)]}, {'rotations': [0, 0]},
      {'confIds': [-1, -1]}, {'offsets': [None]}, {'offsets': [(0, 0)]},
      {'offsets': 'invalid'}, {'rotations': ['90']}, {'confIds': [0.5]},
      {'offsets': [P(float('nan'), 0)]}, {'rotations': [float('inf')]},
    ):
      with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
        self.render([self.carbonyl], [], **kwargs)

  @unittest.skipUnless(hasattr(rdMolDraw2D, 'MolDraw2DCairo'), 'Cairo not enabled')
  def testCairoAndSvgPositionsAgree(self):
    positions = []
    for backend in (rdMolDraw2D.MolDraw2DSVG, rdMolDraw2D.MolDraw2DCairo):
      drawer = backend(480, 320, -1, -1, True)
      # Carbon-only labels avoid platform-specific font metrics in this comparison.
      mol = molecule('CCC', [(0, 0), (1.5, 0), (2.25, 1.3)])
      positions.append(drawer.DrawMechanism([mol], [self.pointArrow()]))
      drawer.FinishDrawing()
      output = drawer.GetDrawingText()
      if backend is rdMolDraw2D.MolDraw2DCairo:
        self.assertIsInstance(output, bytes)
        self.assertTrue(output.startswith(b'\x89PNG\r\n\x1a\n'))
        self.assertGreater(len(output), 100)
    for actual, expected in zip(positions[0][0], positions[1][0]):
      self.assertPoint(actual, expected, 0.5)


if __name__ == '__main__':
  unittest.main()
