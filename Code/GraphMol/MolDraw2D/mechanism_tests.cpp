//
// Copyright (C) 2026 RDKit contributors
// This file is part of the RDKit, covered by the BSD license in license.txt.
//

#include <catch2/catch_all.hpp>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Exceptions.h>

#ifdef RDK_BUILD_CAIRO_SUPPORT
#include <GraphMol/MolDraw2D/MolDraw2DCairo.h>
#endif

#include <array>
#include <cmath>
#include <limits>
#include <memory>

using namespace RDKit;
using Catch::Approx;

namespace {
std::unique_ptr<RWMol> fragment(const std::string &smiles,
                                const std::vector<Point2D> &positions) {
  std::unique_ptr<RWMol> mol(SmilesToMol(smiles));
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == positions.size());
  auto conf = std::make_unique<Conformer>(positions.size());
  conf->set3D(false);
  for (size_t i = 0; i < positions.size(); ++i) {
    conf->setAtomPos(i, {positions[i].x, positions[i].y, 0});
  }
  mol->addConformer(conf.release(), true);
  return mol;
}

void nearPoint(const Point2D &actual, const Point2D &expected) {
  CHECK(actual.x == Approx(expected.x).margin(1.e-6));
  CHECK(actual.y == Approx(expected.y).margin(1.e-6));
}

class RecordingSVG : public MolDraw2DSVG {
 public:
  RecordingSVG(int width, int height)
      : MolDraw2DSVG(width, height, -1, -1, true) {}
  std::vector<std::array<Point2D, 4>> curves;
  std::vector<std::vector<Point2D>> heads;
  std::vector<std::pair<Point2D, Point2D>> dots;
  std::vector<double> widths;

  void drawCubicBezier(const Point2D &start, const Point2D &control1,
                       const Point2D &control2, const Point2D &end,
                       bool rawCoords = false) override {
    REQUIRE(rawCoords);
    curves.push_back({start, control1, control2, end});
    widths.push_back(getDrawLineWidth());
    MolDraw2DSVG::drawCubicBezier(start, control1, control2, end, rawCoords);
  }
  void drawPolygon(const std::vector<Point2D> &points,
                   bool rawCoords = false) override {
    if (getActiveClass().find("mechanism-head") != std::string::npos) {
      REQUIRE(rawCoords);
      heads.push_back(points);
    }
    MolDraw2DSVG::drawPolygon(points, rawCoords);
  }
  void drawLine(const Point2D &start, const Point2D &end,
                bool rawCoords = false) override {
    if (getActiveClass().find("mechanism-head") != std::string::npos) {
      REQUIRE(rawCoords);
      heads.push_back({start, end});
    }
    MolDraw2DSVG::drawLine(start, end, rawCoords);
  }
  void drawEllipse(const Point2D &start, const Point2D &end,
                   bool rawCoords = false) override {
    if (getActiveClass().find("mechanism-lone-pair") != std::string::npos) {
      REQUIRE(rawCoords);
      dots.emplace_back(start, end);
    }
    MolDraw2DSVG::drawEllipse(start, end, rawCoords);
  }
};

MechanismArrow manualArrow() {
  return {MechanismAnchor::point({-1, 0}),
          MechanismAnchor::point({2, 0}),
          {-1, 2},
          {2, 2}};
}

Point2D evaluate(const std::array<Point2D, 4> &curve, double t) {
  const double u = 1 - t;
  return curve[0] * (u * u * u) + curve[1] * (3 * u * u * t) +
         curve[2] * (3 * u * t * t) + curve[3] * (t * t * t);
}
}  // namespace

TEST_CASE("Mechanism batch retains fragment identities and shared transforms",
          "[mechanism]") {
  auto carbonyl = fragment("C=O", {{0, 0}, {1.5, 0}});
  auto donor = fragment("[OH-]", {{0, 0}});
  donor->getAtomWithIdx(0)->setAtomMapNum(7);
  auto source = MechanismAnchor::lonePair(1, 0, {0.0, 0.5});
  auto target = MechanismAnchor::atom(0, 0, {0.0, 0.1});
  MechanismArrow attack(source, target, {-2, 2}, {0, 2});
  MechanismArrow shift(MechanismAnchor::bond(0, 0, 0.5, {0, 0.1}),
                       MechanismAnchor::atom(0, 1, {0, 0.4}), {0.75, 1.0},
                       {1.5, 1.0});
  RecordingSVG drawer(600, 300);
  drawer.setColour({0.3, 0.4, 0.6, 0.7});
  const auto beforeWidth = drawer.lineWidth();
  const auto beforeFill = drawer.fillPolys();
  const auto beforeColour = drawer.colour();
  drawer.setDash({2, 3});
  drawer.setActiveClass("custom");
  auto coords = drawer.drawMechanism({carbonyl.get(), donor.get()},
                                     {attack, shift}, {{0, 0}, {-2, 0}});
  REQUIRE(coords.size() == 2);
  REQUIRE(coords[0].size() == 2);
  REQUIRE(coords[1].size() == 1);
  REQUIRE(drawer.curves.size() == 2);
  REQUIRE(drawer.heads.size() == 2);
  REQUIRE(drawer.dots.size() == 2);
  const auto scale = (coords[0][1] - coords[0][0]).length() / 1.5;
  nearPoint(coords[1][0], coords[0][0] + Point2D{-2 * scale, 0});
  nearPoint(drawer.curves[0][0], coords[1][0] + Point2D{0, -0.5 * scale});
  nearPoint(drawer.curves[0][3], coords[0][0] + Point2D{0, -0.1 * scale});
  nearPoint(drawer.curves[1][0],
            (coords[0][0] + coords[0][1]) / 2 + Point2D{0, -0.1 * scale});
  CHECK(drawer.widths[0] == Approx(attack.lineWidth * scale));
  CHECK(drawer.lineWidth() == beforeWidth);
  CHECK(drawer.fillPolys() == beforeFill);
  CHECK(drawer.colour() == beforeColour);
  CHECK(drawer.dash() == DashPattern{2, 3});
  CHECK(drawer.getActiveClass() == "custom");
  CHECK(carbonyl->getNumAtoms() == 2);
  CHECK(donor->getNumAtoms() == 1);
  CHECK(donor->getAtomWithIdx(0)->getAtomMapNum() == 7);
  CHECK(carbonyl->getConformer().getAtomPos(1).x == 1.5);
  drawer.finishDrawing();
  CHECK(drawer.getDrawingText().find("mechanism-arrow-0") != std::string::npos);
  CHECK(drawer.getDrawingText().find(" C ") != std::string::npos);
}

TEST_CASE("Mechanism scene rotates controls and anchors consistently",
          "[mechanism]") {
  auto mol = fragment("CC", {{0, 0}, {1.5, 0}});
  auto arrow = MechanismArrow(MechanismAnchor::atom(0, 1, {0, 0.2}),
                              MechanismAnchor::point({2, 2}), {0, 3}, {2, 3});
  RecordingSVG drawer(500, 400);
  drawer.drawOptions().centreMoleculesBeforeDrawing = true;
  drawer.drawOptions().rotate = 90;
  const auto coords =
      drawer.drawMechanism({mol.get()}, {arrow}, {{1, 1}}, {90});
  REQUIRE(drawer.curves.size() == 1);
  const double scale = (coords[0][1] - coords[0][0]).length() / 1.5;
  // Fragment (0, 0.2) rotates left; global clockwise rotation points it up.
  nearPoint(drawer.curves[0][0], coords[0][1] + Point2D{0, -0.2 * scale});
  nearPoint(coords[0][1], coords[0][0] + Point2D{1.5 * scale, 0});
  nearPoint(drawer.curves[0][3], coords[0][0] + Point2D{scale, scale});
  CHECK(drawer.drawOptions().centreMoleculesBeforeDrawing);
}

TEST_CASE("Mechanism heads follow tangents with consistent fishhook handedness",
          "[mechanism]") {
  auto mol = fragment("CC", {{0, 0}, {1.5, 0}});
  auto arrow = manualArrow();
  for (auto head : {ElectronArrowHead::Full, ElectronArrowHead::FishhookLeft,
                    ElectronArrowHead::FishhookRight}) {
    RecordingSVG drawer(300, 250);
    arrow.head = head;
    const auto coords = drawer.drawMechanism({mol.get()}, {arrow});
    const double scale = (coords[0][1] - coords[0][0]).length() / 1.5;
    REQUIRE(drawer.heads.size() == 1);
    const auto &curve = drawer.curves[0];
    auto tangent = curve[3] - curve[2];
    tangent.normalize();
    const Point2D left{tangent.y, -tangent.x};
    const auto base = curve[3] - tangent * arrow.headLength * scale;
    const auto &points = drawer.heads[0];
    REQUIRE(points.size() == (head == ElectronArrowHead::Full ? 3 : 2));
    if (head != ElectronArrowHead::FishhookRight) {
      nearPoint(points.front(), base + left * arrow.headWidth * scale);
    }
    if (head != ElectronArrowHead::FishhookLeft) {
      nearPoint(points.back(), base - left * arrow.headWidth * scale);
    }
  }
}

TEST_CASE("Mechanism bounds include curve extrema heads and strokes",
          "[mechanism]") {
  auto mol = fragment("CC", {{0, 0}, {1.5, 0}});
  auto arrow = manualArrow();
  arrow.control1 = {-12, 14};
  arrow.control2 = {15, 14};
  arrow.headLength = 2.0;
  arrow.headWidth = 1.0;
  arrow.lineWidth = 0.3;
  for (auto dimensions :
       {std::pair{-1, -1}, {450, 160}, {-1, 200}, {450, -1}}) {
    RecordingSVG drawer(dimensions.first, dimensions.second);
    drawer.drawOptions().fixedFontSize = 18;
    drawer.drawMechanism({mol.get()}, {arrow});
    REQUIRE(drawer.width() > 0);
    REQUIRE(drawer.height() > 0);
    const double radius = drawer.widths[0] / 2;
    auto contained = [&](const Point2D &point) {
      CHECK(point.x - radius >= -0.1);
      CHECK(point.y - radius >= -0.1);
      CHECK(point.x + radius <= drawer.width() + 0.1);
      CHECK(point.y + radius <= drawer.height() + 0.1);
    };
    for (int i = 0; i <= 1000; ++i) {
      contained(evaluate(drawer.curves[0], i / 1000.0));
    }
    for (auto point : drawer.heads[0]) {
      contained(point);
    }
  }
}

TEST_CASE("Mechanism atom label clipping is optional and manual offsets win",
          "[mechanism]") {
  auto mol = fragment("N", {{0, 0}});
  auto arrow = MechanismArrow(MechanismAnchor::atom(0, 0),
                              MechanismAnchor::point({3, 0}), {0, 2}, {3, 2});
  RecordingSVG automatic(400, 250);
  auto coords = automatic.drawMechanism({mol.get()}, {arrow});
  CHECK(automatic.curves[0][0].y < coords[0][0].y);
  arrow.start.clipToLabel = false;
  RecordingSVG manual(400, 250);
  coords = manual.drawMechanism({mol.get()}, {arrow});
  nearPoint(manual.curves[0][0], coords[0][0]);
  arrow.start.clipToLabel = true;
  arrow.start.offset = {0.1, 0};
  RecordingSVG offset(400, 250);
  coords = offset.drawMechanism({mol.get()}, {arrow});
  CHECK(offset.curves[0][0].x > coords[0][0].x);
  CHECK(offset.curves[0][0].y == Approx(coords[0][0].y));
  arrow.start.offset = {};
  arrow.control1 = {0, 0.001};
  RecordingSVG invalid(400, 250);
  CHECK_THROWS_AS(invalid.drawMechanism({mol.get()}, {arrow}),
                  ValueErrorException);
  CHECK(invalid.curves.empty());
}

TEST_CASE("Mechanism repeated lone-pair anchors do not duplicate dots",
          "[mechanism]") {
  auto mol = fragment("N", {{0, 0}});
  MechanismArrow arrow(MechanismAnchor::lonePair(0, 0, {0, 0.6}),
                       MechanismAnchor::point({3, 0}), {0, 2}, {3, 2});
  auto second = arrow;
  second.end = MechanismAnchor::point({-3, 0});
  second.control2 = {-3, 2};
  RecordingSVG drawer(400, 250);
  drawer.drawMechanism({mol.get()}, {arrow, second});
  CHECK(drawer.dots.size() == 2);
}

TEST_CASE("Mechanism validation precedes rendering", "[mechanism]") {
  auto mol = fragment("CC", {{0, 0}, {1.5, 0}});
  RecordingSVG drawer(300, 250);
  const auto before = drawer.getDrawingText();
  const auto arrow = manualArrow();
  CHECK_THROWS_AS(drawer.drawMechanism({}, {arrow}), ValueErrorException);
  CHECK_THROWS_AS(drawer.drawMechanism({nullptr}, {arrow}),
                  ValueErrorException);
  CHECK_THROWS_AS(drawer.drawMechanism({mol.get()}, {arrow}, {{}, {}}),
                  ValueErrorException);
  for (auto anchor :
       {MechanismAnchor::atom(1, 0), MechanismAnchor::atom(0, 2),
        MechanismAnchor::bond(0, 1), MechanismAnchor::bond(0, 0, 1.1),
        MechanismAnchor::point({std::numeric_limits<double>::infinity(), 0})}) {
    auto invalid = arrow;
    invalid.start = anchor;
    CHECK_THROWS_AS(drawer.drawMechanism({mol.get()}, {invalid}),
                    ValueErrorException);
  }
  auto invalid = arrow;
  invalid.control2 = {2, 0};
  CHECK_THROWS_AS(drawer.drawMechanism({mol.get()}, {invalid}),
                  ValueErrorException);
  invalid = arrow;
  invalid.lineWidth = 0;
  CHECK_THROWS_AS(drawer.drawMechanism({mol.get()}, {invalid}),
                  ValueErrorException);
  invalid = arrow;
  invalid.colour.a = std::numeric_limits<double>::quiet_NaN();
  CHECK_THROWS_AS(drawer.drawMechanism({mol.get()}, {invalid}),
                  ValueErrorException);
  CHECK(drawer.getDrawingText() == before);
}

TEST_CASE("Cubic fallback preserves strokes and handles straight curves",
          "[mechanism]") {
  RecordingSVG drawer(400, 300);
  drawer.setScale(400, 300, {-1, -1}, {3, 2});
  drawer.setFillPolys(true);
  drawer.MolDraw2D::drawCubicBezier({10, 100}, {40, 20}, {140, 20}, {200, 100},
                                    true);
  drawer.MolDraw2D::drawCubicBezier({10, 150}, {50, 150}, {150, 150},
                                    {200, 150}, true);
  CHECK(drawer.fillPolys());
  drawer.finishDrawing();
  const auto svg = drawer.getDrawingText();
  CHECK(svg.find("fill:none") != std::string::npos);
  CHECK(svg.find("nan") == std::string::npos);
  CHECK(svg.find("inf") == std::string::npos);
}

#ifdef RDK_BUILD_CAIRO_SUPPORT
TEST_CASE("Mechanism SVG and Cairo share geometry", "[mechanism]") {
  auto mol = fragment("CC", {{0, 0}, {1.5, 0}});
  const auto arrow = manualArrow();
  RecordingSVG svg(300, 250);
  MolDraw2DCairo cairo(300, 250, -1, -1, true);
  const auto svgCoords = svg.drawMechanism({mol.get()}, {arrow});
  const auto cairoCoords = cairo.drawMechanism({mol.get()}, {arrow});
  for (size_t i = 0; i < svgCoords[0].size(); ++i) {
    nearPoint(cairoCoords[0][i], svgCoords[0][i]);
  }
  cairo.finishDrawing();
  const auto png = cairo.getDrawingText();
  REQUIRE(png.size() > 8);
  CHECK(png.substr(0, 8) == std::string("\x89PNG\r\n\x1a\n", 8));
}
#endif
