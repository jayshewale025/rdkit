//
// Copyright (C) 2026 RDKit contributors
// This file is part of the RDKit, covered by the BSD license in license.txt.
//

#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/DrawMol.h>
#include <GraphMol/MolDraw2D/DrawText.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <Geometry/Transform3D.h>
#include <RDGeneral/Exceptions.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>

namespace RDKit {

MechanismAnchor MechanismAnchor::atom(unsigned int fragment,
                                      unsigned int atomIndex,
                                      const Point2D &offset) {
  MechanismAnchor result;
  result.type = MechanismAnchorType::Atom;
  result.fragment = fragment;
  result.index = atomIndex;
  result.offset = offset;
  return result;
}

MechanismAnchor MechanismAnchor::bond(unsigned int fragment,
                                      unsigned int bondIndex, double fraction,
                                      const Point2D &offset) {
  auto result = atom(fragment, bondIndex, offset);
  result.type = MechanismAnchorType::Bond;
  result.fraction = fraction;
  return result;
}

MechanismAnchor MechanismAnchor::lonePair(unsigned int fragment,
                                          unsigned int atomIndex,
                                          const Point2D &offset, double angle) {
  auto result = atom(fragment, atomIndex, offset);
  result.type = MechanismAnchorType::LonePair;
  result.lonePairAngle = angle;
  return result;
}

MechanismAnchor MechanismAnchor::point(const Point2D &position) {
  MechanismAnchor result;
  result.offset = position;
  return result;
}

namespace {

void requireFinite(const Point2D &point) {
  if (!std::isfinite(point.x) || !std::isfinite(point.y)) {
    throw ValueErrorException("Mechanism coordinates must be finite");
  }
}

void requirePositive(double value, const char *message) {
  if (!std::isfinite(value) || value <= 0.0) {
    throw ValueErrorException(message);
  }
}

Point2D rotatePoint(const Point2D &point, double degrees) {
  const auto radians = std::remainder(degrees, 360.0) * M_PI / 180.0;
  return {std::cos(radians) * point.x - std::sin(radians) * point.y,
          std::sin(radians) * point.x + std::cos(radians) * point.y};
}

Point2D unitTangent(const Point2D &from, const Point2D &to) {
  auto delta = to - from;
  const auto length = std::hypot(delta.x, delta.y);
  if (!std::isfinite(length) || length <= 1.e-8) {
    throw ValueErrorException(
        "Mechanism control points must give nonzero endpoint tangents");
  }
  return delta / length;
}

struct ResolvedAnchor {
  Point2D position;
  int atom = -1;
  bool clip = false;
  bool lonePair = false;
  Point2D dotAxis;
  double dotRadius = 0.0;
};

struct ResolvedArrow {
  ResolvedAnchor start;
  ResolvedAnchor end;
  MechanismArrow specification;
};

// Mechanism styling must not leak into later molecule or primitive drawing.
class MechanismStyle {
 public:
  MechanismStyle(MolDraw2D &drawer, const std::string &className = {})
      : drawer_(drawer),
        oldClass_(drawer.getActiveClass()),
        oldDash_(drawer.dash()),
        oldColour_(drawer.colour()),
        oldWidth_(drawer.lineWidth()),
        oldFill_(drawer.fillPolys()),
        oldScaleWidth_(drawer.drawOptions().scaleBondWidth) {
    if (!className.empty()) {
      drawer_.setActiveClass(oldClass_.empty() ? className
                                               : oldClass_ + " " + className);
      drawer_.setDash({});
    }
  }
  ~MechanismStyle() {
    drawer_.setActiveClass(oldClass_);
    drawer_.setDash(oldDash_);
    drawer_.setColour(oldColour_);
    drawer_.setLineWidth(oldWidth_);
    drawer_.setFillPolys(oldFill_);
    drawer_.drawOptions().scaleBondWidth = oldScaleWidth_;
  }

 private:
  MolDraw2D &drawer_;
  std::string oldClass_;
  DashPattern oldDash_;
  DrawColour oldColour_;
  double oldWidth_;
  bool oldFill_, oldScaleWidth_;
};

class MechanismArrowShape : public MolDraw2D_detail::DrawShape {
 public:
  MechanismArrowShape(const std::array<Point2D, 4> &curve,
                      const MechanismArrow &arrow, size_t index,
                      const Point2D &leftWing, const Point2D &rightWing)
      : DrawShape(std::vector<Point2D>(curve.begin(), curve.end()),
                  arrow.lineWidth, false, arrow.colour, false),
        head_(arrow.head),
        index_(index) {
    points_.push_back(leftWing);
    points_.push_back(rightWing);
  }

  void myDraw(MolDraw2D &drawer) const override {
    const auto name = "mechanism-arrow-" + std::to_string(index_);
    {
      MechanismStyle style(drawer, name + " mechanism-shaft");
      drawer.drawCubicBezier(points_[0], points_[1], points_[2], points_[3],
                             true);
    }
    MechanismStyle style(drawer, name + " mechanism-head");
    std::vector<Point2D> head;
    if (head_ != ElectronArrowHead::FishhookRight) {
      head.push_back(points_[4]);
    }
    head.push_back(points_[3]);
    if (head_ != ElectronArrowHead::FishhookLeft) {
      head.push_back(points_[5]);
    }
    if (head.size() == 2) {
      drawer.drawLine(head[0], head[1], true);
    } else {
      drawer.drawPolygon(head, true);
    }
  }

  void scale(const Point2D &factor) override {
    DrawShape::scale(factor);
    lineWidth_ *= std::abs(factor.x);
  }

  void findExtremes(double &xmin, double &xmax, double &ymin,
                    double &ymax) const override {
    double xlo = std::min(points_[0].x, points_[3].x);
    double xhi = std::max(points_[0].x, points_[3].x);
    double ylo = std::min(points_[0].y, points_[3].y);
    double yhi = std::max(points_[0].y, points_[3].y);
    auto extrema = [](double p0, double p1, double p2, double p3, double &lo,
                      double &hi) {
      const double a = -p0 + 3 * p1 - 3 * p2 + p3;
      const double b = 2 * (p0 - 2 * p1 + p2);
      const double c = p1 - p0;
      auto include = [&](double t) {
        if (t > 0 && t < 1) {
          const double u = 1 - t;
          const double value = u * u * u * p0 + 3 * u * u * t * p1 +
                               3 * u * t * t * p2 + t * t * t * p3;
          lo = std::min(lo, value);
          hi = std::max(hi, value);
        }
      };
      if (std::abs(a) < 1.e-12) {
        if (std::abs(b) > 1.e-12) {
          include(-c / b);
        }
      } else {
        const double discriminant = b * b - 4 * a * c;
        if (discriminant >= 0) {
          const double q =
              -0.5 * (b + std::copysign(std::sqrt(discriminant), b));
          if (q == 0) {
            include(-b / (2 * a));
          } else {
            include(q / a);
            include(c / q);
          }
        }
      }
    };
    extrema(points_[0].x, points_[1].x, points_[2].x, points_[3].x, xlo, xhi);
    extrema(points_[0].y, points_[1].y, points_[2].y, points_[3].y, ylo, yhi);
    for (size_t i = 4; i < 6; ++i) {
      if ((i == 4 && head_ == ElectronArrowHead::FishhookRight) ||
          (i == 5 && head_ == ElectronArrowHead::FishhookLeft)) {
        continue;
      }
      xlo = std::min(xlo, points_[i].x);
      xhi = std::max(xhi, points_[i].x);
      ylo = std::min(ylo, points_[i].y);
      yhi = std::max(yhi, points_[i].y);
    }
    double strokePadding = lineWidth_ / 2;
    if (head_ == ElectronArrowHead::Full) {
      // Existing polygon backends may use a miter join at the V-shaped tip.
      const auto halfBase = (points_[4] - points_[5]).length() / 2;
      strokePadding *= (points_[4] - points_[3]).length() / halfBase;
    }
    xmin = std::min(xmin, xlo - strokePadding);
    xmax = std::max(xmax, xhi + strokePadding);
    ymin = std::min(ymin, ylo - strokePadding);
    ymax = std::max(ymax, yhi + strokePadding);
  }

 private:
  ElectronArrowHead head_;
  size_t index_;
};

class LonePairShape : public MolDraw2D_detail::DrawShape {
 public:
  LonePairShape(const Point2D &first, const Point2D &second, double radius,
                const DrawColour &colour)
      : DrawShape({first, second}, 0.0, false, colour, true), radius_(radius) {}

  void myDraw(MolDraw2D &drawer) const override {
    MechanismStyle style(drawer, "mechanism-lone-pair");
    const Point2D radius{radius_, radius_};
    for (const auto &point : points_) {
      drawer.drawEllipse(point - radius, point + radius, true);
    }
  }
  void scale(const Point2D &factor) override {
    DrawShape::scale(factor);
    radius_ *= std::abs(factor.x);
  }
  void findExtremes(double &xmin, double &xmax, double &ymin,
                    double &ymax) const override {
    for (const auto &point : points_) {
      xmin = std::min(xmin, point.x - radius_);
      xmax = std::max(xmax, point.x + radius_);
      ymin = std::min(ymin, point.y - radius_);
      ymax = std::max(ymax, point.y + radius_);
    }
  }
  bool samePair(const Point2D &first, const Point2D &second, double radius,
                const DrawColour &colour) const {
    return (points_[0] - first).lengthSq() < 1.e-16 &&
           (points_[1] - second).lengthSq() < 1.e-16 && radius_ == radius &&
           lineColour_ == colour;
  }

 private:
  double radius_;
};

class MechanismDrawMol : public MolDraw2D_detail::DrawMol {
 public:
  MechanismDrawMol(const ROMol &mol, int width, int height,
                   std::shared_ptr<MolDrawOptions> options,
                   MolDraw2D_detail::DrawText &text,
                   std::vector<ResolvedArrow> arrows)
      : DrawMol(mol, "", width, height, *options, text, nullptr, nullptr,
                nullptr, nullptr, nullptr, nullptr, true, 0, true),
        options_(std::move(options)),
        arrows_(std::move(arrows)) {}

  void extractAll(double scale) override {
    DrawMol::extractAll(scale);
    std::vector<LonePairShape *> pairs;
    for (size_t i = 0; i < arrows_.size(); ++i) {
      const auto &arrow = arrows_[i];
      auto toFrame = [&](const Point2D &point) {
        return rotatePoint({point.x, -point.y}, drawOptions_.rotate);
      };
      std::array<Point2D, 4> curve{
          toFrame(arrow.start.position), toFrame(arrow.specification.control1),
          toFrame(arrow.specification.control2), toFrame(arrow.end.position)};
      auto clip = [&](const ResolvedAnchor &anchor, Point2D &point,
                      const Point2D &control) {
        if (anchor.clip && atomLabels_[anchor.atom]) {
          const auto padding =
              0.033 * meanBondLength_ +
              std::max(0.0, drawOptions_.additionalAtomLabelPadding);
          const auto direction = unitTangent(point, control);
          double xmin = std::numeric_limits<double>::max();
          double xmax = std::numeric_limits<double>::lowest();
          double ymin = xmin, ymax = xmax;
          atomLabels_[anchor.atom]->findExtremes(xmin, xmax, ymin, ymax);
          const auto outside =
              point + direction * (std::hypot(xmax - xmin, ymax - ymin) +
                                   padding + 1.0);
          MolDraw2D_detail::adjustBondEndForString(
              outside, padding, atomLabels_[anchor.atom]->rects_, point);
          if ((control - point).dotProduct(direction) <= 1.e-8) {
            throw ValueErrorException(
                "Mechanism control point lies inside its atom label; move "
                "it outward or disable clipToLabel");
          }
        }
      };
      clip(arrow.start, curve[0], curve[1]);
      clip(arrow.end, curve[3], curve[2]);
      unitTangent(curve[0], curve[1]);
      const auto tangent = unitTangent(curve[2], curve[3]);
      // The drawing frame reverses handedness relative to Cartesian scene
      // units.
      const Point2D left{tangent.y, -tangent.x};
      const auto base = curve[3] - tangent * arrow.specification.headLength;
      const auto wing = left * arrow.specification.headWidth;
      requireFinite(base + wing);
      requireFinite(base - wing);
      postShapes_.push_back(std::make_unique<MechanismArrowShape>(
          curve, arrow.specification, i, base + wing, base - wing));

      for (const auto *anchor : {&arrow.start, &arrow.end}) {
        if (!anchor->lonePair) {
          continue;
        }
        const auto first = toFrame(anchor->position - anchor->dotAxis);
        const auto second = toFrame(anchor->position + anchor->dotAxis);
        requireFinite(first);
        requireFinite(second);
        const auto &colour = arrow.specification.colour;
        if (std::any_of(pairs.begin(), pairs.end(), [&](const auto *pair) {
              return pair->samePair(first, second, anchor->dotRadius, colour);
            })) {
          continue;
        }
        auto pair = std::make_unique<LonePairShape>(first, second,
                                                    anchor->dotRadius, colour);
        pairs.push_back(pair.get());
        postShapes_.push_back(std::move(pair));
      }
    }
  }

 private:
  // DrawMol holds a reference, so the scene's options snapshot must outlive it.
  std::shared_ptr<MolDrawOptions> options_;
  std::vector<ResolvedArrow> arrows_;
};

double distanceToSegment(const Point2D &point, const Point2D &start,
                         const Point2D &end) {
  const auto delta = end - start;
  const double length = delta.lengthSq();
  if (length < 1.e-16) {
    return (point - start).length();
  }
  const double t =
      std::clamp((point - start).dotProduct(delta) / length, 0.0, 1.0);
  return (point - (start + delta * t)).length();
}

void flattenCubic(const std::array<Point2D, 4> &curve,
                  std::vector<Point2D> &points, unsigned int depth = 0) {
  if (std::max(distanceToSegment(curve[1], curve[0], curve[3]),
               distanceToSegment(curve[2], curve[0], curve[3])) <= 0.05) {
    points.push_back(curve[3]);
    return;
  }
  if (depth == 16) {
    throw ValueErrorException(
        "Bezier curve exceeds subdivision accuracy limit");
  }
  const auto a = (curve[0] + curve[1]) / 2;
  const auto b = (curve[1] + curve[2]) / 2;
  const auto c = (curve[2] + curve[3]) / 2;
  const auto d = (a + b) / 2;
  const auto e = (b + c) / 2;
  const auto middle = (d + e) / 2;
  flattenCubic({curve[0], a, d, middle}, points, depth + 1);
  flattenCubic({middle, e, c, curve[3]}, points, depth + 1);
}

}  // namespace

void MolDraw2D::drawCubicBezier(const Point2D &start, const Point2D &control1,
                                const Point2D &control2, const Point2D &end,
                                bool rawCoords) {
  std::array<Point2D, 4> curve{start, control1, control2, end};
  for (auto &point : curve) {
    requireFinite(point);
    if (!rawCoords) {
      point = getDrawCoords(point);
      requireFinite(point);
    }
  }
  std::vector<Point2D> points{curve[0]};
  flattenCubic(curve, points);
  if (points.size() == 2) {
    drawLine(points[0], points[1], true);
    return;
  }
  struct FillGuard {
    MolDraw2D &drawer;
    bool fill;
    ~FillGuard() { drawer.setFillPolys(fill); }
  } guard{*this, fillPolys()};
  setFillPolys(false);
  drawPolygon(points, true);
}

std::vector<std::vector<Point2D>> MolDraw2D::drawMechanism(
    const std::vector<ROMol *> &mols, const std::vector<MechanismArrow> &arrows,
    const std::vector<Point2D> &offsets, const std::vector<double> &rotations,
    const std::vector<int> &confIds) {
  if (mols.empty()) {
    throw ValueErrorException("DrawMechanism requires at least one fragment");
  }
  if ((!offsets.empty() && offsets.size() != mols.size()) ||
      (!rotations.empty() && rotations.size() != mols.size()) ||
      (!confIds.empty() && confIds.size() != mols.size())) {
    throw ValueErrorException(
        "Mechanism placement lists must match fragment count");
  }
  if (panelWidth() == 0 || panelHeight() == 0) {
    throw ValueErrorException("Mechanism canvas dimensions cannot be zero");
  }
  if (!std::isfinite(drawOptions().rotate)) {
    throw ValueErrorException("Mechanism scene rotation must be finite");
  }
  MechanismStyle restoreStyle(*this);
  RWMol scene;
  std::vector<unsigned int> atomOffsets;
  for (size_t i = 0; i < mols.size(); ++i) {
    if (!mols[i] || !mols[i]->getNumAtoms()) {
      throw ValueErrorException(
          "Mechanism fragments must be nonempty molecules");
    }
    const Point2D offset = offsets.empty() ? Point2D{} : offsets[i];
    const double rotation = rotations.empty() ? 0.0 : rotations[i];
    requireFinite(offset);
    if (!std::isfinite(rotation)) {
      throw ValueErrorException("Mechanism fragment rotations must be finite");
    }
    RWMol fragment(*mols[i]);
    const int confId = confIds.empty() ? -1 : confIds[i];
    if (fragment.getNumConformers()) {
      auto conf = std::make_unique<Conformer>(fragment.getConformer(confId));
      fragment.clearConformers();
      conf->setId(0);
      conf->set3D(false);
      fragment.addConformer(conf.release(), false);
    } else if (confId != -1) {
      throw ValueErrorException(
          "Conformer ID given for a fragment without coordinates");
    }
    if (!fragment.getRingInfo()->isSymmSssr()) {
      MolOps::symmetrizeSSSR(fragment);
    }
    if (drawOptions().prepareMolsBeforeDrawing) {
      MolDraw2DUtils::prepareMolForDrawing(fragment);
    } else if (!fragment.getNumConformers()) {
      RDDepict::compute2DCoords(fragment);
    }
    auto &conf = fragment.getConformer();
    for (auto &point : conf.getPositions()) {
      requireFinite({point.x, point.y});
      const auto placed = rotatePoint({point.x, point.y}, rotation) + offset;
      requireFinite(placed);
      point = {placed.x, placed.y, 0.0};
    }
    RDGeom::Transform3D placement;
    placement.SetRotation(std::remainder(rotation, 360.0) * M_PI / 180.0,
                          RDGeom::Point3D{0.0, 0.0, 1.0});
    placement.SetTranslation({offset.x, offset.y, 0.0});
    MolTransforms::transformMolSubstanceGroups(fragment, placement);
    atomOffsets.push_back(scene.getNumAtoms());
    scene.insertMol(fragment);
  }
  MolOps::symmetrizeSSSR(scene);
  const auto &conf = scene.getConformer();
  auto resolve = [&](const MechanismAnchor &anchor) {
    ResolvedAnchor result;
    requireFinite(anchor.offset);
    if (anchor.type == MechanismAnchorType::Point) {
      result.position = anchor.offset;
      return result;
    }
    if (anchor.fragment >= mols.size()) {
      throw ValueErrorException("Mechanism anchor fragment index out of range");
    }
    const auto &mol = *mols[anchor.fragment];
    const auto first = atomOffsets[anchor.fragment];
    const double rotation =
        rotations.empty() ? 0.0 : rotations[anchor.fragment];
    auto atomPoint = [&](unsigned int atom) {
      const auto &point = conf.getAtomPos(first + atom);
      return Point2D{point.x, point.y};
    };
    if (anchor.type == MechanismAnchorType::Atom ||
        anchor.type == MechanismAnchorType::LonePair) {
      if (anchor.index >= mol.getNumAtoms()) {
        throw ValueErrorException("Mechanism anchor atom index out of range");
      }
      result.position = atomPoint(anchor.index);
      result.atom = first + anchor.index;
      result.clip = anchor.type == MechanismAnchorType::Atom &&
                    anchor.clipToLabel && anchor.offset.lengthSq() == 0.0;
      if (anchor.type == MechanismAnchorType::LonePair) {
        requirePositive(anchor.lonePairRadius,
                        "Lone-pair radius must be positive");
        requirePositive(anchor.lonePairSeparation,
                        "Lone-pair separation must be positive");
        if (!std::isfinite(anchor.lonePairAngle) ||
            anchor.lonePairSeparation <= 2 * anchor.lonePairRadius) {
          throw ValueErrorException(
              "Lone-pair dots must be separated and angle finite");
        }
        result.lonePair = true;
        result.dotRadius = anchor.lonePairRadius;
        result.dotAxis = rotatePoint({anchor.lonePairSeparation / 2, 0.0},
                                     anchor.lonePairAngle + rotation);
      }
    } else if (anchor.type == MechanismAnchorType::Bond) {
      if (anchor.index >= mol.getNumBonds()) {
        throw ValueErrorException("Mechanism anchor bond index out of range");
      }
      if (!std::isfinite(anchor.fraction) || anchor.fraction < 0 ||
          anchor.fraction > 1) {
        throw ValueErrorException(
            "Mechanism bond fraction must be between 0 and 1");
      }
      const auto *bond = mol.getBondWithIdx(anchor.index);
      result.position =
          atomPoint(bond->getBeginAtomIdx()) * (1 - anchor.fraction) +
          atomPoint(bond->getEndAtomIdx()) * anchor.fraction;
    } else {
      throw ValueErrorException("Unknown mechanism anchor type");
    }
    result.position += rotatePoint(anchor.offset, rotation);
    requireFinite(result.position);
    return result;
  };
  std::vector<ResolvedArrow> resolved;
  for (const auto &arrow : arrows) {
    requireFinite(arrow.control1);
    requireFinite(arrow.control2);
    requirePositive(arrow.lineWidth, "Mechanism line width must be positive");
    requirePositive(arrow.headLength, "Mechanism head length must be positive");
    requirePositive(arrow.headWidth, "Mechanism head width must be positive");
    for (double channel :
         {arrow.colour.r, arrow.colour.g, arrow.colour.b, arrow.colour.a}) {
      if (!std::isfinite(channel) || channel < 0 || channel > 1) {
        throw ValueErrorException(
            "Mechanism colour channels must be between 0 and 1");
      }
    }
    if (arrow.head != ElectronArrowHead::Full &&
        arrow.head != ElectronArrowHead::FishhookLeft &&
        arrow.head != ElectronArrowHead::FishhookRight) {
      throw ValueErrorException("Unknown mechanism arrowhead type");
    }
    auto start = resolve(arrow.start);
    auto end = resolve(arrow.end);
    unitTangent(start.position, arrow.control1);
    unitTangent(arrow.control2, end.position);
    resolved.push_back({start, end, arrow});
  }

  auto options = std::make_shared<MolDrawOptions>(drawOptions());
  options->centreMoleculesBeforeDrawing = false;
  options->prepareMolsBeforeDrawing = false;
  options->drawingExtentsInclude |= DrawElement::POSTSHAPES;
  setupTextDrawer();
  auto drawing = std::make_shared<MechanismDrawMol>(
      scene, panelWidth(), panelHeight(), options, *text_drawer_,
      std::move(resolved));
  drawing->setOffsets(x_offset_, y_offset_);
  drawing->createDrawObjects();
  fixVariableDimensions(*drawing);
  drawMols_.push_back(drawing);
  activeMolIdx_ = static_cast<int>(drawMols_.size()) - 1;
  startDrawing();
  drawTheMolecule(*drawing);

  std::vector<std::vector<Point2D>> positions(mols.size());
  for (size_t i = 0; i < mols.size(); ++i) {
    for (unsigned int atom = 0; atom < mols[i]->getNumAtoms(); ++atom) {
      positions[i].push_back(drawing->getDrawCoords(atomOffsets[i] + atom));
    }
  }
  return positions;
}

}  // namespace RDKit
