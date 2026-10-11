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
#include <utility>

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

MechanismAnchor MechanismAnchor::lonePair(unsigned int fragment,
                                          unsigned int atomIndex) {
  auto result = lonePair(fragment, atomIndex, {});
  result.autoPosition = true;
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
  int bond = -1;
  bool clip = false;
  bool adjustBond = false;
  bool lonePair = false;
  bool autoPosition = false;
  Point2D dotAxis;
  double dotRadius = 0.0;
};

struct ResolvedArrow {
  ResolvedAnchor start;
  ResolvedAnchor end;
  MechanismArrow specification;
};

struct Bounds {
  double xmin = std::numeric_limits<double>::max();
  double xmax = std::numeric_limits<double>::lowest();
  double ymin = std::numeric_limits<double>::max();
  double ymax = std::numeric_limits<double>::lowest();

  void include(const Point2D &point) {
    xmin = std::min(xmin, point.x);
    xmax = std::max(xmax, point.x);
    ymin = std::min(ymin, point.y);
    ymax = std::max(ymax, point.y);
  }
  void pad(double padding) {
    xmin -= padding;
    xmax += padding;
    ymin -= padding;
    ymax += padding;
  }
  bool overlaps(const Bounds &other) const {
    return xmin < other.xmax && other.xmin < xmax && ymin < other.ymax &&
           other.ymin < ymax;
  }
  Point2D center() const { return {(xmin + xmax) / 2, (ymin + ymax) / 2}; }
  double distance(const Point2D &point) const {
    return std::hypot(std::max({xmin - point.x, 0.0, point.x - xmax}),
                      std::max({ymin - point.y, 0.0, point.y - ymax}));
  }
};

struct FragmentPlacement {
  Point2D offset;
  double rotation = 0.0;
};

Point2D placedPoint(const Point2D &point, const FragmentPlacement &placement) {
  return rotatePoint(point, placement.rotation) + placement.offset;
}

Point2D atomPosition(const ROMol &mol, unsigned int atom) {
  const auto &position = mol.getConformer().getAtomPos(atom);
  return {position.x, position.y};
}

Point2D anchorPosition(const MechanismAnchor &anchor, const ROMol &mol) {
  if (anchor.type == MechanismAnchorType::Bond) {
    const auto *bond = mol.getBondWithIdx(anchor.index);
    return atomPosition(mol, bond->getBeginAtomIdx()) * (1 - anchor.fraction) +
           atomPosition(mol, bond->getEndAtomIdx()) * anchor.fraction +
           anchor.offset;
  }
  return atomPosition(mol, anchor.index) + anchor.offset;
}

Bounds fragmentBounds(const ROMol &mol, const FragmentPlacement &placement,
                      double padding) {
  Bounds bounds;
  for (const auto &point : mol.getConformer().getPositions()) {
    bounds.include(placedPoint({point.x, point.y}, placement));
  }
  bounds.pad(padding);
  return bounds;
}

// Rigid placement keeps internal depictions and stereochemistry intact. Each
// candidate puts a linked endpoint near an already placed endpoint, then scores
// the entire set of connections, not just the edge used to propose it.
std::vector<FragmentPlacement> arrangeFragments(
    const std::vector<std::unique_ptr<RWMol>> &fragments,
    const std::vector<MechanismArrow> &arrows, double padding) {
  const auto count = fragments.size();
  std::vector<FragmentPlacement> placements(count);
  std::vector<bool> placed(count, false);
  std::vector<Bounds> boxes(count);
  std::vector<unsigned int> degrees(count, 0);
  for (const auto &arrow : arrows) {
    if (arrow.start.fragment != arrow.end.fragment) {
      ++degrees[arrow.start.fragment];
      ++degrees[arrow.end.fragment];
    }
  }
  auto priority = [&](size_t i) {
    return std::make_pair(degrees[i], fragments[i]->getNumAtoms());
  };
  size_t root = 0;
  for (size_t i = 1; i < count; ++i) {
    if (priority(i) > priority(root)) {
      root = i;
    }
  }
  placements[root].offset =
      fragmentBounds(*fragments[root], {}, 0).center() * -1;
  boxes[root] = fragmentBounds(*fragments[root], placements[root], padding);
  placed[root] = true;

  struct Link {
    const MechanismAnchor *local;
    const MechanismAnchor *remote;
  };
  for (size_t step = 1; step < count; ++step) {
    size_t next = count;
    std::vector<Link> links;
    for (size_t i = 0; i < count; ++i) {
      if (placed[i]) {
        continue;
      }
      std::vector<Link> candidateLinks;
      for (const auto &arrow : arrows) {
        if (arrow.start.fragment == i && placed[arrow.end.fragment]) {
          candidateLinks.push_back({&arrow.start, &arrow.end});
        } else if (arrow.end.fragment == i && placed[arrow.start.fragment]) {
          candidateLinks.push_back({&arrow.end, &arrow.start});
        }
      }
      if (next == count || candidateLinks.size() > links.size() ||
          (candidateLinks.size() == links.size() &&
           priority(i) > priority(next))) {
        next = i;
        links = std::move(candidateLinks);
      }
    }
    auto scorePlacement = [&](const FragmentPlacement &candidate,
                              const Bounds &bounds) {
      Bounds total = bounds;
      for (size_t i = 0; i < count; ++i) {
        if (!placed[i]) {
          continue;
        }
        if (bounds.overlaps(boxes[i])) {
          return std::numeric_limits<double>::infinity();
        }
        total.include({boxes[i].xmin, boxes[i].ymin});
        total.include({boxes[i].xmax, boxes[i].ymax});
      }
      double score =
          0.03 * (total.xmax - total.xmin) * (total.ymax - total.ymin);
      for (const auto &link : links) {
        const auto remoteIndex = link.remote->fragment;
        const auto local = placedPoint(
            anchorPosition(*link.local, *fragments[next]), candidate);
        const auto remote =
            placedPoint(anchorPosition(*link.remote, *fragments[remoteIndex]),
                        placements[remoteIndex]);
        const auto separation = remote - local;
        const double length = std::hypot(separation.x, separation.y);
        score += length * length;
        if (length <= 1.e-8) {
          continue;
        }
        // Prefer reactive sites on the facing sides of their fragments.
        const auto direction = separation / length;
        score += 2.0 *
                 std::max(0.0, (bounds.center() - local).dotProduct(direction));
        score += 2.0 * std::max(0.0, (boxes[remoteIndex].center() - remote)
                                         .dotProduct(direction * -1));
      }
      return score;
    };
    double bestScore = std::numeric_limits<double>::infinity();
    FragmentPlacement best;
    Bounds bestBounds;
    if (links.empty()) {
      double right = std::numeric_limits<double>::lowest();
      for (size_t i = 0; i < count; ++i) {
        if (placed[i]) {
          right = std::max(right, boxes[i].xmax);
        }
      }
      const auto local = fragmentBounds(*fragments[next], {}, padding);
      best.offset = {right + 1.0 - local.xmin, -local.center().y};
      bestBounds = fragmentBounds(*fragments[next], best, padding);
    } else {
      for (const auto &pivot : links) {
        const auto remoteIndex = pivot.remote->fragment;
        const auto target =
            placedPoint(anchorPosition(*pivot.remote, *fragments[remoteIndex]),
                        placements[remoteIndex]);
        const auto local = anchorPosition(*pivot.local, *fragments[next]);
        for (unsigned int rotation = 0; rotation < 8; ++rotation) {
          for (unsigned int side = 0; side < 16; ++side) {
            const auto direction = rotatePoint({1, 0}, side * 22.5);
            FragmentPlacement candidate;
            candidate.rotation = rotation * 45.0;
            // Increase the separation until the inflated molecular boxes fit.
            double separation = 1.6;
            for (unsigned int attempt = 0; attempt < 40; ++attempt) {
              candidate.offset = target + direction * separation -
                                 rotatePoint(local, candidate.rotation);
              const auto bounds =
                  fragmentBounds(*fragments[next], candidate, padding);
              const auto score = scorePlacement(candidate, bounds);
              if (std::isfinite(score)) {
                if (score < bestScore) {
                  bestScore = score;
                  best = candidate;
                  bestBounds = bounds;
                }
                break;
              }
              separation *= 1.2;
            }
          }
        }
      }
      if (!std::isfinite(bestScore)) {
        throw ValueErrorException("Unable to arrange mechanism fragments");
      }
    }
    placements[next] = best;
    boxes[next] = bestBounds;
    placed[next] = true;
  }
  return placements;
}

DrawColour automaticColour(size_t index) {
  static const std::array<DrawColour, 8> palette{{{0.0, 0.35, 0.65},
                                                  {0.8, 0.25, 0.0},
                                                  {0.0, 0.5, 0.3},
                                                  {0.55, 0.25, 0.65},
                                                  {0.65, 0.45, 0.0},
                                                  {0.0, 0.5, 0.65},
                                                  {0.75, 0.25, 0.45},
                                                  {0.35, 0.35, 0.35}}};
  if (index < palette.size()) {
    return palette[index];
  }
  const auto hue = std::fmod(index * 0.61803398875, 1.0) * 6;
  const auto x = 0.55 * (1 - std::abs(std::fmod(hue, 2.0) - 1));
  const std::array<DrawColour, 6> colours{{{0.7, x + 0.15, 0.15},
                                           {x + 0.15, 0.7, 0.15},
                                           {0.15, 0.7, x + 0.15},
                                           {0.15, x + 0.15, 0.7},
                                           {x + 0.15, 0.15, 0.7},
                                           {0.7, 0.15, x + 0.15}}};
  return colours[static_cast<size_t>(hue)];
}

double distanceToSegment(const Point2D &point, const Point2D &start,
                         const Point2D &end);

Point2D cubicPoint(const std::array<Point2D, 4> &curve, double t) {
  const double u = 1 - t;
  return curve[0] * (u * u * u) + curve[1] * (3 * u * u * t) +
         curve[2] * (3 * u * t * t) + curve[3] * (t * t * t);
}

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
    std::vector<std::array<Point2D, 4>> previousCurves;
    for (size_t i = 0; i < arrows_.size(); ++i) {
      auto arrow = arrows_[i];
      auto toFrame = [&](const Point2D &point) {
        return rotatePoint({point.x, -point.y}, drawOptions_.rotate);
      };
      for (auto *anchor : {&arrow.start, &arrow.end}) {
        anchor->position = toFrame(anchor->position);
        anchor->dotAxis = toFrame(anchor->dotAxis);
      }
      std::array<Point2D, 4> curve{
          arrow.start.position, toFrame(arrow.specification.control1),
          toFrame(arrow.specification.control2), arrow.end.position};
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
      if (arrow.specification.autoRoute) {
        curve = routeArrow(arrow, previousCurves);
      } else {
        if (arrow.start.autoPosition) {
          arrow.start = attach(arrow.start, unitTangent(curve[0], curve[1]));
          curve[0] = arrow.start.position;
        }
        if (arrow.end.autoPosition) {
          arrow.end = attach(arrow.end, unitTangent(curve[3], curve[2]));
          curve[3] = arrow.end.position;
        }
        clip(arrow.start, curve[0], curve[1]);
        clip(arrow.end, curve[3], curve[2]);
      }
      previousCurves.push_back(curve);
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
        const auto first = anchor->position - anchor->dotAxis;
        const auto second = anchor->position + anchor->dotAxis;
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
  Bounds labelBounds(int atom) const {
    Bounds bounds;
    atomLabels_[atom]->findExtremes(bounds.xmin, bounds.xmax, bounds.ymin,
                                    bounds.ymax);
    return bounds;
  }

  ResolvedAnchor attach(const ResolvedAnchor &source,
                        const Point2D &direction, bool bondTail = false) const {
    auto anchor = source;
    if (!anchor.autoPosition && !anchor.clip && !anchor.adjustBond) {
      return anchor;
    }
    if (bondTail && anchor.adjustBond) {
      const auto type = drawMol_->getBondWithIdx(anchor.bond)->getBondType();
      if (type == Bond::DOUBLE || type == Bond::TRIPLE) {
        double furthest = -1.0;
        for (const auto &shape : bonds_) {
          if (shape->bond_ != anchor.bond + activeBndIdxOffset_ ||
              dynamic_cast<const MolDraw2D_detail::DrawShapeSimpleLine *>(
                  shape.get()) == nullptr) {
            continue;
          }
          const auto segment = shape->points_[1] - shape->points_[0];
          const double denominator =
              direction.x * segment.y - direction.y * segment.x;
          if (std::abs(denominator) < 1.e-12) {
            continue;
          }
          const auto offset = shape->points_[0] - source.position;
          const double alongRay =
              (offset.x * segment.y - offset.y * segment.x) / denominator;
          const double alongStroke =
              (offset.x * direction.y - offset.y * direction.x) / denominator;
          if (alongRay >= -1.e-8 && alongStroke >= -1.e-8 &&
              alongStroke <= 1 + 1.e-8) {
            furthest = std::max(furthest, std::max(0.0, alongRay));
          }
        }
        if (furthest >= 0) {
          // Leave from the last rendered stroke, not between parallel lines.
          anchor.position += direction * furthest;
          requireFinite(anchor.position);
          return anchor;
        }
      }
    }
    double distance = 0.0;
    if (anchor.atom >= 0 && atomLabels_[anchor.atom]) {
      const auto bounds = labelBounds(anchor.atom);
      double enter = std::numeric_limits<double>::lowest();
      double exit = std::numeric_limits<double>::max();
      for (const auto &slab :
           {std::array<double, 4>{source.position.x, direction.x, bounds.xmin,
                                  bounds.xmax},
            std::array<double, 4>{source.position.y, direction.y, bounds.ymin,
                                  bounds.ymax}}) {
        if (std::abs(slab[1]) < 1.e-12) {
          if (slab[0] < slab[2] || slab[0] > slab[3]) {
            exit = -1;
            break;
          }
        } else {
          const double a = (slab[2] - slab[0]) / slab[1];
          const double b = (slab[3] - slab[0]) / slab[1];
          enter = std::max(enter, std::min(a, b));
          exit = std::min(exit, std::max(a, b));
        }
      }
      if (exit >= std::max(enter, 0.0)) {
        distance = exit;
      }
    }
    distance += 0.08 + std::max(0.0, drawOptions_.additionalAtomLabelPadding);
    if (anchor.autoPosition) {
      distance += anchor.dotRadius + anchor.dotAxis.length();
      anchor.dotAxis =
          Point2D{-direction.y, direction.x} * anchor.dotAxis.length();
    }
    anchor.position += direction * distance;
    requireFinite(anchor.position);
    return anchor;
  }

  double attachmentPenalty(const ResolvedAnchor &anchor,
                           const Point2D &direction) const {
    double score = 0.0;
    if (anchor.bond >= 0) {
      const auto *bond = drawMol_->getBondWithIdx(anchor.bond);
      const auto delta =
          atCds_[bond->getEndAtomIdx()] - atCds_[bond->getBeginAtomIdx()];
      if (delta.lengthSq() > 1.e-12) {
        score += 4 * std::pow(direction.dotProduct(delta / delta.length()), 4);
      }
    }
    if (anchor.atom >= 0) {
      const auto *atom = drawMol_->getAtomWithIdx(anchor.atom);
      for (const auto *neighbor : drawMol_->atomNeighbors(atom)) {
        const auto delta = atCds_[neighbor->getIdx()] - atCds_[anchor.atom];
        if (delta.lengthSq() > 1.e-12) {
          score += 3 * std::pow(std::max(0.0, direction.dotProduct(
                                                  delta / delta.length())),
                                4);
        }
      }
    }
    return score;
  }

  std::array<Point2D, 4> routeArrow(
      ResolvedArrow &arrow,
      const std::vector<std::array<Point2D, 4>> &previous) const {
    std::vector<std::pair<int, Bounds>> labels;
    for (size_t i = 0; i < atomLabels_.size(); ++i) {
      if (atomLabels_[i]) {
        labels.emplace_back(static_cast<int>(i),
                            labelBounds(static_cast<int>(i)));
      }
    }
    auto delta = arrow.end.position - arrow.start.position;
    const bool loop = delta.lengthSq() < 1.e-12;
    const auto forward =
        loop ? Point2D{1, 0}
             : unitTangent(arrow.start.position, arrow.end.position);
    double bestScore = std::numeric_limits<double>::infinity();
    std::array<Point2D, 4> best;
    ResolvedAnchor bestStart, bestEnd;
    for (unsigned int departure = 0; departure < 16; ++departure) {
      const auto startDirection = rotatePoint(forward, departure * 22.5);
      const auto start = attach(arrow.start, startDirection, true);
      for (unsigned int arrival = 0; arrival < 16; ++arrival) {
        const auto endDirection = rotatePoint(forward * -1, arrival * 22.5);
        if (loop && std::abs(startDirection.dotProduct(endDirection)) > 0.8) {
          continue;
        }
        const auto end = attach(arrow.end, endDirection);
        const auto separation =
            std::max(0.8, (end.position - start.position).length());
        for (double factor : {0.35, 0.65, 1.0}) {
          const double handle = separation * factor + 0.2;
          const std::array<Point2D, 4> curve{
              start.position, start.position + startDirection * handle,
              end.position + endDirection * handle, end.position};
          double score = attachmentPenalty(arrow.start, startDirection) +
                         attachmentPenalty(arrow.end, endDirection) +
                         0.8 * (2 - startDirection.dotProduct(forward) +
                                endDirection.dotProduct(forward));
          const auto chord = curve[3] - curve[0];
          const double chordLength = chord.length();
          if (chordLength > 1.e-8) {
            auto side = [&](const Point2D &point) {
              const auto offset = point - curve[0];
              return (chord.x * offset.y - chord.y * offset.x) / chordLength;
            };
            const double firstSide = side(curve[1]),
                         secondSide = side(curve[2]);
            if (std::max(std::abs(firstSide), std::abs(secondSide)) <
                0.12 * separation) {
              score += 4.0;
            }
            if (firstSide * secondSide < -1.e-6) {
              score += 2.0;
            }
          }
          auto obstaclePenalty = [&](const Point2D &point, double t) {
            double penalty = 0;
            for (const auto &label : labels) {
              // Explicit anchors may intentionally lie inside a label.
              if ((label.first == start.atom && !start.clip &&
                   !start.autoPosition && t < 0.12) ||
                  (label.first == end.atom && !end.clip && !end.autoPosition &&
                   t > 0.88)) {
                continue;
              }
              const double distance = label.second.distance(point);
              if (distance < 0.07) {
                penalty += distance == 0 ? 250.0 : 2.0 * (0.07 - distance);
              }
            }
            for (const auto *bond : drawMol_->bonds()) {
              const auto index = static_cast<int>(bond->getIdx());
              if ((index == start.bond && t < 0.12) ||
                  (index == end.bond && t > 0.88)) {
                continue;
              }
              const auto distance =
                  distanceToSegment(point, atCds_[bond->getBeginAtomIdx()],
                                    atCds_[bond->getEndAtomIdx()]);
              if (distance < 0.14) {
                penalty += 4.0 * (0.14 - distance) / 0.14;
              }
            }
            return penalty;
          };
          auto last = curve[0];
          for (unsigned int sample = 0; sample <= 24; ++sample) {
            const double t = sample / 24.0;
            const auto point = cubicPoint(curve, t);
            score += 0.12 * (point - last).length();
            score += obstaclePenalty(point, t);
            if (!loop && chordLength > 1.e-8) {
              const auto derivative =
                  (curve[1] - curve[0]) * ((1 - t) * (1 - t)) +
                  (curve[2] - curve[1]) * (2 * t * (1 - t)) +
                  (curve[3] - curve[2]) * (t * t);
              score += 4 * std::max(0.0, -derivative.dotProduct(chord) /
                                             (chordLength * chordLength));
              if (derivative.length() < 0.02 * chordLength) {
                score += 10.0;
              }
            }
            last = point;
            if (sample > 2 && sample < 22) {
              for (const auto &earlier : previous) {
                for (unsigned int k = 1; k < 8; ++k) {
                  if ((point - cubicPoint(earlier, k / 8.0)).lengthSq() <
                      0.01) {
                    score += 0.25;
                  }
                }
              }
            }
            if (score >= bestScore) {
              break;
            }
          }
          // Score the arrowhead and dots too, not just the shaft centerline.
          const auto tangent = endDirection * -1;
          const Point2D wing{tangent.y * arrow.specification.headWidth,
                             -tangent.x * arrow.specification.headWidth};
          const auto base =
              end.position - tangent * arrow.specification.headLength;
          if (arrow.specification.head != ElectronArrowHead::FishhookRight) {
            score += obstaclePenalty(base + wing, 1);
          }
          if (arrow.specification.head != ElectronArrowHead::FishhookLeft) {
            score += obstaclePenalty(base - wing, 1);
          }
          for (const auto *anchor : {&start, &end}) {
            if (anchor->lonePair) {
              score += obstaclePenalty(anchor->position - anchor->dotAxis,
                                       anchor == &start ? 0 : 1);
              score += obstaclePenalty(anchor->position + anchor->dotAxis,
                                       anchor == &start ? 0 : 1);
            }
          }
          if (score < bestScore) {
            bestScore = score;
            best = curve;
            bestStart = start;
            bestEnd = end;
          }
        }
      }
    }
    if (!std::isfinite(bestScore)) {
      throw ValueErrorException(
          "Unable to choose a finite mechanism arrow route");
    }
    arrow.start = bestStart;
    arrow.end = bestEnd;
    return best;
  }

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
    const std::vector<int> &confIds, bool autoLayout) {
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
  for (const auto *mol : mols) {
    if (!mol || !mol->getNumAtoms()) {
      throw ValueErrorException(
          "Mechanism fragments must be nonempty molecules");
    }
  }
  bool arrange =
      autoLayout && !arrows.empty() && offsets.empty() && rotations.empty();
  for (const auto &arrow : arrows) {
    arrange = arrange && arrow.autoRoute;
    for (const auto *anchor : {&arrow.start, &arrow.end}) {
      requireFinite(anchor->offset);
      if (anchor->type == MechanismAnchorType::Point) {
        arrange = false;
        continue;
      }
      if (anchor->fragment >= mols.size()) {
        throw ValueErrorException(
            "Mechanism anchor fragment index out of range");
      }
      const auto &mol = *mols[anchor->fragment];
      if (anchor->type == MechanismAnchorType::Bond) {
        if (anchor->index >= mol.getNumBonds()) {
          throw ValueErrorException("Mechanism anchor bond index out of range");
        }
        if (!std::isfinite(anchor->fraction) || anchor->fraction < 0 ||
            anchor->fraction > 1) {
          throw ValueErrorException(
              "Mechanism bond fraction must be between 0 and 1");
        }
      } else if (anchor->type == MechanismAnchorType::Atom ||
                 anchor->type == MechanismAnchorType::LonePair) {
        if (anchor->index >= mol.getNumAtoms()) {
          throw ValueErrorException("Mechanism anchor atom index out of range");
        }
      } else {
        throw ValueErrorException("Unknown mechanism anchor type");
      }
    }
  }
  RWMol scene;
  std::vector<unsigned int> atomOffsets;
  std::vector<std::unique_ptr<RWMol>> fragments;
  std::vector<FragmentPlacement> placements(mols.size());
  for (size_t i = 0; i < mols.size(); ++i) {
    const Point2D offset = offsets.empty() ? Point2D{} : offsets[i];
    const double rotation = rotations.empty() ? 0.0 : rotations[i];
    requireFinite(offset);
    if (!std::isfinite(rotation)) {
      throw ValueErrorException("Mechanism fragment rotations must be finite");
    }
    fragments.push_back(std::make_unique<RWMol>(*mols[i]));
    auto &fragment = *fragments.back();
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
    for (const auto &point : fragment.getConformer().getPositions()) {
      requireFinite({point.x, point.y});
    }
    placements[i] = {offset, rotation};
  }
  if (arrange) {
    double padding = 0.65;
    if (drawOptions().fixedFontSize > 0) {
      padding =
          std::max(padding, drawOptions().fixedFontSize * 0.65 /
                                std::max(20.0, drawOptions().scalingFactor));
    }
    placements = arrangeFragments(fragments, arrows, padding);
  }
  for (size_t i = 0; i < fragments.size(); ++i) {
    auto &fragment = *fragments[i];
    const auto &position = placements[i];
    for (auto &point : fragment.getConformer().getPositions()) {
      const auto placed = placedPoint({point.x, point.y}, position);
      requireFinite(placed);
      point = {placed.x, placed.y, 0.0};
    }
    RDGeom::Transform3D placement;
    placement.SetRotation(
        std::remainder(position.rotation, 360.0) * M_PI / 180.0,
        RDGeom::Point3D{0.0, 0.0, 1.0});
    placement.SetTranslation({position.offset.x, position.offset.y, 0.0});
    MolTransforms::transformMolSubstanceGroups(fragment, placement);
    atomOffsets.push_back(scene.getNumAtoms());
    scene.insertMol(fragment);
  }
  MolOps::symmetrizeSSSR(scene);
  const auto &conf = scene.getConformer();
  auto resolve = [&](const MechanismAnchor &anchor) {
    ResolvedAnchor result;
    if (anchor.type == MechanismAnchorType::Point) {
      result.position = anchor.offset;
      return result;
    }
    const auto &mol = *mols[anchor.fragment];
    const auto first = atomOffsets[anchor.fragment];
    const double rotation = placements[anchor.fragment].rotation;
    auto atomPoint = [&](unsigned int atom) {
      const auto &point = conf.getAtomPos(first + atom);
      return Point2D{point.x, point.y};
    };
    if (anchor.type == MechanismAnchorType::Atom ||
        anchor.type == MechanismAnchorType::LonePair) {
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
        result.autoPosition = anchor.autoPosition;
        result.dotRadius = anchor.lonePairRadius;
        result.dotAxis = rotatePoint({anchor.lonePairSeparation / 2, 0.0},
                                     anchor.lonePairAngle + rotation);
      }
    } else if (anchor.type == MechanismAnchorType::Bond) {
      const auto *bond = mol.getBondWithIdx(anchor.index);
      result.bond = scene
                        .getBondBetweenAtoms(first + bond->getBeginAtomIdx(),
                                             first + bond->getEndAtomIdx())
                        ->getIdx();
      result.adjustBond = anchor.offset.lengthSq() == 0;
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
    if (!arrow.autoRoute) {
      requireFinite(arrow.control1);
      requireFinite(arrow.control2);
    }
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
    if (!arrow.autoRoute) {
      unitTangent(start.position, arrow.control1);
      unitTangent(arrow.control2, end.position);
    }
    auto specification = arrow;
    if (arrow.autoColour) {
      specification.colour = automaticColour(resolved.size());
    }
    resolved.push_back({start, end, specification});
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
