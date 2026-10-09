//
// Copyright (C) 2026 RDKit contributors
// This file is part of the RDKit, covered by the BSD license in license.txt.
//

#ifndef RDKIT_MOLDRAW2DMECHANISM_H
#define RDKIT_MOLDRAW2DMECHANISM_H

#include <RDGeneral/export.h>
#include <Geometry/point.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>

namespace RDKit {

enum class MechanismAnchorType {
  Atom,
  Bond,
  LonePair,
  Point
};
enum class ElectronArrowHead {
  Full,
  FishhookLeft,
  FishhookRight
};

//! A user-specified electron-flow endpoint, not a chemical inference.
struct RDKIT_MOLDRAW2D_EXPORT MechanismAnchor {
  MechanismAnchorType type = MechanismAnchorType::Point;
  unsigned int fragment = 0;
  unsigned int index = 0;
  //! Atom/bond/lone-pair offsets use the fragment's unrotated coordinate frame.
  //! For Point anchors this is an absolute scene coordinate.
  RDGeom::Point2D offset;
  //! Position along a bond, from its begin atom (0) to its end atom (1).
  double fraction = 0.5;
  //! Trim zero-offset atom endpoints to labels, toward the adjacent control
  //! point. Nonzero offsets are explicit manual placements and are never
  //! trimmed.
  bool clipToLabel = true;
  //! Lone-pair dot-axis angle in degrees, counterclockwise in the fragment
  //! frame.
  double lonePairAngle = 0.0;
  double lonePairSeparation = 0.12;
  double lonePairRadius = 0.02;

  static MechanismAnchor atom(unsigned int fragment, unsigned int atomIndex,
                              const RDGeom::Point2D &offset = {});
  static MechanismAnchor bond(unsigned int fragment, unsigned int bondIndex,
                              double fraction = 0.5,
                              const RDGeom::Point2D &offset = {});
  //! Also draw two dots at the explicitly supplied position.
  static MechanismAnchor lonePair(unsigned int fragment, unsigned int atomIndex,
                                  const RDGeom::Point2D &offset,
                                  double angle = 0.0);
  static MechanismAnchor point(const RDGeom::Point2D &position);
};

//! A cubic electron-flow arrow. Control points and dimensions use scene units.
//! Arrowheads follow the final curve tangent. Fishhook sides are defined
//! looking along that tangent in the scene's Cartesian (positive Y up)
//! coordinate frame.
struct RDKIT_MOLDRAW2D_EXPORT MechanismArrow {
  MechanismAnchor start;
  MechanismAnchor end;
  RDGeom::Point2D control1;
  RDGeom::Point2D control2;
  ElectronArrowHead head = ElectronArrowHead::Full;
  double headLength = 0.16;
  //! Half-width of a full arrowhead.
  double headWidth = 0.06;
  double lineWidth = 0.02;
  DrawColour colour{0.0, 0.0, 0.0, 1.0};

  MechanismArrow(const MechanismAnchor &start, const MechanismAnchor &end,
                 const RDGeom::Point2D &control1,
                 const RDGeom::Point2D &control2,
                 ElectronArrowHead head = ElectronArrowHead::Full)
      : start(start),
        end(end),
        control1(control1),
        control2(control2),
        head(head) {}
};

}  // namespace RDKit

#endif
