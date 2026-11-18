// Copyright (C) 2026 RDKit contributors
// This file is part of the RDKit, covered by the BSD license in license.txt.

#include "minilib.h"

#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>
#include <boost/json.hpp>

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace bj = boost::json;
using namespace RDKit;

[[noreturn]] void invalid(const std::string &message) {
  throw std::invalid_argument("Mechanism: " + message);
}

const bj::object &object(const bj::value &value, const char *name) {
  if (!value.is_object()) {
    invalid(std::string(name) + " must be an object");
  }
  return value.as_object();
}

const bj::value &required(const bj::object &obj, const char *key) {
  const auto *value = obj.if_contains(key);
  if (!value) {
    invalid(std::string("missing '") + key + "'");
  }
  return *value;
}

std::string string(const bj::value &value, const char *name) {
  if (!value.is_string()) {
    invalid(std::string(name) + " must be a string");
  }
  return std::string(value.as_string());
}

double number(const bj::value &value, const char *name) {
  if (!value.is_number()) {
    invalid(std::string(name) + " must be a finite number");
  }
  const auto result = bj::value_to<double>(value);
  if (!std::isfinite(result)) {
    invalid(std::string(name) + " must be a finite number");
  }
  return result;
}

unsigned int index(const bj::value &value, const char *name) {
  const auto result = number(value, name);
  if (result < 0 || result > std::numeric_limits<unsigned int>::max() ||
      result != std::floor(result)) {
    invalid(std::string(name) + " must be a non-negative integer");
  }
  return static_cast<unsigned int>(result);
}

bool boolean(const bj::object &obj, const char *key, bool fallback) {
  const auto *value = obj.if_contains(key);
  if (!value) {
    return fallback;
  }
  if (!value->is_bool()) {
    invalid(std::string(key) + " must be a boolean");
  }
  return value->as_bool();
}

RDGeom::Point2D point(const bj::value &value, const char *name) {
  if (!value.is_array() || value.as_array().size() != 2) {
    invalid(std::string(name) + " must be [x, y]");
  }
  return {number(value.as_array()[0], name),
          number(value.as_array()[1], name)};
}

void optionalNumber(const bj::object &obj, const char *key, double &target) {
  if (const auto *value = obj.if_contains(key)) {
    target = number(*value, key);
  }
}

MechanismAnchor anchor(const bj::value &value,
                       const std::vector<std::unique_ptr<RWMol>> &mols) {
  const auto &obj = object(value, "anchor");
  const auto type = string(required(obj, "type"), "anchor.type");
  if (type == "point") {
    return MechanismAnchor::point(point(required(obj, "position"), "position"));
  }
  if (type != "atom" && type != "bond" && type != "lonePair") {
    invalid("anchor.type must be atom, bond, lonePair or point");
  }
  const auto fragment = index(required(obj, "fragment"), "fragment");
  const auto idx = index(required(obj, "index"), "index");
  if (fragment >= mols.size()) {
    invalid("fragment is out of range");
  }
  const auto count = type == "bond" ? mols[fragment]->getNumBonds()
                                     : mols[fragment]->getNumAtoms();
  if (idx >= count) {
    invalid(type + " index is out of range");
  }
  const auto *offsetValue = obj.if_contains("offset");
  const auto offset = offsetValue ? point(*offsetValue, "offset") : RDGeom::Point2D{};
  auto result = type == "atom" ? MechanismAnchor::atom(fragment, idx, offset)
                : type == "bond" ? MechanismAnchor::bond(fragment, idx, 0.5, offset)
                : offsetValue ? MechanismAnchor::lonePair(fragment, idx, offset)
                              : MechanismAnchor::lonePair(fragment, idx);
  result.clipToLabel = boolean(obj, "clipToLabel", result.clipToLabel);
  result.autoPosition = boolean(obj, "autoPosition", result.autoPosition);
  optionalNumber(obj, "fraction", result.fraction);
  optionalNumber(obj, "lonePairAngle", result.lonePairAngle);
  optionalNumber(obj, "lonePairRadius", result.lonePairRadius);
  optionalNumber(obj, "lonePairSeparation", result.lonePairSeparation);
  return result;
}
}  // namespace

std::string get_mechanism_svg(const std::string &spec_json,
                              const std::string &details_json) {
  using namespace RDKit;
  // Keep the parsed value alive while references into it are used.
  const auto spec = bj::parse(spec_json);
  const auto &obj = object(spec, "specification");
  const auto &inputs = required(obj, "molecules");
  if (!inputs.is_array() || inputs.as_array().empty()) {
    invalid("molecules must be a non-empty array of SMILES strings");
  }
  const auto &arrowInputs = required(obj, "arrows");
  if (!arrowInputs.is_array()) {
    invalid("arrows must be an array");
  }
  SmilesParserParams params;
  params.removeHs = false;  // moving beta-H atoms must retain their indices
  const bool kekulize = boolean(obj, "kekulize", false);
  const bool showAtomMaps = boolean(obj, "showAtomMaps", true);
  std::vector<std::unique_ptr<RWMol>> mols;
  std::vector<ROMol *> views;
  for (const auto &input : inputs.as_array()) {
    auto mol = std::unique_ptr<RWMol>(SmilesToMol(string(input, "molecule"), params));
    if (!mol) {
      invalid("invalid molecule at fragment " + std::to_string(mols.size()));
    }
    if (kekulize) {
      MolOps::Kekulize(*mol, true);
    }
    if (!showAtomMaps) {
      for (auto *atom : mol->atoms()) {
        atom->setAtomMapNum(0);
      }
    }
    views.push_back(mol.get());
    mols.push_back(std::move(mol));
  }
  std::vector<MechanismArrow> arrows;
  for (const auto &input : arrowInputs.as_array()) {
    const auto &entry = object(input, "arrow");
    const auto start = anchor(required(entry, "start"), mols);
    const auto end = anchor(required(entry, "end"), mols);
    const auto *headValue = entry.if_contains("head");
    const auto headName = headValue ? string(*headValue, "head") : "full";
    ElectronArrowHead head;
    if (headName == "full") {
      head = ElectronArrowHead::Full;
    } else if (headName == "fishhookLeft") {
      head = ElectronArrowHead::FishhookLeft;
    } else if (headName == "fishhookRight") {
      head = ElectronArrowHead::FishhookRight;
    } else {
      invalid("head must be full, fishhookLeft or fishhookRight");
    }
    MechanismArrow arrow(start, end, head);
    const auto *c1 = entry.if_contains("control1");
    const auto *c2 = entry.if_contains("control2");
    if (bool(c1) != bool(c2)) {
      invalid("manual curves need both control1 and control2");
    }
    if (c1) {
      arrow.control1 = point(*c1, "control1");
      arrow.control2 = point(*c2, "control2");
      arrow.autoRoute = false;
    }
    arrow.autoRoute = boolean(entry, "autoRoute", arrow.autoRoute);
    arrow.autoColour = boolean(entry, "autoColour", arrow.autoColour);
    optionalNumber(entry, "lineWidth", arrow.lineWidth);
    optionalNumber(entry, "headLength", arrow.headLength);
    optionalNumber(entry, "headWidth", arrow.headWidth);
    if (const auto *colour = entry.if_contains("colour")) {
      if (!colour->is_array() || colour->as_array().size() != 3) {
        invalid("colour must be [red, green, blue]");
      }
      double rgb[3];
      for (unsigned int i = 0; i < 3; ++i) {
        rgb[i] = number(colour->as_array()[i], "colour");
        if (rgb[i] < 0 || rgb[i] > 1) {
          invalid("colour components must be between 0 and 1");
        }
      }
      arrow.colour = DrawColour(rgb[0], rgb[1], rgb[2]);
      arrow.autoColour = false;
    }
    arrows.push_back(arrow);
  }

  const bj::value detailsValue = details_json.empty() ? bj::value(bj::object{})
                                                     : bj::parse(details_json);
  const auto &details = object(detailsValue, "details");
  auto dimension = [&](const char *name, int fallback) {
    const auto *value = details.if_contains(name);
    const auto size = value ? index(*value, name) : static_cast<unsigned int>(fallback);
    if (!size || size > 4096) {
      invalid(std::string(name) + " must be an integer from 1 to 4096");
    }
    return static_cast<int>(size);
  };
  MolDraw2DSVG drawer(dimension("width", 640), dimension("height", 320), -1, -1,
                      boolean(details, "noFreetype", false));
  if (!details_json.empty()) {
    MolDraw2DUtils::updateDrawerParamsFromJSON(drawer, details_json);
  }
  std::vector<RDGeom::Point2D> offsets;
  if (const auto *values = obj.if_contains("offsets")) {
    if (!values->is_array()) {
      invalid("offsets must be an array of [x, y]");
    }
    for (const auto &value : values->as_array()) {
      offsets.push_back(point(value, "offset"));
    }
  }
  std::vector<double> rotations;
  if (const auto *values = obj.if_contains("rotations")) {
    if (!values->is_array()) {
      invalid("rotations must be an array of angles in degrees");
    }
    for (const auto &value : values->as_array()) {
      rotations.push_back(number(value, "rotation"));
    }
  }
  drawer.drawMechanism(views, arrows, offsets, rotations, {},
                       boolean(obj, "autoLayout", true));
  drawer.finishDrawing();
  return drawer.getDrawingText();
}
