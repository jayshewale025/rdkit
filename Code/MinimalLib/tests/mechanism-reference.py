"""Draw the same fixtures through Python for comparison with the WASM outputs.

Requires the mechanism fork, not stock RDKit. Run after mechanism.js.
"""

import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D as draw

HERE = Path(__file__).resolve().parent
OUT = HERE / "mechanism-output"


def render(spec):
    params = Chem.SmilesParserParams()
    params.removeHs = False
    mols = [Chem.MolFromSmiles(smiles, params) for smiles in spec["molecules"]]
    for mol in mols:
        assert mol is not None
        if spec.get("kekulize"):
            Chem.Kekulize(mol, clearAromaticFlags=True)
    A = draw.MechanismAnchor

    def anchor(value):
        kind = value["type"]
        if kind == "point":
            from rdkit.Geometry import Point2D

            return A.Point(Point2D(*value["position"]))
        fn = {"atom": A.Atom, "bond": A.Bond, "lonePair": A.LonePair}[kind]
        result = fn(value["fragment"], value["index"])
        if kind == "lonePair":
            result.lonePairRadius = 0.045
            result.lonePairSeparation = 0.17
        return result

    from rdkit.Geometry import Point2D

    heads = {
        "full": draw.ElectronArrowHead.Full,
        "fishhookLeft": draw.ElectronArrowHead.FishhookLeft,
        "fishhookRight": draw.ElectronArrowHead.FishhookRight,
    }
    arrows = []
    for entry in spec["arrows"]:
        arrow = draw.MechanismArrow(anchor(entry["start"]), anchor(entry["end"]),
                                    head=heads[entry.get("head", "full")])
        arrow.lineWidth, arrow.headLength, arrow.headWidth = 0.05, 0.26, 0.1
        if "control1" in entry:
            arrow.control1 = Point2D(*entry["control1"])
            arrow.control2 = Point2D(*entry["control2"])
            arrow.autoRoute = False
        if "colour" in entry:
            arrow.colour = tuple(entry["colour"])
            arrow.autoColour = False
        arrows.append(arrow)
    drawer = draw.MolDraw2DSVG(640, 320)
    opts = drawer.drawOptions()
    opts.fixedBondLength = 34
    opts.bondLineWidth = 1.6
    opts.clearBackground = False
    opts.padding = 0.06
    opts.fontFile = "BuiltinTelexRegular"
    drawer.DrawMechanism(
        mols, arrows,
        offsets=[Point2D(*p) for p in spec.get("offsets", [])],
        rotations=spec.get("rotations", []),
        autoLayout=spec.get("autoLayout", True),
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def main():
    assert hasattr(draw.MolDraw2DSVG, "DrawMechanism"), "Install the mechanism fork's wheel"
    OUT.mkdir(exist_ok=True)
    for spec in json.loads((HERE / "mechanism-cases.json").read_text()):
        (OUT / f"python-{spec['name']}.svg").write_text(render(spec))
        print(f"Python reference: {spec['name']}")


if __name__ == "__main__":
    main()
