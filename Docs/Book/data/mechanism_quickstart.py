"""Draw carbonyl addition with the installed fork; no coordinates required."""

import argparse
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("mechanism-demo"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    anchor, arrow = rdMolDraw2D.MechanismAnchor, rdMolDraw2D.MechanismArrow
    aldehyde = Chem.MolFromSmiles("CC=O")
    hydroxide = Chem.MolFromSmiles("[OH-]")
    carbonyl = aldehyde.GetBondBetweenAtoms(1, 2)
    arrows = [
        arrow(anchor.LonePair(1, 0), anchor.Atom(0, 1)),
        arrow(anchor.Bond(0, carbonyl.GetIdx()), anchor.Atom(0, 2)),
    ]
    for extension, backend in (
        ("svg", rdMolDraw2D.MolDraw2DSVG),
        ("png", rdMolDraw2D.MolDraw2DCairo),
    ):
        drawer = backend(640, 320)
        drawer.drawOptions().fixedFontSize = 22
        drawer.drawOptions().fixedBondLength = 65
        drawer.drawOptions().useBWAtomPalette()
        drawer.DrawMechanism([aldehyde, hydroxide], arrows)
        drawer.FinishDrawing()
        data = drawer.GetDrawingText()
        destination = args.output / f"carbonyl-addition.{extension}"
        if extension == "svg":
            destination.write_text(data, encoding="utf-8")
        else:
            destination.write_bytes(data)
        print(destination.resolve())


if __name__ == "__main__":
    main()
