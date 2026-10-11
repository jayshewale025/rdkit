"""Run installed examples without relying on a source checkout."""

import argparse
import json
from pathlib import Path
import subprocess
import sys


def doctor():
    import numpy
    import rdkit
    from rdkit import Chem, RDConfig
    from rdkit.Chem import AllChem, Descriptors
    from rdkit.Chem.Draw import rdMolDraw2D as draw

    anchor, arrow = draw.MechanismAnchor, draw.MechanismArrow
    mols = [Chem.MolFromSmiles("CC=O"), Chem.MolFromSmiles("[OH-]")]
    arrows = [
        arrow(anchor.LonePair(1, 0), anchor.Atom(0, 1)),
        arrow(anchor.Bond(0, 1), anchor.Atom(0, 2)),
    ]
    for backend in (draw.MolDraw2DSVG, draw.MolDraw2DCairo):
        drawer = backend(640, 320)
        positions = drawer.DrawMechanism(mols, arrows)
        drawer.FinishDrawing()
        assert [len(group) for group in positions] == [3, 1]
        output = drawer.GetDrawingText()
        if isinstance(output, str):
            assert output.count("mechanism-shaft") == 2
        else:
            assert output.startswith(b"\x89PNG\r\n\x1a\n")
    assert all(mol.GetNumConformers() == 0 for mol in mols)
    assert Descriptors.MolWt(Chem.MolFromSmiles("CCO")) > 46
    probe = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert AllChem.EmbedMolecule(probe, randomSeed=17) == 0
    assert (Path(RDConfig.RDDataDir) / "BaseFeatures.fdef").is_file()
    print(json.dumps({
        "status": "ok",
        "python": sys.version.split()[0],
        "rdkit": rdkit.__version__,
        "rdkit_location": rdkit.__file__,
        "data": RDConfig.RDDataDir,
        "numpy": numpy.__version__,
        "mechanism_svg": True,
        "mechanism_png": True,
    }, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("doctor", "demo", "gallery", "question6"))
    args, rest = parser.parse_known_args()
    if args.command == "doctor":
        if rest:
            parser.error("doctor takes no additional arguments")
        doctor()
        return 0
    scripts = {
        "demo": "mechanism_quickstart.py",
        "gallery": "mechanism_automatic.py",
        "question6": "mechanism_question6.py",
    }
    script = Path(__file__).parent / "examples" / scripts[args.command]
    return subprocess.call([sys.executable, str(script), *rest])


if __name__ == "__main__":
    raise SystemExit(main())
