"""Render the anthrone synthesis with the installed native mechanism API."""

import argparse
import html
import io
import json
from pathlib import Path
import sys
import textwrap
from collections import Counter

from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem, RDConfig
from rdkit.Chem.Draw import rdMolDraw2D as D
from mechanism_common import Scene, render_scene

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=Path("question-6"))
parser.add_argument("--compact", action="store_true",
                    help="Render a one-page textbook scheme instead of 16 panels")
args = parser.parse_args()
OUT = args.output
A, Arrow = D.MechanismAnchor, D.MechanismArrow


def molecule(smiles):
    params = Chem.SmilesParserParams()
    params.removeHs = False
    mol = Chem.MolFromSmiles(smiles, params)
    if mol is None:
        raise ValueError(f"Invalid mechanism structure: {smiles}")
    Chem.Kekulize(mol, clearAromaticFlags=True)
    return mol


def index(mol, number):
    matches = [atom.GetIdx() for atom in mol.GetAtoms()
               if atom.GetAtomMapNum() == number]
    if len(matches) != 1:
        raise ValueError(f"Atom map {number} is not unique")
    return matches[0]


def atom(mol, number, fragment=0):
    return A.Atom(fragment, index(mol, number))


def pair(mol, number, fragment=0):
    return A.LonePair(fragment, index(mol, number))


def bond(mol, first, second, fragment=0):
    found = mol.GetBondBetweenAtoms(index(mol, first), index(mol, second))
    if found is None:
        raise ValueError(f"Missing bond {first}-{second}")
    return A.Bond(fragment, found.GetIdx())


def changed(mol, bonds=(), charges=(), remove=()):
    rw = Chem.RWMol(mol)
    for first, second, order in bonds:
        i, j = index(rw, first), index(rw, second)
        existing = rw.GetBondBetweenAtoms(i, j)
        if order is None:
            rw.RemoveBond(i, j)
        elif existing is None:
            rw.AddBond(i, j, order)
        else:
            existing.SetBondType(order)
    for number, charge in charges:
        rw.GetAtomWithIdx(index(rw, number)).SetFormalCharge(charge)
    for number in remove:
        rw.RemoveAtom(index(rw, number))
    result = rw.GetMol()
    Chem.SanitizeMol(result)
    Chem.Kekulize(result, clearAromaticFlags=True)
    return result


def clean(mol):
    result = Chem.Mol(mol)
    for item in result.GetAtoms():
        item.SetAtomMapNum(0)
    return result


panels = []


def panel(slug, title, description, mols, arrows=(), outcome=None):
    assert all(mol.GetNumConformers() == 0 for mol in mols)
    assert all(arrow.autoRoute for arrow in arrows)
    panels.append((Scene(slug, title, description, [clean(m) for m in mols],
                         list(arrows), width=1000, height=450), outcome))


S, DB = Chem.BondType.SINGLE, Chem.BondType.DOUBLE
phthalic = molecule(
    "[O:2]=[C:1]1[O:3][C:4](=[O:5])[c:6]2[cH:7][cH:8][cH:9][cH:10][c:11]21")
alcl3 = molecule("[Al:20]([Cl:21])([Cl:22])[Cl:23]")
activated = changed(
    Chem.CombineMols(phthalic, alcl3), [(3, 20, S)], [(3, 1), (20, -1)])
open_acyl = changed(activated, [(1, 3, None)], [(1, 1), (3, 0)])
panel("01-activation", "1. Lewis-acid activation of phthalic anhydride",
      "An anhydride oxygen donates to AlCl3. This drawing selects the "
      "bridging-oxygen adduct as a useful Lewis-acid activation model.",
      [phthalic, alcl3], [Arrow(pair(phthalic, 3), atom(alcl3, 20, 1))],
      activated)
panel("02-opening", "2. Open the activated anhydride",
      "The acyl-O bond pair moves to oxygen, giving an acylium/aluminate "
      "ion pair. O=C+ is one resonance contributor of the acylium group.",
      [activated], [Arrow(bond(activated, 1, 3), atom(activated, 3))],
      open_acyl)

benzene = molecule("[c:31]1([H:37])[cH:32][cH:33][cH:34][cH:35][cH:36]1")
attack_bond = next(b for b in benzene.GetAtomWithIdx(index(benzene, 31)).GetBonds()
                   if b.GetBondType() == DB)
neighbor = attack_bond.GetOtherAtomIdx(index(benzene, 31))
positive = benzene.GetAtomWithIdx(neighbor).GetAtomMapNum()
arenium = changed(Chem.CombineMols(open_acyl, benzene),
                  [(1, 31, S), (31, positive, S)], [(1, 0), (positive, 1)])
panel("03-friedel-crafts", "3. Benzene attacks the acylium carbon",
      "A benzene pi pair forms the new C-C bond. The product is a sigma "
      "complex (arenium ion), not an aromatic substitution product yet.",
      [open_acyl, benzene],
      [Arrow(bond(benzene, 31, positive, 1), atom(open_acyl, 1))], arenium)

acid_complex = changed(arenium, [(31, 37, None), (3, 37, S),
                                (31, positive, DB)], [(positive, 0), (3, 1)])
panel("04-rearomatization", "4. Lose H+ and restore the benzene ring",
      "The aluminate-bound carboxylate oxygen is shown as the proton "
      "acceptor. C-H electrons restore the missing aromatic pi bond. "
      "This is formal proton-transfer bookkeeping; a solvent/counterion "
      "relay can replace the direct transfer.",
      [arenium], [Arrow(pair(arenium, 3), atom(arenium, 37)),
                   Arrow(bond(arenium, 31, 37), bond(arenium, 31, positive))],
      acid_complex)

acid_a = molecule("[O:2]=[C:1]([c:31]1[cH:32][cH:33][cH:34][cH:35][cH:36]1)"
                  "[c:11]1[c:6]([C:4](=[O:5])[O:3][H:37])"
                  "[cH:7][cH:8][cH:9][cH:10]1")
acid_b = molecule("[CH2:1]([c:31]1[cH:32][cH:33][cH:34][cH:35][cH:36]1)"
                  "[c:11]1[c:6]([C:4](=[O:5])[O:3][H:37])"
                  "[cH:7][cH:8][cH:9][cH:10]1")
panel("05-A", "5. H3O+ workup releases A: 2-benzoylbenzoic acid",
      "Water destroys the aluminium complexes and protonates the carboxylate. "
      "The isolated product has an ortho-COOH group and a ketone. "
      "Detailed aluminium hydrolysis is omitted.", [acid_a])
panel("06-clemmensen", "6. Zn-Hg / HCl, heat: A becomes B",
      "Clemmensen reduction replaces the ketone C=O by CH2; COOH remains. "
      "Net carbonyl reduction: R2C=O + 4 H+ + 4 e- -> R2CH2 + H2O. "
      "This is a heterogeneous metal-surface step, not an asserted sequence "
      "of solution-phase two-electron arrows.", [acid_a], outcome=acid_b)

thionyl = molecule("[O:42]=[S:40]([Cl:41])[Cl:43]")
oxysulfonium = changed(Chem.CombineMols(acid_b, thionyl),
                      [(3, 40, S), (40, 41, None)], [(3, 1), (41, -1)])
panel("07-thionyl", "7. B + SOCl2: activate the carboxylic acid",
      "One standard arrow-pushing representation uses the acid OH oxygen "
      "to attack sulfur as chloride leaves. Proton transfer gives an acyl "
      "chlorosulfite. Carbonyl-oxygen attack/proton-transfer representations "
      "are also possible.", [acid_b, thionyl],
      [Arrow(pair(acid_b, 3), atom(thionyl, 40, 1)),
       Arrow(bond(thionyl, 40, 41, 1), atom(thionyl, 41, 1))], oxysulfonium)

chlorosulfite_hcl = changed(oxysulfonium, [(3, 37, None), (41, 37, S)],
                          [(3, 0), (41, 0)])
panel("08-deprotonation", "8. Form the neutral acyl chlorosulfite",
      "Chloride removes the OH proton, giving HCl and R-C(=O)-O-S(=O)Cl.",
      [oxysulfonium],
      [Arrow(pair(oxysulfonium, 41), atom(oxysulfonium, 37)),
       Arrow(bond(oxysulfonium, 3, 37), atom(oxysulfonium, 3))],
      chlorosulfite_hcl)
chlorosulfite = changed(chlorosulfite_hcl, remove=(37, 41))
chloride = molecule("[Cl-:44]")
tetrahedral = changed(Chem.CombineMols(chlorosulfite, chloride),
                      [(4, 44, S), (4, 5, S)], [(44, 0), (5, -1)])
panel("09-chloride-attack", "9. Chloride attacks the activated acyl carbon",
      "Chloride in the reaction medium attacks C=O. The pi pair moves to "
      "oxygen, giving the tetrahedral acyl-substitution intermediate. "
      "HCl and other spectator species are omitted here.",
      [chlorosulfite, chloride],
      [Arrow(pair(chloride, 44, 1), atom(chlorosulfite, 4)),
       Arrow(bond(chlorosulfite, 4, 5), atom(chlorosulfite, 5))], tetrahedral)

chlorosulfite_anion = changed(tetrahedral, [(4, 5, DB), (4, 3, None)],
                             [(5, 0), (3, -1)])
panel("10-collapse", "10. Collapse to the acid chloride",
      "Oxygen reforms C=O and the C-O bond breaks. The departing "
      "chlorosulfite anion subsequently fragments into SO2 and chloride.",
      [tetrahedral],
      [Arrow(pair(tetrahedral, 5), bond(tetrahedral, 4, 5)),
       Arrow(bond(tetrahedral, 4, 3), atom(tetrahedral, 3))],
      chlorosulfite_anion)
leaving = molecule("[O-:3][S:40](=[O:42])[Cl:43]")
panel("11-gas-loss", "11. Leaving-group fragmentation releases SO2",
      "A lone pair forms the second S=O bond while chloride leaves. "
      "SO2 and HCl formation drive the overall acid-to-acid-chloride conversion.",
      [leaving], [Arrow(pair(leaving, 3), bond(leaving, 3, 40)),
                  Arrow(bond(leaving, 40, 43), atom(leaving, 43))],
      molecule("[O:3]=[S:40]=[O:42].[Cl-:43]"))

acyl_chloride = changed(chlorosulfite_anion, remove=(3, 40, 42, 43))
alcl3_2 = molecule("[Al:50]([Cl:51])([Cl:52])[Cl:53]")
acyl_chloride_adduct = changed(Chem.CombineMols(acyl_chloride, alcl3_2),
                               [(44, 50, S)], [(44, 1), (50, -1)])
panel("12-acid-chloride-activation", "12. AlCl3 activates the acid chloride",
      "Chloride donates to AlCl3. Heterolytic C-Cl cleavage then generates "
      "the intramolecular acylating electrophile.",
      [acyl_chloride, alcl3_2],
      [Arrow(pair(acyl_chloride, 44), atom(alcl3_2, 50, 1))],
      acyl_chloride_adduct)
acylium_with_counterion = changed(acyl_chloride_adduct, [(4, 44, None)],
                                  [(4, 1), (44, 0)])
panel("13-acylium", "13. Generate the tethered acylium ion",
      "The acyl C-Cl bond pair moves to chloride. AlCl4- is the counterion.",
      [acyl_chloride_adduct],
      [Arrow(bond(acyl_chloride_adduct, 4, 44), atom(acyl_chloride_adduct, 44))],
      acylium_with_counterion)
acylium = changed(acylium_with_counterion, remove=(44, 50, 51, 52, 53))

# Make the ortho hydrogen explicit for the rearomatization panel.
rw = Chem.RWMol(acylium)
site = rw.GetAtomWithIdx(index(rw, 32))
site.SetNumExplicitHs(0)
h = Chem.Atom(1)
h.SetAtomMapNum(60)
hi = rw.AddAtom(h)
rw.AddBond(site.GetIdx(), hi, S)
acylium = rw.GetMol()
Chem.SanitizeMol(acylium)
Chem.Kekulize(acylium, clearAromaticFlags=True)
ring_bond = next(b for b in acylium.GetAtomWithIdx(index(acylium, 32)).GetBonds()
                 if b.GetBondType() == DB)
other = ring_bond.GetOtherAtomIdx(index(acylium, 32))
positive = acylium.GetAtomWithIdx(other).GetAtomMapNum()
cyclic_arenium = changed(acylium, [(4, 32, S), (32, positive, S)],
                        [(4, 0), (positive, 1)])
panel("14-ring-closure", "14. Intramolecular Friedel-Crafts acylation",
      "The ortho pi bond of the benzyl ring attacks the tethered acylium "
      "carbon, closing a SIX-membered central ring. AlCl4- is omitted in "
      "this panel for clarity. No molecular coordinates were supplied.",
      [acylium], [Arrow(bond(acylium, 32, positive), atom(acylium, 4))],
      cyclic_arenium)

alcl4 = molecule("[Cl:44][Al-:50]([Cl:51])([Cl:52])[Cl:53]")
anthrone = changed(cyclic_arenium, [(32, positive, DB)],
                   [(positive, 0)], remove=(60,))
panel("15-final-rearomatization", "15. Rearomatize, then aqueous workup",
      "The counterion removes H+ as the C-H pair restores aromaticity. "
      "The Al-Cl pair forms H-Cl, releasing AlCl3. H3O+ workup releases "
      "the ketone from its Lewis-acid complex. The retained bridge is CH2.",
      [cyclic_arenium, alcl4],
      [Arrow(bond(alcl4, 44, 50, 1), atom(cyclic_arenium, 60)),
       Arrow(bond(cyclic_arenium, 32, 60), bond(cyclic_arenium, 32, positive))],
      anthrone)
panel("16-C", "C = ANTHRONE (anthracen-9(10H)-one) -- option (b)",
      "The central ring contains ONE carbonyl and ONE CH2 bridge. "
      "Option (a), anthraquinone, has two carbonyls and ignores the "
      "Clemmensen reduction in the question.", [anthrone])

for mol, formula in ((acid_a, "C14H10O3"), (acid_b, "C14H12O2"),
                     (anthrone, "C14H10O")):
    counts = Counter(atom.GetSymbol() for atom in mol.GetAtoms())
    counts["H"] += sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())
    actual = "".join(symbol + (str(counts[symbol]) if counts[symbol] != 1 else "")
                     for symbol in ("C", "H", "O") if counts[symbol])
    assert actual == formula, (actual, formula)
    assert Chem.GetFormalCharge(mol) == 0
assert sorted(len(r) for r in anthrone.GetRingInfo().AtomRings()) == [6, 6, 6]
assert sum(b.GetBondType() == DB and
           {b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum()} == {6, 8}
           for b in anthrone.GetBonds()) == 1

OUT.mkdir(exist_ok=True)


def render_compact():
    font = str(Path(RDConfig.RDDataDir) / "Fonts" / "Telex-Regular.ttf")
    title_font = ImageFont.truetype(font, 34)
    heading_font = ImageFont.truetype(font, 25)
    label_font = ImageFont.truetype(font, 23)
    note_font = ImageFont.truetype(font, 20)
    sheet = Image.new("RGB", (1600, 1470), "white")
    pen = ImageDraw.Draw(sheet)
    colour = "#202020"

    def text(x, y, value, font=label_font, centered=False):
        if centered:
            x -= pen.textlength(value, font=font) / 2
        bounds = pen.textbbox((x, y), value, font=font)
        assert 0 <= bounds[0] and bounds[2] < sheet.width, value
        assert 0 <= bounds[1] and bounds[3] < sheet.height, value
        pen.text((x, y), value, font=font, fill=colour)

    def reaction(x1, x2, y, above=(), below=()):
        pen.line((x1, y, x2, y), fill=colour, width=3)
        pen.polygon([(x2, y), (x2 - 12, y - 6), (x2 - 12, y + 6)],
                    fill=colour)
        for row, line in enumerate(above):
            text((x1 + x2) / 2, y - 28 * (len(above) - row) - 8,
                 line, note_font, True)
        for row, line in enumerate(below):
            text((x1 + x2) / 2, y + 12 + 26 * row, line, note_font, True)

    rendered = []

    def diagram(slug, mols, arrows, x, y, width=430, height=310):
        for arrow in arrows:
            arrow.lineWidth = 0.035
            arrow.headLength = 0.20
            arrow.headWidth = 0.08
        scene = Scene(slug, "", "", [clean(mol) for mol in mols], list(arrows),
                      width=width, height=height)
        before = [mol.ToBinary() for mol in mols]
        assert all(mol.GetNumConformers() == 0 for mol in mols)
        assert all(arrow.autoRoute for arrow in arrows)
        for backend in ("svg", "png"):
            data, _ = render_scene(scene, backend, font_size=23)
            path = OUT / f"compact-{slug}.{backend}"
            if backend == "svg":
                path.write_text(data, encoding="utf-8")
            else:
                path.write_bytes(data)
                sheet.paste(Image.open(io.BytesIO(data)).convert("RGB"), (x, y))
        assert before == [mol.ToBinary() for mol in mols]
        rendered.append({"id": slug, "arrows": len(arrows)})

    # O-Al speciation is suppressed, leaving the formal carboxylate/acylium pair.
    simple_acyl = changed(open_acyl, charges=[(3, -1)], remove=(20, 21, 22, 23))
    simple_arenium = changed(arenium, charges=[(3, -1)], remove=(20, 21, 22, 23))
    benzene_neighbor = benzene.GetAtomWithIdx(neighbor).GetAtomMapNum()
    first_attack = Arrow(bond(benzene, 31, benzene_neighbor, 1),
                        atom(simple_acyl, 1))
    first_restore = Arrow(bond(simple_arenium, 31, 37),
                          bond(simple_arenium, 31, benzene_neighbor))

    # Either ortho site is equivalent; use the other one for the compact scheme.
    closing = changed(acyl_chloride, charges=[(4, 1)], remove=(44,))
    rw = Chem.RWMol(closing)
    site = rw.GetAtomWithIdx(index(rw, 36))
    site.SetNumExplicitHs(0)
    hydrogen = Chem.Atom(1)
    hydrogen.SetAtomMapNum(60)
    hydrogen_index = rw.AddAtom(hydrogen)
    rw.AddBond(site.GetIdx(), hydrogen_index, S)
    closing = rw.GetMol()
    Chem.SanitizeMol(closing)
    Chem.Kekulize(closing, clearAromaticFlags=True)
    pi_bond = next(b for b in closing.GetAtomWithIdx(index(closing, 36)).GetBonds()
                   if b.GetBondType() == DB)
    partner = closing.GetAtomWithIdx(
        pi_bond.GetOtherAtomIdx(index(closing, 36))).GetAtomMapNum()
    sigma = changed(closing, [(4, 36, S), (36, partner, S)],
                    [(4, 0), (partner, 1)])
    final = changed(sigma, [(36, partner, DB)], [(partner, 0)], remove=(60,))
    final_check, expected = clean(final), clean(anthrone)
    Chem.SanitizeMol(final_check)
    Chem.SanitizeMol(expected)
    assert Chem.MolToSmiles(final_check) == Chem.MolToSmiles(expected)

    text(40, 25, "Question 6 | Condensed mechanism", title_font)
    text(40, 72, "Benzene + phthalic anhydride  ->  A  ->  B  ->  C", note_font)

    text(40, 115, "1. Friedel-Crafts acylation", heading_font)
    text(1015, 119, "AlCl3 activation; then aqueous workup", note_font)
    diagram("acylation", [simple_acyl, benzene], [first_attack], 25, 155)
    diagram("first-sigma", [simple_arenium], [first_restore], 580, 155)
    diagram("A", [acid_a], [], 1140, 155)
    reaction(460, 565, 305)
    reaction(1015, 1125, 305, ("-H+",), ("H3O+",))
    text(240, 477, "Activated anhydride + benzene", note_font, True)
    text(795, 477, "Sigma complex", label_font, True)
    text(1355, 477, "A: 2-benzoylbenzoic acid", note_font, True)

    pen.line((40, 530, 1560, 530), fill="#cccccc", width=1)
    text(40, 550, "2. Reduce the ketone; activate the acid", heading_font)
    text(235, 720, "A", title_font, True)
    reaction(320, 565, 740, ("Zn-Hg, HCl", "heat"), ("C=O -> CH2",))
    diagram("B", [acid_b], [], 580, 590)
    reaction(1015, 1125, 740, ("SOCl2",), ("-SO2", "-HCl"))
    diagram("acid-chloride", [acyl_chloride], [], 1140, 590)
    text(795, 907, "B: 2-benzylbenzoic acid", note_font, True)
    text(1355, 907, "Acid chloride", label_font, True)

    pen.line((40, 955, 1560, 955), fill="#cccccc", width=1)
    text(40, 975, "3. Intramolecular acylation (AlCl3)", heading_font)
    diagram("cyclization", [closing],
            [Arrow(bond(closing, 36, partner), atom(closing, 4))], 25, 1015)
    diagram("cyclic-sigma", [sigma],
            [Arrow(bond(sigma, 36, 60), bond(sigma, 36, partner))], 580, 1015)
    diagram("C", [anthrone], [], 1140, 1015)
    reaction(460, 565, 1165, ("6-membered", "ring closure"))
    reaction(1015, 1125, 1165, ("-H+",), ("H3O+",))
    text(240, 1330, "Acylium ion from the acid chloride", note_font, True)
    text(795, 1330, "Cyclic sigma complex", label_font, True)
    text(1355, 1330, "C: anthrone -- (b)", heading_font, True)
    text(40, 1395, "Aluminium complexes and proton-transfer partners are omitted.",
         note_font)
    text(40, 1423, "Clemmensen reduction is shown as a net step. C retains one C=O and one CH2.",
         note_font)

    sheet.save(OUT / "textbook-mechanism.png")
    sheet.save(OUT / "textbook-mechanism.pdf", resolution=180)
    (OUT / "textbook.html").write_text(
        '<!doctype html><html><head><meta charset="utf-8">'
        '<title>Question 6 - condensed mechanism</title>'
        '<style>body{margin:16px auto;max-width:1250px;font:16px system-ui}'
        'img{width:100%}a{margin-right:20px}@media print{nav{display:none}'
        'body{margin:0}}@page{size:landscape;margin:8mm}</style></head><body>'
        '<nav><a href="textbook-mechanism.pdf">One-page PDF</a>'
        '<a href="textbook-mechanism.png">PNG</a>'
        '</nav>'
        '<img src="textbook-mechanism.png" alt="Condensed question 6 mechanism">'
        '</body></html>', encoding="utf-8")
    (OUT / "compact-manifest.json").write_text(
        json.dumps(rendered, indent=2), encoding="utf-8")
    print(f"One-page mechanism: {OUT / 'textbook.html'}")
    print(f"PNG: {OUT / 'textbook-mechanism.png'}")
    print(f"PDF: {OUT / 'textbook-mechanism.pdf'}")


if args.compact:
    render_compact()
    raise SystemExit(0)

font = str(Path(RDConfig.RDDataDir) / "Fonts" / "Telex-Regular.ttf")
regular = ImageFont.truetype(font, 22)
bold = ImageFont.truetype(font, 26)
small = ImageFont.truetype(font, 19)
cards, sections, manifest = [], [], []
for scene, outcome in panels:
    entry = {"id": scene.slug, "title": scene.title,
             "description": scene.description, "arrows": len(scene.arrows)}
    snapshots = [m.ToBinary() for m in scene.mols]
    pictures = []
    for backend in ("svg", "png"):
        payload, positions = render_scene(scene, backend, font_size=24)
        assert all(len(pos) == mol.GetNumAtoms()
                   for pos, mol in zip(positions, scene.mols))
        path = OUT / f"{scene.slug}.{backend}"
        if backend == "svg":
            path.write_text(payload, encoding="utf-8")
        else:
            path.write_bytes(payload)
            pictures.append(Image.open(io.BytesIO(payload)).convert("RGB"))
    assert snapshots == [m.ToBinary() for m in scene.mols]
    product_html = ""
    if outcome is not None:
        result = Scene(scene.slug + "-result", "Intermediate / result", "",
                       [clean(outcome)], [], width=650, height=330)
        for backend in ("svg", "png"):
            payload, _ = render_scene(result, backend, font_size=24)
            path = OUT / f"{result.slug}.{backend}"
            if backend == "svg":
                path.write_text(payload, encoding="utf-8")
            else:
                path.write_bytes(payload)
                pictures.append(Image.open(io.BytesIO(payload)).convert("RGB"))
        product_html = (f'<div class="result"><b>Intermediate / result</b>'
                        f'<img src="{result.slug}.svg"></div>')
    sections.append(
        f'<section><h2>{html.escape(scene.title)}</h2>'
        f'<p>{html.escape(scene.description)}</p><div class="reaction">'
        f'<img src="{scene.slug}.svg">{product_html}</div></section>')
    card = Image.new("RGB", (1700, 630), "white")
    draw = ImageDraw.Draw(card)
    draw.text((25, 16), scene.title, font=bold, fill="#172b4d")
    for line_no, line in enumerate(textwrap.wrap(scene.description, 145)):
        draw.text((25, 55 + line_no * 25), line, font=small, fill="#444444")
    card.paste(pictures[0], (15, 165))
    if len(pictures) == 2:
        draw.line((1020, 350, 1060, 350), fill="#444444", width=3)
        draw.polygon([(1060, 350), (1048, 344), (1048, 356)], fill="#444444")
        product = pictures[1].resize((610, 310))
        card.paste(product, (1075, 200))
        draw.text((1090, 535), "Intermediate / result", font=small, fill="#555555")
    cards.append(card)
    manifest.append(entry)

overview = Image.new("RGB", (1700, 610), "white")
draw = ImageDraw.Draw(overview)
draw.text((30, 18), "QUESTION 6: A -> B -> C = ANTHRONE, option (b)",
          font=bold, fill="#172b4d")
overview_labels = [
    "A: 2-benzoylbenzoic acid", "B: 2-benzylbenzoic acid", "C: anthrone"]
for pos, (mol, label) in enumerate(zip((acid_a, acid_b, anthrone), overview_labels)):
    scene = Scene("overview", label, "", [clean(mol)], [], width=500, height=320)
    data, _ = render_scene(scene, "png", font_size=25)
    overview.paste(Image.open(io.BytesIO(data)).convert("RGB"), (20 + 575 * pos, 130))
    draw.text((25 + 575 * pos, 480), label, font=regular, fill="#172b4d")
draw.text((480, 85), "Zn-Hg, HCl, heat", font=small, fill="#444444")
draw.text((1050, 85), "1. SOCl2   2. AlCl3   3. H3O+", font=small, fill="#444444")
for x in (530, 1100):
    draw.line((x, 270, x + 40, 270), fill="#444444", width=3)
    draw.polygon([(x + 40, 270), (x + 28, 264), (x + 28, 276)], fill="#444444")
draw.text((30, 550), "Starting materials: benzene + phthalic anhydride; "
          "1. AlCl3, 2. H3O+ gives A.", font=regular, fill="#444444")
overview.save(OUT / "overview.png")

for page in range(4):
    selected = cards[page * 4:(page + 1) * 4]
    sheet = Image.new("RGB", (1700, 630 * len(selected)), "#dddddd")
    for row, card in enumerate(selected):
        sheet.paste(card, (0, row * 630))
    sheet.save(OUT / f"mechanism-page-{page + 1}.png")
full = Image.new("RGB", (1700, 610 + 630 * len(cards)), "white")
full.paste(overview, (0, 0))
for row, card in enumerate(cards):
    full.paste(card, (0, 610 + row * 630))
full.save(OUT / "complete-mechanism.png")

page = """<!doctype html><html><head><meta charset="utf-8">
<title>Question 6 - anthrone mechanism</title><style>
body{font:17px/1.5 system-ui;margin:30px auto;max-width:1700px;color:#172b4d}
h1{margin-bottom:8px}section{border-top:2px solid #dce4ed;padding:18px 0}
.reaction{display:flex;align-items:center;gap:20px}.reaction>img{width:60%}
.result{width:38%;text-align:center}.result img{width:100%}
.note{background:#fff5d9;padding:15px}p{max-width:1250px}
@media(max-width:800px){.reaction{display:block}.reaction>img,.result{width:100%}}
</style></head><body><h1>Question 6: anthrone is option (b)</h1>
<p>Endpoint-only native RDKit mechanism rendering. Molecules and intermediates
are supplied explicitly; RDKit is not predicting chemistry. Every electron
arrow supplies only start/end anchors. No conformers, molecule transforms,
or curve control points are authored.</p>
<p class="note">Scope: detailed organic electron-flow steps, with selected
Lewis-acid adducts and one arenium resonance contributor. Aluminium speciation
and aqueous workup are abbreviated. Clemmensen reduction is an overall
metal-surface transformation, not a claimed definitive elementary mechanism.</p>
<img src="overview.png" style="width:100%">
""" + "\n".join(sections) + "</body></html>"
(OUT / "index.html").write_text(page, encoding="utf-8")
(OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
print(f"Rendered {len(panels)} steps: {OUT / 'index.html'}")
print(f"Overview: {OUT / 'overview.png'}")
print(f"Complete mechanism: {OUT / 'complete-mechanism.png'}")
