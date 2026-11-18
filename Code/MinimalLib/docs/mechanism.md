# Reaction mechanisms in RDKit.js (fork prototype)

This fork adds `RDKit.get_mechanism_svg(spec_json, details_json?)`. It calls the
**same C++ `MolDraw2D::drawMechanism` implementation as Python**. SVG is generated
inside WebAssembly in the browser; it needs no Python server, screenshot service,
or network request after the module has loaded.

This API draws the electron flow **you specify**. It does not predict mechanisms,
infer products, or prove that a mechanism is chemically correct. Automatic 2D
placement is schematic; it does not guarantee SN2 backside or E2 anti-periplanar
geometry. Use explicitly placed fragments/curves when that geometry matters.

## Get a prototype build

In this fork's **Actions → Mechanism WebAssembly**, download `rdkit-mechanism-wasm`
after a successful run. It contains:

* `dist/RDKit_minimal.js` and `dist/RDKit_minimal.wasm` — serve these together.
* `dist/RDKit_minimal.d.ts` — generated TypeScript declarations (including the API).
* A scoped npm tarball named `@jayshewale025/rdkit-mechanism` (not published to npm).
* `SHA256SUMS` for the binaries and tarball.

These are fork artifacts, **not** the upstream `@rdkit/rdkit` package. Existing
`get_mol` and `get_rxn` APIs remain available. No npm credentials are needed to
build or test; CI never publishes into the upstream npm namespace.

The build pins Emscripten 6.0.13 and includes FreeType with RDKit's built-in fonts.
For Python/browser comparisons, pin the same fork revision, drawing settings and
`fontFile: "BuiltinTelexRegular"`. Floating-point and font-library differences can
still cause small differences; there is no promise of byte-identical Python SVG.

## Example: SN2

```html
<script src="/rdkit/RDKit_minimal.js"></script>
<div id="mechanism"></div>
<script>
  initRDKitModule({
    locateFile: (name) => `/rdkit/${name}`,
  }).then((RDKit) => {
    const specification = {
      molecules: ["[OH-]", "CBr"],
      arrows: [
        {
          start: { type: "lonePair", fragment: 0, index: 0 },
          end: { type: "atom", fragment: 1, index: 0 },
        },
        {
          start: { type: "bond", fragment: 1, index: 0 },
          end: { type: "atom", fragment: 1, index: 1 },
        },
      ],
    };
    const options = {
      width: 640,
      height: 320,
      fixedBondLength: 34,
      bondLineWidth: 1.6,
      clearBackground: false,
      fontFile: "BuiltinTelexRegular",
    };
    try {
      const svg = RDKit.get_mechanism_svg(
        JSON.stringify(specification), JSON.stringify(options),
      );
      document.getElementById("mechanism").innerHTML = svg;
    } catch (error) {
      console.error(error.message); // ordinary JavaScript Error
    }
  });
</script>
```

There are **no returned C++ handles to delete**. Each call owns and frees its
molecules, arrows and drawer; its result is a plain SVG string. Cache that string
in the application to avoid drawing again on every UI update.

## JSON specification

* `molecules`: non-empty array of SMILES strings. Molecules are parsed separately
  in this order; fragment numbers and atom/bond indices are **zero-based**.
  Explicit hydrogen atoms are retained (`removeHs=false`), important for E2.
* `arrows`: array of `{start, end, ...options}`. An empty array draws the fragments.
* `kekulize`: default `false`; `true` converts aromatic flags to explicit alternating
  bonds before drawing. The application must reference the intended Kekulé bond.
* `showAtomMaps`: default `true`; `false` hides map labels without changing indices.
* `autoLayout`: default `true`; has the same conditions/limitations as the Python API.
* `offsets`: optional array of `[x, y]`, one per fragment, in molecular scene units.
* `rotations`: optional array of counterclockwise angles in degrees, one per fragment.

This is deliberately a low-level RDKit API, independent of any tutoring app's
schema. If an app stores atom-map references such as `lp:7` or `bond:2-3`, resolve
them to the molecule's atom/bond indices before calling it. Do not treat atom-map
numbers as RDKit indices.

### Anchors

| Type | Required fields | Meaning |
|---|---|---|
| `atom` | `fragment`, `index` | Atom index in that fragment |
| `bond` | `fragment`, `index` | Bond index in that fragment |
| `lonePair` | `fragment`, `index` | Lone pair on that atom; position selected automatically |
| `point` | `position: [x,y]` | Explicit position in scene units |

Atom/bond/lone-pair anchors accept `offset: [dx,dy]` in the fragment's unrotated
frame. An explicit offset for a lone pair selects manual placement. Additional
fields mirror the native API: `fraction` (position along a bond), `clipToLabel`,
`autoPosition`, `lonePairAngle`, `lonePairRadius`, `lonePairSeparation`.

### Arrows

`head` is `"full"` (electron pair, default), `"fishhookLeft"` or `"fishhookRight"`
(one electron). Optional fields:

* `lineWidth`, `headLength`, `headWidth`: scene units, same as the Python API.
  For screen reading, try `0.05`, `0.26`, `0.1` respectively.
* `colour: [r,g,b]`: components 0–1; overrides automatic colouring.
* `autoColour`, `autoRoute`: native switches, both default `true`.
* `control1: [x,y]`, `control2: [x,y]`: supplying both selects a manual cubic curve
  unless `autoRoute` is explicitly `true`.

### Drawing options

The second JSON string accepts standard RDKit `MolDrawOptions`, as used by
`get_svg_with_highlights`: bond length/width, fonts, palette, padding, background,
etc. `width` and `height` default to 640×320 and must be integers from 1 to 4096.
`noFreetype: true` opts out of outline fonts; keep it `false` for consistent
cross-browser glyphs. This API returns the full canvas; crop the SVG in the app
if desired.

Bad JSON, malformed SMILES, invalid indices, invalid arrow settings or invalid
canvas dimensions throw a catchable JS `Error`. A bad call does not poison the
module; subsequent valid calls continue working.

## Tests and visual review

```bash
node Code/MinimalLib/tests/tests.js              # unchanged upstream API suite
node Code/MinimalLib/tests/mechanism.js          # 8 scenes, invalid inputs, repeated draws
NODE_PATH=/path/to/node_modules \
  node Code/MinimalLib/tests/mechanism-browser.js # actual Chromium; needs Playwright
```

The `mechanism-wasm-examples` artifact includes SVGs, desktop and phone-width
screenshots, and load/draw timings. Phone-width screenshots are **not a benchmark
of a physical low-end phone**. Test a real target device before migrating an app.

Optional Python parity review (requires this fork's Python wheel):

```bash
python Code/MinimalLib/tests/mechanism-reference.py
node Code/MinimalLib/tests/mechanism-compare.js
```

The reference uses the same scene fixtures and options and reports numeric SVG
path differences. CI compares against the Python `mechanism-v1` wheel (whose
native drawing code is unchanged in this PR). If the drawing engine changes
later, update the baseline deliberately. It never makes LLM calls.
