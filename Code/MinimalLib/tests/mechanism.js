// Native mechanism API exercised through the real compiled WebAssembly module.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const initRDKitModule = require('../build/dist/RDKit_minimal.js');
const cases = require('./mechanism-cases.json');

const output = path.join(__dirname, 'mechanism-output');
const options = {
  width: 640, height: 320, fixedBondLength: 34, bondLineWidth: 1.6,
  clearBackground: false, padding: 0.06, fontFile: 'BuiltinTelexRegular',
};

// Match the production Python renderer's readable arrow settings.
function styled(spec) {
  const copy = JSON.parse(JSON.stringify(spec));
  for (const arrow of copy.arrows) {
    Object.assign(arrow, { lineWidth: 0.05, headLength: 0.26, headWidth: 0.1 });
    if (arrow.start.type === 'lonePair') {
      Object.assign(arrow.start, { lonePairRadius: 0.045, lonePairSeparation: 0.17 });
    }
  }
  return copy;
}

async function main() {
  const rdkit = await initRDKitModule();
  assert.equal(typeof rdkit.get_mechanism_svg, 'function');
  fs.mkdirSync(output, { recursive: true });
  for (const spec of cases) {
    const svg = rdkit.get_mechanism_svg(JSON.stringify(styled(spec)), JSON.stringify(options));
    assert.match(svg, /<svg/);
    assert.match(svg, /class='atom-\d+'/); // FreeType glyph paths, not browser-dependent text
    assert.equal((svg.match(/mechanism-shaft/g) || []).length, spec.arrows.length);
    assert.equal((svg.match(/mechanism-head/g) || []).length, spec.arrows.length);
    assert.doesNotMatch(svg, /NaN|Infinity|<script/);
    assert.match(svg, /viewBox='0 0 640 320'/);
    fs.writeFileSync(path.join(output, `${spec.name}.svg`), svg);
    console.log(`PASS ${spec.name} (${svg.length} bytes)`);
  }
  // Default options/one-argument overload, no arrows, and unmodified molecule API.
  assert.match(rdkit.get_mechanism_svg(JSON.stringify(cases[2])), /<svg/);
  assert.match(rdkit.get_mechanism_svg(JSON.stringify({ molecules: ['CCO'], arrows: [] })), /<svg/);
  const mol = rdkit.get_mol('c1ccccc1O');
  assert.ok(mol);
  assert.match(mol.get_svg(), /<svg/);
  mol.delete();

  function bad(spec, detail = '{}', message = /./) {
    assert.throws(() => rdkit.get_mechanism_svg(
      typeof spec === 'string' ? spec : JSON.stringify(spec), detail,
    ), (e) => e instanceof Error && message.test(e.message));
  }
  const base = { molecules: ['CBr'], arrows: [{
    start: { type: 'bond', fragment: 0, index: 0 },
    end: { type: 'atom', fragment: 0, index: 1 },
  }] };
  const changeArrow = (props) => ({ ...base, arrows: [{ ...base.arrows[0], ...props }] });
  bad('{');
  bad([]);
  bad({ molecules: [], arrows: [] }, '{}', /non-empty/);
  bad({ molecules: ['C1CC'], arrows: [] }, '{}', /invalid molecule/);
  bad({ molecules: ['CC'], arrows: {} }, '{}', /arrows/);
  bad(changeArrow({ start: { type: 'banana', fragment: 0, index: 0 } }), '{}', /anchor.type/);
  bad(changeArrow({ start: { type: 'bond', fragment: 8, index: 0 } }), '{}', /fragment.*range/);
  bad(changeArrow({ start: { type: 'atom', fragment: 0, index: 99 } }), '{}', /index.*range/);
  bad(changeArrow({ start: { type: 'atom', fragment: -1, index: 0 } }), '{}', /integer/);
  bad(changeArrow({ start: { type: 'atom', fragment: 0.5, index: 0 } }), '{}', /integer/);
  bad(changeArrow({ head: 'banana' }), '{}', /head/);
  bad(changeArrow({ lineWidth: -1 }));
  bad(changeArrow({ control1: [0, 1] }), '{}', /both control/);
  bad(changeArrow({ colour: [9, 0, 0] }), '{}', /colour/);
  bad({ ...base, offsets: [[0, 0], [1, 1]] });
  bad(base, '{');
  bad(base, '{"width":0}', /width/);
  bad(base, '{"width":100000}', /width/);
  bad(base, '{"height":-1}', /integer/);
  // The runtime remains usable after exceptions, and repeated calls own/free their C++ data.
  const spec = JSON.stringify(styled(cases[2]));
  const details = JSON.stringify(options);
  const expected = rdkit.get_mechanism_svg(spec, details);
  for (let i = 0; i < 200; ++i) {
    assert.equal(rdkit.get_mechanism_svg(spec, details), expected);
  }
  console.log('PASS invalid input, unchanged molecule API and 200 deterministic repeated draws');
}

if (require.main === module) main().catch((e) => { console.error(e); process.exitCode = 1; });
module.exports = { options, styled };
