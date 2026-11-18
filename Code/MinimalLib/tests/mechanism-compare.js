// Compare native Python vs WASM path geometry. Different FreeType builds may
// introduce small differences; fail rather than claim parity if tolerance is exceeded.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const cases = require('./mechanism-cases.json');
const output = path.join(__dirname, 'mechanism-output');
const tolerance = Number(process.env.MECHANISM_PARITY_TOLERANCE || 1);

function paths(svg) {
  return [...svg.matchAll(/<path\b[^>]*\bd='([^']*)'[^>]*>/gs)].map((m) => ({
    commands: (m[1].match(/[A-Za-z]/g) || []).join(''),
    coords: (m[1].match(/-?\d+(?:\.\d+)?/g) || []).map(Number),
  }));
}

const report = [];
for (const spec of cases) {
  const native = paths(fs.readFileSync(path.join(output, `python-${spec.name}.svg`), 'utf8'));
  const wasm = paths(fs.readFileSync(path.join(output, `${spec.name}.svg`), 'utf8'));
  assert.equal(wasm.length, native.length, `${spec.name}: number of SVG paths`);
  let max = 0;
  for (let i = 0; i < native.length; ++i) {
    assert.equal(wasm[i].commands, native[i].commands, `${spec.name}: path ${i} commands`);
    assert.equal(wasm[i].coords.length, native[i].coords.length, `${spec.name}: path ${i} coordinates`);
    wasm[i].coords.forEach((v, j) => { max = Math.max(max, Math.abs(v - native[i].coords[j])); });
  }
  report.push({ name: spec.name, paths: native.length, maxDifferencePx: max });
  console.log(`${spec.name}: ${native.length} paths, largest coordinate difference ${max.toFixed(3)} px`);
  assert.ok(max <= tolerance, `${spec.name}: ${max} exceeds ${tolerance}px tolerance; inspect both SVGs`);
}
fs.writeFileSync(path.join(output, 'python-parity.json'), JSON.stringify(report, null, 2));
console.log(`PASS Python/WASM path geometry within ${tolerance}px for all ${cases.length} examples`);
