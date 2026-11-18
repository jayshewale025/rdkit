// Exercise the compiled module in Chromium, not just Node. No Python/backend server involved.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const http = require('node:http');
const path = require('node:path');
const { chromium } = require('playwright');
const { options, styled } = require('./mechanism');
const cases = require('./mechanism-cases.json').map(styled);
const dist = path.resolve(__dirname, '../build/dist');
const output = path.join(__dirname, 'mechanism-output');

async function main() {
  const server = http.createServer((req, res) => {
    if (req.url === '/RDKit_minimal.js' || req.url === '/RDKit_minimal.wasm') {
      res.setHeader('content-type', req.url.endsWith('.wasm') ? 'application/wasm' : 'text/javascript');
      fs.createReadStream(path.join(dist, req.url.slice(1))).pipe(res);
    } else if (req.url === '/') {
      res.end('<!doctype html><meta charset="utf-8"><script src="/RDKit_minimal.js"></script><main id="gallery"></main>');
    } else {
      res.writeHead(404).end();
    }
  });
  await new Promise((resolve) => server.listen(0, '127.0.0.1', resolve));
  const browser = await chromium.launch();
  try {
    const page = await browser.newPage({ viewport: { width: 900, height: 800 } });
    await page.goto(`http://127.0.0.1:${server.address().port}/`);
    const result = await page.evaluate(async ({ cases, options }) => {
      const start = performance.now();
      const rdkit = await initRDKitModule();
      const loadMs = performance.now() - start;
      const t = performance.now();
      const svgs = cases.map((spec) => rdkit.get_mechanism_svg(JSON.stringify(spec), JSON.stringify(options)));
      const drawMs = performance.now() - t;
      document.getElementById('gallery').innerHTML = svgs.map((svg, i) =>
        `<section><h2>${cases[i].name}</h2>${svg}</section>`).join('');
      return { svgs, loadMs, drawMs };
    }, { cases, options });
    for (let i = 0; i < cases.length; ++i) {
      // Same WASM in Node and a real browser must produce identical SVG.
      assert.equal(result.svgs[i], fs.readFileSync(path.join(output, `${cases[i].name}.svg`), 'utf8'));
    }
    await page.screenshot({ path: path.join(output, 'browser-gallery.png'), fullPage: true });
    await page.setViewportSize({ width: 390, height: 844 });
    await page.addStyleTag({ content: 'svg{max-width:100%;height:auto}' });
    await page.screenshot({ path: path.join(output, 'browser-mobile.png'), fullPage: true });
    fs.writeFileSync(path.join(output, 'browser-metrics.json'), JSON.stringify({
      loadMs: result.loadMs, drawMs: result.drawMs, diagrams: cases.length,
      note: 'GitHub runner desktop Chromium; mobile screenshot tests layout, not physical phone performance.',
    }, null, 2));
    console.log(`PASS browser/Node SVG equality: ${cases.length} diagrams; load ${result.loadMs.toFixed(1)} ms, draw ${result.drawMs.toFixed(1)} ms`);
  } finally {
    await browser.close();
    await new Promise((resolve) => server.close(resolve));
  }
}

main().catch((e) => { console.error(e); process.exitCode = 1; });
