"use strict";

// Examples/embedded/*.js must match the example files they copy
// (regenerate with: python tools/make_embedded_examples.py).
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const zlib = require("node:zlib");

const ROOT = path.join(__dirname, "..");

test("embedded example copies match the example files", () => {
  const dir = path.join(ROOT, "Examples", "embedded");
  const scripts = fs.readdirSync(dir).filter((f) => f.endsWith(".js"));
  const sandbox = { window: {} };
  vm.createContext(sandbox);
  for (const f of scripts) vm.runInContext(fs.readFileSync(path.join(dir, f), "utf8"), sandbox);
  const store = sandbox.window.__3dsExamples;
  const urls = Object.keys(store);
  assert.equal(urls.length, scripts.length);
  for (const url of urls) {
    const copy = zlib.gunzipSync(Buffer.from(store[url], "base64"));
    assert.ok(copy.equals(fs.readFileSync(path.join(ROOT, url))), `${url} differs from its embedded copy`);
  }
  // the page lists exactly these, and each is behind an example button
  const html = fs.readFileSync(path.join(ROOT, "index.html"), "utf8");
  const listed = JSON.parse(/var EMBEDDED_EXAMPLES=(\[[^\]]*\])/.exec(html)[1]);
  assert.deepEqual([...listed].sort(), [...urls].sort());
  const buttons = new Set([...html.matchAll(/data-example="([^"]+)"/g)].flatMap((m) => m[1].split(/\s+/)));
  for (const url of listed) assert.ok(buttons.has(url), `${url} has no example button`);
});
