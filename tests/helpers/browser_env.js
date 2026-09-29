"use strict";

// Loads the app's classic browser scripts (js/*.js) into one isolated global
// scope so Node tests can exercise them without a browser.
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");

const ROOT = path.resolve(__dirname, "..", "..");

function loadScripts(files) {
  const ctx = {
    console, process, require, __dirname: ROOT, __filename: path.join(ROOT, "index.html"),
    TextDecoder, TextEncoder, WebAssembly, performance, URL, setTimeout, clearTimeout,
    ArrayBuffer, Uint8Array, Int8Array, Int32Array, Uint32Array, Float32Array, Float64Array,
    BigInt, BigInt64Array, BigUint64Array, Date, Math, Map, Set,
  };
  ctx.self = ctx;
  ctx.window = ctx;
  ctx.globalThis = ctx;
  vm.createContext(ctx);
  for (const file of files) {
    vm.runInContext(fs.readFileSync(path.join(ROOT, file), "utf8"), ctx, { filename: file });
  }
  return ctx;
}

async function loadH5(extraFiles = []) {
  const ctx = loadScripts(["js/h5wasm.js", ...extraFiles]);
  await ctx.h5wasm.ready;
  return ctx;
}

let fixtureSequence = 0;

// Writes /entry/data with the given datasets ({name: {data, shape, dtype}})
// and attributes, returning the file bytes.
function writeEntryData(h5, datasets, attributes = {}) {
  const file = `/tmp/test_fixture_${Date.now()}_${fixtureSequence++}.h5`;
  const f = new h5.File(file, "w");
  try {
    const data = f.create_group("entry").create_group("data");
    for (const [name, def] of Object.entries(datasets)) data.create_dataset({ name, ...def });
    for (const [name, value] of Object.entries(attributes)) data.create_attribute(name, value);
    f.flush();
  } finally {
    f.close();
  }
  const bytes = Uint8Array.from(h5.FS.readFile(file));
  h5.FS.unlink(file);
  return bytes;
}

module.exports = { ROOT, loadScripts, loadH5, writeEntryData };
