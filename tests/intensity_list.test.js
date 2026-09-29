"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const { loadScripts } = require("./helpers/browser_env");

const ctx = loadScripts(["js/intensity_list_reader.js"]);
const R = ctx.IntensityListReader;

// Spinteract / Scatty style rows: fixed-width Fortran numbers, a weight column.
function row(h, k, l, v) {
  return `  ${h.toFixed(16)}  ${k.toFixed(16)}  ${l.toFixed(16)}  ${v.toExponential(16)}  1.00000000`;
}

const H = [-0.5, 0, 0.5], K = [-1, -0.5, 0, 0.5], L = [0, 0.04];
const value = (h, k, l) => 100 * h + 10 * k + l + 1000;
function gridText({ skip = () => false, shuffle = false } = {}) {
  const lines = [];
  for (const h of H) for (const k of K) for (const l of L) if (!skip(h, k, l)) lines.push(row(h, k, l, value(h, k, l)));
  if (shuffle) lines.reverse();
  return lines.join("\n") + "\n";
}

test("detects intensity lists and nothing else", () => {
  assert.equal(R.detect(gridText()), true);
  assert.equal(R.detect("  -2.4  -3.0  -3.0  1.604  1.0\n  -2.36  -3.0  -3.0  1.412  1.0\n  -2.32  -3.0  -3.0  1.408  1.0\n"), true);
  assert.equal(R.detect("(Version 6f format configuration file)\nAtoms:\n1 Mn 0 0 0 1 0 0 0\n"), false);
  assert.equal(R.detect("LAMMPS data file\n\n8 atoms\n1 atom types\n"), false);
  assert.equal(R.detect("TITLE MnO\nCELL 4.4 4.4 4.4 90 90 90\nSITE 0 0 0\n"), false);
  assert.equal(R.detect("0.0917 864 1\n864\n0.1 0.2 0.3\n0.1 0.2 0.3\n"), false);
});

test("places points on their hkl grid with l fastest; missing points are NaN", () => {
  const g = R.parseText(gridText({ skip: (h, k, l) => h === 0 && k === 0.5, shuffle: true }));
  assert.deepEqual(Array.from(g.shape), [3, 4, 2]);
  assert.deepEqual(Array.from(g.h), H);
  assert.deepEqual(Array.from(g.k), K);
  assert.deepEqual(Array.from(g.l), L);
  for (let i = 0; i < 3; i++) for (let j = 0; j < 4; j++) for (let m = 0; m < 2; m++) {
    const v = g.I[(i * 4 + j) * 2 + m];
    if (H[i] === 0 && K[j] === 0.5) assert.ok(Number.isNaN(v));
    else assert.ok(Math.abs(v - value(H[i], K[j], L[m])) < 1e-9);
  }
  assert.equal(g.masked, 2);
  assert.equal(g.filled, 22);
  assert.equal(g.min, value(-0.5, -1, 0));
  assert.equal(g.max, value(0.5, 0.5, 0.04));
});

test("chunked input gives the same grid, duplicates are averaged", () => {
  const text = gridText() + row(0.5, 0.5, 0.04, 0) + "\n";
  const whole = R.parseText(text);
  const parser = R.createParser();
  for (let i = 0; i < text.length; i += 7) parser.push(text.slice(i, i + 7));
  const chunked = R.buildGrid(parser.finish());
  assert.deepEqual(Array.from(chunked.I), Array.from(whole.I));
  assert.equal(whole.duplicates, 1);
  assert.ok(Math.abs(whole.I[whole.I.length - 1] - value(0.5, 0.5, 0.04) / 2) < 1e-9);
});

test("a plane gives a singleton axis; gaps of whole planes stay on the grid", () => {
  const plane = R.parseText([row(0, 0, 0, 1), row(0.1, 0, 0, 2), row(0, 0.1, 0, 3), row(0.1, 0.1, 0, 4)].join("\n"));
  assert.deepEqual(Array.from(plane.shape), [2, 2, 1]);
  const gap = R.parseText([row(0, 0, 0, 1), row(0.3, 0, 0, 2), row(0.1, 0, 0, 3)].join("\n"));
  assert.deepEqual(Array.from(gap.h), [0, 0.1, 0.2, 0.3]);
  assert.ok(Number.isNaN(gap.I[2]));
});

test("irregular axes and non-numeric files are refused", () => {
  assert.throws(() => R.parseText([row(0, 0, 0, 1), row(0.1, 0, 0, 2), row(0.25, 0, 0, 3)].join("\n")), /regular grid/);
  assert.throws(() => R.parseText("TITLE x\n"), /no data lines/);
});

test("reads Spinteract configurations and recognises Scatty settings files", () => {
  const spinteract = "TITLE MnO\n\nCELL 4.4344 4.4344 4.4344 90 90 90 \nCENTRING F\nSITE 0.0 0.0 0.0\n\nSPIN_DIMENSION 3\nX_AXIS 6.0 0.0 0.0 151\n";
  assert.equal(R.detectConfig(spinteract), "spinteract");
  assert.deepEqual(JSON.parse(JSON.stringify(R.parseConfig(spinteract))), { title: "MnO", cell: [4.4344, 4.4344, 4.4344, 90, 90, 90] });
  assert.equal(R.detectConfig("NAME hk0\nCENTRE 0 0 0\nX_AXIS 5 0 0 50\nY_AXIS 0 5 0 50\nRADIATION N\n"), "scatty-settings");
  assert.equal(R.detectConfig(" Name: vol\n Window type: lanczos\n"), "scatty-info");
  assert.equal(R.detectConfig("TITLE mno\nCELL 4.4 4.4 4.4 90 90 90\nSITE 0 0 0\nBOX 6 6 6\nSPIN 1 0 0 0 1 0 0\n"), null);
});
