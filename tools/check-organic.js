#!/usr/bin/env node
/* Checks the organic (Functional Group Flashcards) data in index.html.
   Run: node tools/check-organic.js
   - group drawings: bond angles match the house style (109.5° at sp3 atoms, 120° at sp2, 180° at sp)
   - every group has a describing label (alt) that never names a group
   - every SLIPS pair names real groups; a phenol item never offers "arene" as a wrong answer
   - Spot It molecules: bond lengths 1, no atoms too close, sane valence, groups exist, atoms in range */
const fs = require('fs'), path = require('path');
const html = fs.readFileSync(path.join(__dirname, '..', 'index.html'), 'utf8');
const js = html.match(/<script>([\s\S]*?)<\/script>/g).map(x => x.slice(8, -9)).join('\n');
const grab = (start, end) => { const a = js.indexOf(start); if (a < 0) throw new Error('missing ' + start); const b = js.indexOf(end, a); return js.slice(a, b + end.length); };
const code = [grab("const R = 'R'", '\n];'), grab('const ALT =', 'const N = GROUPS.length;'), grab('const SLIPS = [', '\n];'),
  grab('const RA = ', '\n];'), 'return {GROUPS, SLIPS, SPOT_MOLS};'].join('\n');
const {GROUPS, SLIPS, SPOT_MOLS} = new Function(code)();
const byId = Object.fromEntries(GROUPS.map(g => [g.id, g]));
let fails = 0; const bad = m => { fails++; console.log('FAIL ' + m); };
const deg = (a, b, c) => { const v1 = [a[0] - b[0], a[1] - b[1]], v2 = [c[0] - b[0], c[1] - b[1]]; return Math.acos((v1[0] * v2[0] + v1[1] * v2[1]) / Math.hypot(...v1) / Math.hypot(...v2)) * 180 / Math.PI; };

/* 1. group drawings */
for (const g of GROUPS) {
  if (!g.alt) bad(g.id + ': no alt label');
  else for (const o of GROUPS) if (new RegExp('\\b' + o.name.toLowerCase().replace(/s$/, '') + '\\b').test(g.alt.toLowerCase()) && !/^(alkane|alkene|alkyne)$/.test(o.id)) bad(`${g.id}: alt names a group (${o.name})`);
  const d = g.draw; if (d.benzene) continue;
  const P = d.atoms.map(a => [a[1], a[2]]);
  d.atoms.forEach((a, i) => {
    const nb = d.bonds.filter(b => b[0] === i || b[1] === i); if (nb.length < 2) return;
    const maxo = Math.max(...nb.map(b => b[2])), want = maxo === 3 ? 180 : maxo === 2 ? 120 : 109.5;
    for (let x = 0; x < nb.length; x++) for (let y = x + 1; y < nb.length; y++) {
      const j = nb[x][0] === i ? nb[x][1] : nb[x][0], k = nb[y][0] === i ? nb[y][1] : nb[y][0], an = deg(P[j], P[i], P[k]);
      if (nb.length === 2 && Math.abs(an - want) > 1.5) bad(`${g.id}: angle at atom ${i} (${a[0] || 'C'}) is ${an.toFixed(1)}°, want ${want}°`);
      if (nb.length === 3 && Math.abs(an - 120) > 1.5) bad(`${g.id}: angle at atom ${i} is ${an.toFixed(1)}°, want 120°`);
    }
  });
  for (const [i, j] of d.bonds) { const L = Math.hypot(P[i][0] - P[j][0], P[i][1] - P[j][1]); if (L < 0.95 || L > 1.25) bad(`${g.id}: bond ${i}-${j} length ${L.toFixed(2)}`); }
}
/* 2. slips and distractors */
for (const x of SLIPS) for (const id of x.pair) if (!byId[id]) bad(`slip ${x.id}: no group ${id}`);
if (!byId.phenol) bad('no phenol group');
if (!/if \(g\.id === 'phenol'\) return pickDistractorsBase\(g, n \+ 1\)\.filter\(x => x\.id !== 'arene'\)/.test(js)) bad('phenol distractor filter missing');

/* 3. Spot It molecules */
const VAL = {C: 4, O: 2, N: 3, S: 2, H: 1};
const el = t => t === '' ? 'C' : t.replace(/^H(?=[A-Z])/, '')[0];
const hOf = t => { if (el(t) === 'H') return 0; const m = t.match(/H(\d?)/); return m ? +(m[1] || 1) : 0; };
for (const m of SPOT_MOLS) {
  const P = m.atoms.map(a => [a[1], a[2]]), n = P.length;
  for (const [i, j] of m.bonds) { if (i >= n || j >= n) { bad(`${m.id}: bond ${i}-${j} out of range`); continue; } const L = Math.hypot(P[i][0] - P[j][0], P[i][1] - P[j][1]); if (Math.abs(L - 1) > 0.02) bad(`${m.id}: bond ${i}-${j} length ${L.toFixed(3)}`); }
  const bonded = new Set(m.bonds.map(([i, j]) => i + '-' + j).concat(m.bonds.map(([i, j]) => j + '-' + i)));
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { const L = Math.hypot(P[i][0] - P[j][0], P[i][1] - P[j][1]); if (!bonded.has(i + '-' + j) && L < 1.4) bad(`${m.id}: atoms ${i} and ${j} only ${L.toFixed(2)} apart`); }
  m.atoms.forEach((a, i) => {
    const sum = m.bonds.filter(b => b[0] === i || b[1] === i).reduce((t, b) => t + b[2], 0) + hOf(a[0]), v = VAL[el(a[0])];
    if (el(a[0]) === 'C' ? sum > 4 : sum !== v) bad(`${m.id}: atom ${i} (${a[0] || 'C'}) has ${sum} bonds, valence ${v}`);
    if (a[0] !== '' && !a[3] && !/^H$/.test(a[0])) bad(`${m.id}: atom ${i} (${a[0]}) has no colour class`);
  });
  for (const [gid, sets] of Object.entries(m.groups)) { if (!byId[gid]) bad(`${m.id}: no group ${gid}`); for (const set of sets) for (const i of set) if (i >= n) bad(`${m.id}: ${gid} atom ${i} out of range`); }
  if (m.groups.phenol) for (const [o, c] of m.groups.phenol) if (!bonded.has(o + '-' + c) || c > 5) bad(`${m.id}: phenol O ${o} is not on ring carbon ${c}`);
}
const tasks = SPOT_MOLS.reduce((t, m) => t + Object.keys(m.groups).length, 0);
console.log(`${GROUPS.length} groups, ${SLIPS.length} slips, ${SPOT_MOLS.length} Spot It molecules, ${tasks} Spot It tasks`);
console.log(fails ? `${fails} problem(s)` : 'all checks pass'); process.exit(fails ? 1 : 0);
