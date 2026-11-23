/* QE Studio front end. Plain JS + Plotly. */
"use strict";
const $ = (s, el = document) => el.querySelector(s);
const $$ = (s, el = document) => [...el.querySelectorAll(s)];
const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const SERIES = () => [1, 2, 3, 4, 5, 6, 7, 8].map((i) => css(`--s${i}`));
const esc = (s) => String(s ?? "").replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
const fmt = (v, d = 4) => (v === null || v === undefined || Number.isNaN(v) ? "–" : Number(v).toFixed(d));
const RY = 13.605693122994;

function toast(msg, ms = 3500) {
  const t = $("#toast");
  t.textContent = msg; t.hidden = false;
  clearTimeout(t._h); t._h = setTimeout(() => (t.hidden = true), ms);
}
async function api(url, opts = {}) {
  const r = await fetch(url, opts);
  const ct = r.headers.get("content-type") || "";
  const body = ct.includes("json") ? await r.json() : await r.text();
  if (!r.ok) throw new Error((body && body.error) || body || r.statusText);
  return body;
}

/* ------------------------------------------------------------------ plot helpers */
const CONFIG = { responsive: true, displaylogo: false, showSendToCloud: false, showEditInChartStudio: false,
  toImageButtonOptions: { format: "png", scale: 3 },
  modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d", "sendDataToCloud"] };
function axis(extra = {}) {
  return { gridcolor: css("--grid"), linecolor: css("--axis"), zeroline: false, tickfont: { color: css("--muted"), size: 11 },
    title: { font: { color: css("--text-secondary"), size: 12 } }, ...extra, ...(extra.title ? { title: { text: extra.title, font: { color: css("--text-secondary"), size: 12 } } } : {}) };
}
function layout(extra = {}) {
  return { paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
    font: { family: 'system-ui, -apple-system, "Segoe UI", sans-serif', color: css("--text-secondary"), size: 12 },
    margin: { l: 60, r: 16, t: 36, b: 48 }, hovermode: "closest",
    hoverlabel: { bgcolor: css("--card"), bordercolor: css("--border"), font: { color: css("--text-primary") } },
    legend: { orientation: "h", y: 1.0, yanchor: "bottom", x: 0, font: { color: css("--text-secondary") } }, ...extra };
}
function alpha(hex, a) {
  const n = parseInt(hex.slice(1), 16);
  return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${a})`;
}

/* Band structure (+ optional DOS/PDOS panel sharing the energy axis) */
function plotBands(el, b, ref, refLabel, dos, pdos, emin, emax, showPdos) {
  const col = SERIES();
  const traces = [];
  const nspin = b.eig.length;
  b.eig.forEach((bandsArr, s) => {
    bandsArr.forEach((band, i) => {
      traces.push({ x: b.x, y: band.map((e) => e - ref), mode: "lines", type: "scatter",
        line: { color: col[s], width: 2, dash: s === 1 ? "dot" : "solid" },
        name: nspin === 2 ? (s === 0 ? "spin up" : "spin down") : "bands", legendgroup: `s${s}`,
        showlegend: nspin === 2 && i === 0,
        hovertemplate: `band ${i + 1}${nspin === 2 ? (s ? " ↓" : " ↑") : ""}<br>E − ${refLabel} = %{y:.3f} eV<extra></extra>` });
    });
  });
  const ann = [];
  const e = b.edges || {};
  if (e.metal === false && e.vbm_index !== undefined) {
    const xv = b.x[e.vbm_index], xc = b.x[e.cbm_index];
    traces.push({ x: [xv, xc], y: [e.vbm_ev - ref, e.cbm_ev - ref], mode: "markers", type: "scatter", showlegend: false,
      marker: { size: 10, color: css("--text-primary"), line: { color: css("--card"), width: 2 } },
      text: ["VBM", "CBM"], hovertemplate: "%{text}: %{y:.3f} eV<extra></extra>" });
    ann.push({ x: xv, y: e.vbm_ev - ref, text: "VBM", showarrow: false, yshift: -14, font: { size: 11, color: css("--text-primary") } });
    ann.push({ x: xc, y: e.cbm_ev - ref, text: `CBM · E<sub>g</sub> = ${e.gap_ev.toFixed(3)} eV (${e.gap_type})`, showarrow: false, yshift: 14,
      xanchor: "left", font: { size: 11, color: css("--text-primary") } });
  }
  const hasDos = dos && dos.energy && dos.energy.length;
  const xmax = b.x[b.x.length - 1];
  const lay = layout({
    xaxis: axis({ range: [0, xmax], tickvals: b.ticks.map((t) => t.x), ticktext: b.ticks.map((t) => t.label),
      gridcolor: css("--axis"), domain: [0, hasDos ? 0.76 : 1], tickfont: { size: 13, color: css("--text-primary") } }),
    yaxis: axis({ range: [emin, emax], title: `E − ${refLabel} (eV)` }),
    shapes: [{ type: "line", xref: "paper", x0: 0, x1: 1, y0: 0, y1: 0, line: { color: css("--muted"), width: 1, dash: "dash" } }],
    annotations: ann, showlegend: true,
  });
  if (hasDos) {
    const E = dos.energy.map((v) => v - ref);
    const tot = dos.dos || dos.dos_up.map((u, i) => u + dos.dos_dw[i]);
    const inWin = E.map((v, i) => (v >= emin && v <= emax ? tot[i] : 0));
    const maxd = Math.max(...inWin, 1e-6);
    traces.push({ x: tot, y: E, xaxis: "x2", yaxis: "y", type: "scatter", mode: "lines", name: "total DOS",
      line: { color: css("--muted"), width: 1.5 }, fill: "tozerox", fillcolor: alpha(css("--muted").startsWith("#") ? css("--muted") : "#898781", 0.18),
      hovertemplate: "DOS %{x:.3f} states/eV<br>E − " + refLabel + " = %{y:.2f} eV<extra></extra>" });
    if (showPdos && pdos && pdos.projections) {
      const Ep = pdos.energy.map((v) => v - ref);
      Object.entries(pdos.projections).slice(0, 7).forEach(([k, v], i) => {
        const y = v.total || v.up.map((u, j) => u + v.dw[j]);
        traces.push({ x: y, y: Ep, xaxis: "x2", yaxis: "y", type: "scatter", mode: "lines", name: k,
          line: { color: col[(i + 1) % 8], width: 2 }, hovertemplate: `${k}: %{x:.3f}<extra></extra>` });
      });
    }
    lay.xaxis2 = axis({ domain: [0.79, 1], range: [0, maxd * 1.08], title: "DOS (states/eV)", showgrid: false });
    lay.showlegend = true;
  }
  Plotly.react(el, traces, lay, CONFIG);
}

function plotPdos(el, pdos, dos, ref, refLabel, emin, emax) {
  const col = SERIES();
  const traces = [];
  if (dos && dos.energy) {
    const E = dos.energy.map((v) => v - ref);
    if (dos.dos) traces.push({ x: E, y: dos.dos, name: "total", mode: "lines", line: { color: css("--muted"), width: 1.5 }, fill: "tozeroy", fillcolor: "rgba(137,135,129,0.15)" });
    else {
      traces.push({ x: E, y: dos.dos_up, name: "total ↑", mode: "lines", line: { color: css("--muted"), width: 1.5 }, fill: "tozeroy", fillcolor: "rgba(137,135,129,0.15)" });
      traces.push({ x: E, y: dos.dos_dw.map((v) => -v), name: "total ↓", mode: "lines", line: { color: css("--muted"), width: 1.5 }, fill: "tozeroy", fillcolor: "rgba(137,135,129,0.15)", showlegend: false });
    }
  }
  if (pdos && pdos.projections) {
    const E = pdos.energy.map((v) => v - ref);
    Object.entries(pdos.projections).slice(0, 7).forEach(([k, v], i) => {
      const c = col[(i + 1) % 8];
      if (v.total) traces.push({ x: E, y: v.total, name: k, mode: "lines", line: { color: c, width: 2 } });
      else {
        traces.push({ x: E, y: v.up, name: k, legendgroup: k, mode: "lines", line: { color: c, width: 2 } });
        traces.push({ x: E, y: v.dw.map((x) => -x), name: k + " ↓", legendgroup: k, showlegend: false, mode: "lines", line: { color: c, width: 2 } });
      }
    });
  }
  traces.forEach((t) => (t.hovertemplate = `${t.name}: %{y:.3f} states/eV<br>E − ${refLabel} = %{x:.2f} eV<extra></extra>`));
  Plotly.react(el, traces, layout({
    xaxis: axis({ range: [emin, emax], title: `E − ${refLabel} (eV)` }), yaxis: axis({ title: "states / eV" }),
    shapes: [{ type: "line", yref: "paper", y0: 0, y1: 1, x0: 0, x1: 0, line: { color: css("--muted"), width: 1, dash: "dash" } }],
    hovermode: "x unified",
  }), CONFIG);
}

function plotScf(el, series) {
  // series: [{name, cycles:[[{iter,energy,accuracy}]]}]
  const col = SERIES();
  const traces = [];
  series.forEach((s, i) => {
    const pts = [];
    let n = 0;
    s.cycles.forEach((c) => c.forEach((p) => { n += 1; if (p.accuracy) pts.push({ n, ...p }); }));
    traces.push({ x: pts.map((p) => p.n), y: pts.map((p) => p.accuracy), name: s.name, mode: "lines+markers",
      line: { color: i < 8 ? col[i] : css("--muted"), width: 2 }, marker: { size: 8, line: { color: css("--card"), width: 2 } },
      customdata: pts.map((p) => p.energy),
      hovertemplate: `${esc(s.name)}<br>iteration %{x}<br>accuracy %{y:.2e} Ry<br>E = %{customdata:.8f} Ry<extra></extra>`,
      showlegend: series.length > 1 && i < 8 });
  });
  Plotly.react(el, traces, layout({ xaxis: axis({ title: "SCF iteration (cumulative over ionic steps)" }),
    yaxis: axis({ title: "estimated scf accuracy (Ry)", type: "log", exponentformat: "power" }) }), CONFIG);
}

function plotConv(el, conv, xTitle, xKey = "value") {
  // the reference (largest value) has ΔE = 0 by definition; it cannot sit on a log axis, so it is listed in the table only
  const rows = (conv.rows || []).filter((r) => r.delta_mev_atom !== undefined && Math.abs(r.delta_mev_atom) > 1e-9);
  const thr = conv.threshold_mev_atom || 1;
  const xs = rows.map((r) => r[xKey]);
  const traces = [{ x: xs, y: rows.map((r) => Math.abs(r.delta_mev_atom)), mode: "lines+markers", name: "|ΔE|",
    line: { color: css("--s1"), width: 2 }, marker: { size: 9, line: { color: css("--card"), width: 2 } },
    text: rows.map((r) => r.label), customdata: rows.map((r) => r.delta_mev_atom),
    hovertemplate: "%{text}<br>ΔE = %{customdata:.3f} meV/atom<extra></extra>", showlegend: false }];
  if (conv.recommended !== null && conv.recommended !== undefined) {
    const r = rows.find((q) => q.value === conv.recommended);
    if (r) traces.push({ x: [r[xKey]], y: [Math.abs(r.delta_mev_atom)], mode: "markers", showlegend: false,
      marker: { size: 14, color: "rgba(0,0,0,0)", line: { color: css("--good"), width: 3 } }, hovertemplate: "recommended: " + esc(r.label) + "<extra></extra>" });
  }
  Plotly.react(el, traces, layout({
    xaxis: axis({ title: xTitle, tickvals: xs, ticktext: rows.map((r) => (xKey === "value" ? r.label.split("= ")[1] : r[xKey])) }),
    yaxis: axis({ title: "|E − E_ref| (meV/atom)", type: "log", exponentformat: "power" }),
    shapes: [{ type: "rect", xref: "paper", x0: 0, x1: 1, y0: 1e-4, y1: thr, fillcolor: alpha(css("--good").startsWith("#") ? css("--good") : "#0ca30c", 0.08), line: { width: 0 } },
      { type: "line", xref: "paper", x0: 0, x1: 1, y0: thr, y1: thr, line: { color: css("--good"), width: 1, dash: "dash" } }],
    annotations: [{ xref: "paper", x: 1, y: Math.log10(thr), text: `${thr} meV/atom`, showarrow: false, xanchor: "right", yshift: 9, font: { size: 11, color: css("--good-text") } }],
  }), CONFIG);
}

function convTable(conv) {
  const rows = conv.rows || [];
  return `<table class="tbl"><thead><tr><th>Setting</th><th class="num">k-points</th><th class="num">E (Ry)</th><th class="num">ΔE (meV/atom)</th><th class="num">P (kbar)</th><th class="num">time (s)</th></tr></thead><tbody>` +
    rows.map((r) => `<tr><td>${esc(r.label)}${r.value === conv.recommended ? " ✓" : ""}${r.ok ? "" : " ⚠"}</td><td class="num">${r.nks ?? "–"}</td><td class="num">${fmt(r.energy_ry, 6)}</td><td class="num">${fmt(r.delta_mev_atom, 3)}</td><td class="num">${fmt(r.pressure_kbar, 2)}</td><td class="num">${fmt(r.wall_s, 1)}</td></tr>`).join("") +
    `</tbody></table><div class="muted small" style="margin-top:6px">${conv.recommended_label ? `Converged to ${conv.threshold_mev_atom} meV/atom from <b>${esc(conv.recommended_label)}</b> (relative to the largest value tested).` : "Not converged within the tested range at this threshold. Test larger values."}</div>`;
}

function plotStructure(el, st) {
  if (!st || !st.positions) { el.innerHTML = '<div class="empty">No structure</div>'; return; }
  const col = SERIES();
  const elems = [...new Set(st.symbols)];
  const c = st.cell;
  const corners = [];
  for (const i of [0, 1]) for (const j of [0, 1]) for (const k of [0, 1])
    corners.push([0, 1, 2].map((d) => i * c[0][d] + j * c[1][d] + k * c[2][d]));
  const edges = [[0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [2, 3], [2, 6], [3, 7], [4, 5], [4, 6], [5, 7], [6, 7]];
  const ex = [], ey = [], ez = [];
  edges.forEach(([a, b]) => { ex.push(corners[a][0], corners[b][0], null); ey.push(corners[a][1], corners[b][1], null); ez.push(corners[a][2], corners[b][2], null); });
  const traces = [{ type: "scatter3d", mode: "lines", x: ex, y: ey, z: ez, line: { color: css("--axis"), width: 3 }, hoverinfo: "skip", showlegend: false }];
  elems.forEach((e, i) => {
    const idx = st.symbols.map((s, j) => (s === e ? j : -1)).filter((j) => j >= 0);
    traces.push({ type: "scatter3d", mode: "markers", name: e, x: idx.map((j) => st.positions[j][0]), y: idx.map((j) => st.positions[j][1]), z: idx.map((j) => st.positions[j][2]),
      marker: { size: 9, color: col[i % 8], line: { color: css("--card"), width: 2 } }, text: idx.map((j) => `${e} #${j + 1}`),
      hovertemplate: "%{text}<br>(%{x:.3f}, %{y:.3f}, %{z:.3f}) Å<extra></extra>" });
  });
  const ax = { showbackground: false, gridcolor: css("--grid"), zeroline: false, color: css("--muted"), title: { text: "" } };
  Plotly.react(el, traces, layout({ margin: { l: 0, r: 0, t: 10, b: 0 }, scene: { xaxis: ax, yaxis: ax, zaxis: ax, aspectmode: "data" },
    showlegend: elems.length > 1 }), CONFIG);
}

/* ------------------------------------------------------------------ navigation */
let current = { view: "new", job: null, poll: null };
function show(view) {
  $$(".view").forEach((v) => v.classList.toggle("active", v.id === `view-${view}`));
  $$(".tab").forEach((t) => t.classList.toggle("active", t.dataset.view === view || (view === "job" && t.dataset.view === "runs")));
  current.view = view;
  if (view === "runs") loadRuns();
  window.dispatchEvent(new Event("resize"));
}
$$(".tab").forEach((t) => t.addEventListener("click", () => { location.hash = t.dataset.view; }));
window.addEventListener("hashchange", route);
function route() {
  const h = location.hash.slice(1) || "new";
  if (h.startsWith("job/")) openJob(h.slice(4));
  else { stopPoll(); show(h); }
}

/* ------------------------------------------------------------------ status */
let SYS = {};
async function loadStatus() {
  try {
    SYS = await api("/api/status");
    const ok = !!SYS.qe_version;
    $("#sys").innerHTML = `<span class="dot" style="background:${ok ? css("--good") : css("--critical")}"></span>` +
      (ok ? `QE ${esc(SYS.qe_version)} · ${SYS.physical_cores} cores` : "pw.x not found – run setup_wsl.sh") +
      (SYS.running ? ` · running 1${SYS.queued ? ` (+${SYS.queued} queued)` : ""}` : "");
    if (!$("#o-np").dataset.touched) $("#o-np").value = SYS.physical_cores;
    $("#np-hint").textContent = `${SYS.physical_cores} physical cores`;
    $("#runs-dir").textContent = SYS.public ? "Runs started from this browser" : SYS.runs_dir;
    applyPublicMode();
  } catch (e) { $("#sys").textContent = "server not reachable"; }
}
$("#o-np").addEventListener("input", (e) => (e.target.dataset.touched = 1));

function applyPublicMode() {
  if (!SYS.public) return;
  document.body.classList.add("public");
  const L = SYS.limits || {};
  const b = $("#public-note");
  if (b) {
    b.hidden = false;
    b.innerHTML = `Shared server: up to ${L.atoms} atoms, ecutwfc ≤ ${L.ecutwfc} Ry, ≤ ${L.kpoints} k-points per grid, ` +
      `${Math.round((L.step_seconds || 0) / 60)} min per step. Runs are deleted after ${L.keep_days} days. Don't upload confidential structures.`;
  }
}
const MY_RUNS_KEY = "qe-studio-my-runs";
const myRuns = () => JSON.parse(localStorage.getItem(MY_RUNS_KEY) || "[]");
function rememberRun(id) { localStorage.setItem(MY_RUNS_KEY, JSON.stringify([...new Set([id, ...myRuns()])].slice(0, 100))); }

/* ------------------------------------------------------------------ new calculation */
let INPUT = { name: "scf.in", pseudos: [] };
let inspectTimer = null;
const editor = $("#editor");

function readFiles(list) {
  const files = [...list];
  const upf = files.filter((f) => /\.upf$/i.test(f.name));
  const inp = files.find((f) => !/\.upf$/i.test(f.name));
  if (upf.length) INPUT.pseudos = [...INPUT.pseudos.filter((p) => !upf.some((u) => u.name === p.name)), ...upf];
  if (inp) {
    INPUT.name = inp.name;
    inp.text().then((t) => { editor.value = t; runInspect(); });
  } else runInspect();
  $("#file-names").textContent = [INPUT.name, ...INPUT.pseudos.map((p) => p.name)].join(" · ");
}
$("#file-in").addEventListener("change", (e) => readFiles(e.target.files));
for (const [dropId, handler] of [["#drop", readFiles], ["#drop2", (f) => analyzeFiles(f)]]) {
  const d = $(dropId);
  d.addEventListener("dragover", (e) => { e.preventDefault(); d.classList.add("over"); });
  d.addEventListener("dragleave", () => d.classList.remove("over"));
  d.addEventListener("drop", (e) => { e.preventDefault(); d.classList.remove("over"); handler(e.dataTransfer.files); });
}
editor.addEventListener("input", () => { clearTimeout(inspectTimer); inspectTimer = setTimeout(runInspect, 600); });
$("#load-example").addEventListener("click", async () => {
  const ex = await api("/api/example");
  INPUT = { name: ex.name, pseudos: [] };
  editor.value = ex.text;
  $("#file-names").textContent = ex.name;
  $("#o-kpath").value = "LGXUG";
  runInspect();
});

function formData() {
  const fd = new FormData();
  fd.append("text", editor.value);
  fd.append("name", INPUT.name);
  INPUT.pseudos.forEach((p) => fd.append("pseudos", p, p.name));
  return fd;
}

let LAST_INSPECT = null;
async function runInspect() {
  if (!editor.value.trim()) return;
  let r;
  try { r = await api("/api/inspect", { method: "POST", body: formData() }); }
  catch (e) { toast(e.message); return; }
  LAST_INSPECT = r;
  $("#inspect-empty").hidden = true; $("#inspect").hidden = false;
  const msgs = [
    ...(r.errors || []).map((m) => ["error", "✕", m]),
    ...(r.warnings || []).map((m) => ["warn", "!", m]),
    ...(r.notes || []).map((m) => ["info", "i", m]),
  ];
  if (!r.errors?.length) msgs.unshift(["ok", "✓", "Input looks valid for pw.x."]);
  $("#msgs").innerHTML = msgs.map(([c, ic, m]) => `<div class="msg ${c}"><span class="ic">${ic}</span><span>${esc(m)}</span></div>`).join("");
  const s = r.structure || {};
  const kv = [
    ["Formula", s.formula ? `${s.formula} (${s.natoms} atoms)` : "–"],
    ["Lattice", s.lattice || "–"],
    ["a, b, c (Å)", s.abc ? s.abc.slice(0, 3).map((x) => x.toFixed(4)).join(", ") : "–"],
    ["α, β, γ (°)", s.abc ? s.abc.slice(3).map((x) => x.toFixed(2)).join(", ") : "–"],
    ["Volume", s.volume ? `${s.volume} Å³` : "–"],
    ["Calculation", `${r.calculation || "–"} · prefix '${esc(r.prefix)}'`],
    ["Cutoffs", r.ecutwfc ? `ecutwfc ${r.ecutwfc} Ry · ecutrho ${r.ecutrho} Ry` : "–"],
    ["k-points", r.kpoints ? r.kpoints.text + (r.kpoints.spacing ? ` · spacing ≈ ${r.kpoints.spacing} Å⁻¹` : "") : "–"],
    ["Occupations", `${r.occupations}${r.nspin === 2 ? " · spin-polarised" : ""}${r.noncolin ? " · non-collinear" : ""}`],
    ["Pseudopotentials", (r.pseudo_status || []).map((p) => `${esc(p.name)} <span class="muted small">(${esc(p.source)})</span>`).join("<br>")],
  ];
  $("#kv").innerHTML = kv.map(([k, v]) => `<dt>${k}</dt><dd>${typeof v === "string" && v.includes("<") ? v : esc(v)}</dd>`).join("");
  $("#kv").outerHTML = `<dl class="kv" id="kv">${$("#kv").innerHTML}</dl>`;
  if (r.kpath) {
    $("#o-kpath").placeholder = r.kpath.path;
    $("#kpath-hint").textContent = `${r.kpath.lattice}; points: ${r.kpath.points.join(" ")}`;
  }
  if (r.nscf_grid) $("#o-nscf").placeholder = r.nscf_grid.join(" ");
  plotStructure($("#struct3d"), s);
  $("#run").disabled = !!r.errors?.length;
}

$("#run").addEventListener("click", async () => {
  const opts = {
    bands: $("#o-bands").checked, dos: $("#o-dos").checked, conv: $("#o-conv").checked,
    kpath: $("#o-kpath").value.trim(), kdensity: Number($("#o-kdens").value) || 40,
    nscf_grid: $("#o-nscf").value.trim() ? $("#o-nscf").value.trim().split(/[\s,x×]+/).map(Number) : null,
    ecut_list: $("#o-ecut").value.trim(), k_list: $("#o-klist").value.trim(), conv_thr_mev: Number($("#o-thr").value) || 1,
    np: Number($("#o-np").value) || SYS.physical_cores, download_pseudos: $("#o-dl").checked,
  };
  const fd = formData();
  fd.append("options", JSON.stringify(opts));
  $("#run").disabled = true;
  try {
    const r = await api("/api/jobs", { method: "POST", body: fd });
    rememberRun(r.id);
    location.hash = `job/${r.id}`;
  } catch (e) { toast(e.message, 6000); }
  finally { $("#run").disabled = false; }
});

/* ------------------------------------------------------------------ runs */
const STATUS_TXT = { done: "done", done_with_errors: "done with errors", failed: "failed", running: "running", queued: "queued", cancelled: "cancelled" };
async function loadRuns() {
  const list = await api("/api/jobs" + (SYS.public ? "?ids=" + encodeURIComponent(myRuns().join(",")) : ""));
  $("#runs-empty").hidden = list.length > 0;
  $("#runs-tbl tbody").innerHTML = list.map((j) => `<tr data-id="${esc(j.id)}">
    <td>${esc(j.id)}</td><td>${esc(j.input_name)}</td>
    <td><span class="badge ${j.status}">${STATUS_TXT[j.status] || j.status}${j.status === "running" && j.current_step ? " · " + esc(j.current_step) : ""}</span></td>
    <td class="num">${fmt(j.energy_ry, 6)}</td><td class="num">${j.is_metal ? "metal" : fmt(j.gap_ev, 3)}</td>
    <td>${j.created ? new Date(j.created * 1000).toLocaleString() : ""}</td></tr>`).join("");
  $$("#runs-tbl tbody tr").forEach((tr) => tr.addEventListener("click", () => (location.hash = `job/${tr.dataset.id}`)));
  if (current.view === "runs" && list.some((j) => j.status === "running" || j.status === "queued")) setTimeout(() => current.view === "runs" && loadRuns(), 3000);
}

/* ------------------------------------------------------------------ job view */
function stopPoll() { clearTimeout(current.poll); current.poll = null; }
let JOB = null, RESULTS_FOR = null;
async function openJob(id) {
  stopPoll();
  if (current.job !== id) { RESULTS_FOR = null; $("#results").hidden = true; }
  current.job = id;
  show("job");
  let st;
  try { st = await api(`/api/jobs/${encodeURIComponent(id)}`); }
  catch (e) { toast(e.message); current.poll = setTimeout(() => current.job === id && openJob(id), 3000); return; }
  if (!st.steps) { current.poll = setTimeout(() => current.job === id && openJob(id), 1000); return; }
  JOB = st;
  $("#job-title").textContent = `${st.prefix} · ${st.id}`;
  const b = $("#job-status");
  b.className = `badge ${st.status}`; b.textContent = STATUS_TXT[st.status] || st.status;
  const active = st.status === "running" || st.status === "queued";
  $("#job-cancel").disabled = !active; $("#job-delete").disabled = active; $("#job-reanalyse").disabled = active;
  $("#job-steps").innerHTML = st.steps.map((s) => {
    const ic = { done: "✓", failed: "✕", running: "●", skipped: "–", pending: "○" }[s.status] || "○";
    return `<li class="${s.status}"><span class="ic">${ic}</span><span>${esc(s.label)}${s.error ? `<span class="err">${esc(s.error).slice(0, 600)}</span>` : ""}</span>
      <span class="t">${s.seconds !== undefined ? s.seconds.toFixed(1) + " s" : s.status === "running" ? "running…" : ""}</span></li>`;
  }).join("") + (st.error ? `<li class="failed"><span class="ic">✕</span><span class="err">${esc(st.error)}</span></li>` : "");
  $("#job-log").textContent = st.log || "";
  $("#job-live").hidden = !st.live; $("#job-live").textContent = st.live || "";
  if (st.results && RESULTS_FOR !== st.id + st.finished) { RESULTS_FOR = st.id + st.finished; renderResults(st); }
  if (active) current.poll = setTimeout(() => current.view === "job" && current.job === id && openJob(id), 1500);
  else loadStatus();
}
$("#back-runs").addEventListener("click", () => (location.hash = "runs"));
$("#job-cancel").addEventListener("click", async () => { await api(`/api/jobs/${current.job}/cancel`, { method: "POST" }); toast("Cancelling…"); });
$("#job-open").addEventListener("click", async () => { try { await api(`/api/jobs/${current.job}/open`, { method: "POST" }); } catch (e) { toast(e.message); } });
$("#job-zip").addEventListener("click", () => (location.href = `/api/jobs/${current.job}/zip`));
$("#job-reanalyse").addEventListener("click", async () => { await api(`/api/jobs/${current.job}/reanalyse`, { method: "POST" }); RESULTS_FOR = null; openJob(current.job); toast("Re-analysed"); });
$("#job-delete").addEventListener("click", async () => {
  if (!confirm(`Delete run ${current.job} and all its files? This cannot be undone.`)) return;
  await api(`/api/jobs/${current.job}`, { method: "DELETE" });
  location.hash = "runs";
});

function tile(label, value, detail) {
  return `<div class="tile"><div class="l">${label}</div><div class="v">${value}</div><div class="d">${detail || ""}</div></div>`;
}
function renderResults(st) {
  const R = st.results, S = R.summary || {};
  $("#results").hidden = false;
  const gb = S.gap_bands || {};
  const gapTile = S.is_metal ? tile("Band gap", "Metal", `E<sub>F</sub> = ${fmt(S.fermi_nscf_ev ?? S.fermi_ev, 3)} eV`)
    : S.gap_ev !== undefined && S.gap_ev !== null ? tile("Band gap", `${S.gap_ev.toFixed(3)} <small>eV</small>`,
      `${S.gap_type || ""}${gb.vbm_where ? ` · VBM ${esc(gb.vbm_where)}, CBM ${esc(gb.cbm_where)}` : ""}${S.direct_gap_ev ? ` · direct gap ${S.direct_gap_ev.toFixed(3)} eV` : ""}`)
      : tile("Band gap", "–", "needs a band structure or DOS run");
  const nat = R.structure?.natoms || 1;
  $("#tiles").innerHTML = [
    gapTile,
    tile("Total energy", `${fmt(S.energy_ry, 6)} <small>Ry</small>`, `${fmt(S.energy_ev, 4)} eV · ${fmt(S.energy_ev_atom, 4)} eV/atom`),
    S.homo_ev !== null && S.homo_ev !== undefined
      ? tile("HOMO / LUMO (SCF grid)", `${fmt(S.homo_ev, 3)} / ${fmt(S.lumo_ev, 3)} <small>eV</small>`, S.lumo_ev ? `HOMO–LUMO gap on SCF k-grid ${(S.lumo_ev - S.homo_ev).toFixed(3)} eV` : "")
      : tile("Fermi energy", `${fmt(S.fermi_ev, 3)} <small>eV</small>`, S.occupations),
    tile("SCF", `${S.scf_iterations ?? "–"} <small>iterations</small>`, `${fmt(S.wall_s, 1)} s wall · ${S.mpi_procs || 1} MPI · QE ${esc(S.qe_version)}`),
    tile("Pressure", `${fmt(S.pressure_kbar, 2)} <small>kbar</small>`, `total force ${fmt(S.total_force, 5)} Ry/bohr`),
    ...(S.total_mag !== null && S.total_mag !== undefined ? [tile("Magnetisation", `${fmt(S.total_mag, 2)} <small>μB/cell</small>`, `absolute ${fmt(S.abs_mag, 2)} μB/cell`)] : []),
  ].join("");
  const checks = [...(R.checks || [])];
  if (!checks.length) checks.push({ level: "ok", msg: "All outputs end with JOB DONE and every SCF converged." });
  $("#checks").innerHTML = checks.map((c) => `<div class="msg ${c.level}"><span class="ic">${{ error: "✕", warn: "!", info: "i", ok: "✓" }[c.level]}</span><span>${esc(c.msg)}</span></div>`).join("") + "<div style='height:8px'></div>";
  drawBandsAndDos();
  // SCF
  plotScf($("#p-scf"), [{ name: "scf", cycles: R.scf?.scf_cycles || [] }]);
  // convergence
  const hasConv = R.conv_ecut || R.conv_k;
  $("#conv-row").hidden = !hasConv;
  if (R.conv_ecut) { plotConv($("#p-conv-ecut"), R.conv_ecut, "ecutwfc (Ry)"); $("#t-conv-ecut").innerHTML = convTable(R.conv_ecut); }
  if (R.conv_k) { plotConv($("#p-conv-k"), R.conv_k, "k-point grid"); $("#t-conv-k").innerHTML = convTable(R.conv_k); }
  plotStructure($("#p-struct"), R.structure);
  // details
  const st0 = R.structure || {};
  const det = [
    ["Formula", `${st0.formula} · ${nat} atoms · V = ${fmt(st0.volume_A3, 3)} Å³`],
    ["Cell a, b, c (Å)", st0.abc ? st0.abc.slice(0, 3).map((x) => x.toFixed(4)).join(", ") : "–"],
    ["XC functional", S.xc || "–"],
    ["Cutoffs", `${S.ecutwfc} / ${S.ecutrho} Ry`],
    ["Electrons", S.nelec], ["Irreducible k-points (SCF)", S.nks],
    ["Energy reference in plots", `${S.energy_reference_label} = ${fmt(S.energy_reference_ev, 4)} eV`],
    ...(S.gap_scf_grid && S.gap_scf_grid.gap_ev !== undefined ? [["Gap on SCF grid", `${fmt(S.gap_scf_grid.gap_ev, 3)} eV`]] : []),
    ...(S.gap_nscf && S.gap_nscf.gap_ev !== undefined ? [["Gap on NSCF grid", `${fmt(S.gap_nscf.gap_ev, 3)} eV`]] : []),
    ...(S.gap_bands && S.gap_bands.gap_ev !== undefined ? [["Gap on band path", `${fmt(S.gap_bands.gap_ev, 3)} eV`]] : []),
  ];
  let html = `<dl class="kv">${det.map(([k, v]) => `<dt>${k}</dt><dd>${esc(v)}</dd>`).join("")}</dl>`;
  const low = R.projwfc?.lowdin || [];
  if (low.length) {
    const keys = [...new Set(low.flatMap((a) => Object.keys(a).filter((k) => !["atom", "charge"].includes(k) && k.length === 1)))];
    html += `<div class="muted small" style="margin:8px 0 4px">Löwdin charges (spilling ${fmt(R.projwfc.spilling, 4)})</div><table class="tbl"><thead><tr><th>Atom</th><th class="num">charge</th>${keys.map((k) => `<th class="num">${k}</th>`).join("")}</tr></thead><tbody>` +
      low.map((a) => `<tr><td>${a.atom} ${esc(st0.symbols ? st0.symbols[a.atom - 1] : "")}</td><td class="num">${fmt(a.charge, 4)}</td>${keys.map((k) => `<td class="num">${fmt(a[k], 4)}</td>`).join("")}</tr>`).join("") + "</tbody></table>";
  }
  $("#details").innerHTML = html;
  // files
  $("#win-dir").textContent = st.win_dir;
  $("#files").innerHTML = st.files.filter((f) => !f.startsWith("pseudos/")).map((f) =>
    `<a href="/api/jobs/${encodeURIComponent(st.id)}/file/${f.split("/").map(encodeURIComponent).join("/")}" target="_blank">${esc(f)}</a>`).join("");
}
function drawBandsAndDos() {
  if (!JOB || !JOB.results) return;
  const R = JOB.results, S = R.summary || {};
  const ref = S.energy_reference_ev ?? 0, refL = S.energy_reference_label || "E_F";
  const emin = Number($("#e-min").value), emax = Number($("#e-max").value);
  if (R.bands) {
    plotBands($("#p-bands"), R.bands, ref, refL, R.dos, R.pdos, emin, emax, $("#show-pdos").checked);
    $("#bands-note").textContent = `Energies from pw.x (${R.bands.nbnd} bands, ${R.bands.x.length} k-points), shifted so that ${refL} = 0. ` +
      (R.bandsx_gnu ? `bands.x wrote ${R.bandsx_gnu.file} too.` : "");
  } else $("#p-bands").innerHTML = '<div class="empty">Band structure was not calculated.</div>';
  if (R.pdos || R.dos) plotPdos($("#p-pdos"), R.pdos, R.dos, ref, refL, emin, emax);
  else $("#p-pdos").innerHTML = '<div class="empty">DOS was not calculated.</div>';
}
["#e-min", "#e-max", "#show-pdos"].forEach((s) => $(s).addEventListener("change", drawBandsAndDos));

/* ------------------------------------------------------------------ analyse outputs */
$("#files-an").addEventListener("change", (e) => analyzeFiles(e.target.files));
$("#an-go").addEventListener("click", async () => {
  const fd = new FormData(); fd.append("folder", $("#an-folder").value);
  try { renderAnalysis(await api("/api/analyze", { method: "POST", body: fd })); } catch (e) { toast(e.message, 6000); }
});
async function analyzeFiles(list) {
  const fd = new FormData();
  [...list].forEach((f) => fd.append("files", f, f.name));
  try { renderAnalysis(await api("/api/analyze", { method: "POST", body: fd })); } catch (e) { toast(e.message, 6000); }
}
const AN_COLS = [
  ["file", "File"], ["program", "Program"], ["calculation", "Type"], ["job_done", "JOB DONE"], ["converged", "SCF converged"],
  ["energy_ry", "Energy (Ry)", 6], ["fermi_ev", "E_F (eV)", 3], ["homo_ev", "HOMO (eV)", 3], ["lumo_ev", "LUMO (eV)", 3], ["gap", "Gap (eV)", 3],
  ["total_force", "Force (Ry/bohr)", 5], ["pressure_kbar", "P (kbar)", 2], ["total_mag", "Mag (μB)", 2], ["ecutwfc", "ecutwfc", 1],
  ["nks", "k-points", 0], ["scf_iterations", "SCF iter", 0], ["wall_s", "Wall (s)", 1], ["version", "QE"],
];
let AN = null;
function anValue(it, key) {
  if (key === "converged") return it.program === "PWSCF" ? !it.not_converged : null;
  if (key === "gap") return it.gap ? (it.gap.metal ? "metal" : it.gap.gap_ev) : it.homo_lumo_gap_ev ?? null;
  return it[key];
}
function renderAnalysis(res) {
  AN = res;
  const items = res.items || [];
  $("#an-results").hidden = false;
  if (!items.length) { $("#an-tbl").innerHTML = "<tr><td class='empty'>No QE output files found.</td></tr>"; return; }
  $("#an-tbl").innerHTML = `<thead><tr>${AN_COLS.map(([, h, d]) => `<th class="${d !== undefined ? "num" : ""}">${h}</th>`).join("")}</tr></thead><tbody>` +
    items.map((it) => `<tr>${AN_COLS.map(([k, , d]) => {
      const v = anValue(it, k);
      const txt = typeof v === "boolean" ? (v ? "✓" : "✕") : d !== undefined && typeof v === "number" ? fmt(v, d) : esc(v ?? "–");
      return `<td class="${d !== undefined ? "num" : ""}">${txt}</td>`;
    }).join("")}</tr>`).join("") + "</tbody>";
  const checks = [];
  items.forEach((it) => {
    if (!it.job_done) checks.push(["error", `${it.file}: does not end with JOB DONE (crashed, still running or killed).`]);
    if (it.not_converged) checks.push(["error", `${it.file}: convergence NOT achieved.`]);
    (it.errors || []).forEach((e) => checks.push(["error", `${it.file}: ${e}`]));
    if (it.relax && !it.relax.converged) checks.push(["warn", `${it.file}: relaxation did not reach the force/energy criteria.`]);
  });
  if (!checks.length) checks.push(["ok", `All ${items.length} files end with JOB DONE and converged.`]);
  $("#an-checks").innerHTML = "<div style='height:10px'></div>" + checks.map(([c, m]) => `<div class="msg ${c}"><span class="ic">${{ error: "✕", warn: "!", ok: "✓" }[c]}</span><span>${esc(m)}</span></div>`).join("");

  const pw = items.filter((it) => it.program === "PWSCF");
  const withScf = pw.filter((it) => it.scf_cycles && it.scf_cycles.some((c) => c.length));
  if (withScf.length) plotScf($("#a-scf"), withScf.map((it) => ({ name: it.file, cycles: it.scf_cycles })));
  else $("#a-scf").innerHTML = '<div class="empty">No SCF iterations in these files.</div>';

  // energy comparison: group by (nat, nelec), biggest group
  const groups = {};
  pw.filter((it) => it.energy_ry !== null && it.energy_ry !== undefined).forEach((it) => { const k = `${it.nat}|${it.nelec}`; (groups[k] = groups[k] || []).push(it); });
  const g = Object.values(groups).sort((a, b) => b.length - a.length)[0] || [];
  if (g.length >= 2) {
    const nat = g[0].nat || 1;
    const emin = Math.min(...g.map((it) => it.final_energy_ry ?? it.energy_ry));
    Plotly.react($("#a-energy"), [{ type: "bar", x: g.map((it) => it.file), y: g.map((it) => ((it.final_energy_ry ?? it.energy_ry) - emin) * RY * 1000 / nat),
      marker: { color: css("--s1") }, customdata: g.map((it) => it.final_energy_ry ?? it.energy_ry),
      hovertemplate: "%{x}<br>E = %{customdata:.6f} Ry<br>ΔE = %{y:.2f} meV/atom<extra></extra>" }],
      layout({ yaxis: axis({ title: `E − E_min (meV/atom), ${nat} atoms, ${g[0].nelec} e⁻` }), xaxis: axis({ automargin: true }), bargap: 0.4 }), CONFIG);
  } else $("#a-energy").innerHTML = '<div class="empty">Needs at least two total energies for the same system (same atoms and electrons) to compare.</div>';

  // relaxation
  const rel = pw.filter((it) => it.relax && it.relax.energy_ry.length > 1);
  $("#a-relax-card").hidden = !rel.length;
  if (rel.length) {
    const col = SERIES();
    Plotly.react($("#a-relax"), rel.slice(0, 8).map((it, i) => {
      const E = it.relax.energy_ry, ef = E[E.length - 1];
      return { x: E.map((_, j) => j + 1), y: E.map((e) => (e - ef) * RY * 1000), name: it.file, mode: "lines+markers",
        line: { color: col[i], width: 2 }, marker: { size: 8 }, customdata: E.map((_, j) => it.relax.force[j] ?? null),
        hovertemplate: `${esc(it.file)}<br>step %{x}<br>E − E_final = %{y:.2f} meV<br>total force %{customdata:.5f} Ry/bohr<extra></extra>` };
    }), layout({ xaxis: axis({ title: "ionic step" }), yaxis: axis({ title: "E − E_final (meV/cell)" }), showlegend: rel.length > 1 }), CONFIG);
  }
  // ecut / k convergence among files of one group
  const ecuts = new Set(g.map((it) => it.ecutwfc));
  const nks = new Set(g.map((it) => it.nks));
  $("#a-conv-card").hidden = !(g.length >= 3 && (ecuts.size >= 3 || nks.size >= 3));
  if (!$("#a-conv-card").hidden) {
    const byEcut = ecuts.size >= nks.size;
    const pts = [...g].sort((a, b) => (byEcut ? a.ecutwfc - b.ecutwfc : a.nks - b.nks));
    const ref = pts[pts.length - 1].energy_ry, nat = pts[0].nat || 1;
    const rows = pts.map((it) => ({ value: byEcut ? it.ecutwfc : it.nks, label: `${it.file} = ${byEcut ? it.ecutwfc : it.nks}`, delta_mev_atom: (it.energy_ry - ref) * RY * 1000 / nat }));
    plotConv($("#a-conv"), { rows, threshold_mev_atom: 1 }, byEcut ? "ecutwfc (Ry)" : "irreducible k-points");
  }
  // band structures
  const bandItems = items.filter((it) => it.band_plot);
  $("#a-bands").innerHTML = bandItems.map((it, i) => `<div class="card bandcard"><div class="card-h"><h2>Bands · ${esc(it.file)}</h2>
    <span class="muted small">${it.band_plot.edges && it.band_plot.edges.metal === false ? `gap ${it.band_plot.edges.gap_ev.toFixed(3)} eV (${it.band_plot.edges.gap_type}) · VBM ${esc(it.band_plot.edges.vbm_where || "")}, CBM ${esc(it.band_plot.edges.cbm_where || "")}` : ""}
    ${it.band_plot.ticks.some((t) => t.label) ? "" : " · add the bands .in file to get k-point labels"}</span></div><div class="plot" id="ab-${i}"></div></div>`).join("");
  bandItems.forEach((it, i) => {
    const bp = it.band_plot;
    plotBands($(`#ab-${i}`), { x: bp.x, eig: bp.eig, ticks: bp.ticks, edges: bp.edges }, bp.reference_ev ?? 0, bp.reference_label, null, null, -6, 6, false);
  });
}
$("#an-csv").addEventListener("click", () => {
  if (!AN) return;
  const rows = [AN_COLS.map(([, h]) => h), ...AN.items.map((it) => AN_COLS.map(([k]) => { const v = anValue(it, k); return v === null || v === undefined ? "" : v; }))];
  const csv = rows.map((r) => r.map((c) => `"${String(c).replace(/"/g, '""')}"`).join(",")).join("\n");
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([csv], { type: "text/csv" }));
  a.download = "qe_summary.csv"; a.click();
});

/* ------------------------------------------------------------------ boot */
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => {
  if (LAST_INSPECT) plotStructure($("#struct3d"), LAST_INSPECT.structure);
  if (JOB && JOB.results) { RESULTS_FOR = null; renderResults(JOB); }
  if (AN) renderAnalysis(AN);
});
loadStatus();
setInterval(loadStatus, 10000);
route();
