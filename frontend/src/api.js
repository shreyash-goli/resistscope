// Thin client for the ResistScope FastAPI backend.
// In dev, Vite proxies /api/* -> http://localhost:8000 (see vite.config.js).
const BASE = "/api";

async function getJSON(path) {
  const r = await fetch(`${BASE}${path}`);
  if (!r.ok) {
    const err = new Error((await r.text()) || `HTTP ${r.status}`);
    err.status = r.status;
    throw err;
  }
  return r.json();
}

async function postJSON(path, body) {
  const r = await fetch(`${BASE}${path}`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  if (!r.ok) {
    const err = new Error((await r.text()) || `HTTP ${r.status}`);
    err.status = r.status;
    throw err;
  }
  return r.json();
}

// Most endpoints are target-aware (?target=HIV1_PR | HIV1_RT). Default PI.
const q = (target) => (target ? `?target=${encodeURIComponent(target)}` : "");

export const health = (target) => getJSON(`/health${q(target)}`);
export const listTargets = () => getJSON("/targets");
export const listDrugs = (target) => getJSON(`/drugs${q(target)}`);
export const getDrug = (abbrev, target) => getJSON(`/drug/${abbrev}${q(target)}`);
export const getBenchmark = (target) => getJSON(`/benchmark${q(target)}`);
export const validationPlotUrl = (target) => `${BASE}/validation/plot${q(target)}`;

// Bring-your-own-target flow.
export const intakeTarget = (body) => postJSON("/targets/intake", body);
export const assembleTarget = (body) => postJSON("/targets/assemble", body);
export const saveTarget = (body) => postJSON("/targets/save", body);

// Live triage of a custom SMILES (slow: minutes of docking).
export async function triage(smiles, target = "HIV1_PR", subset = "primary") {
  const r = await fetch(`${BASE}/triage`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ smiles, target, subset }),
  });
  if (!r.ok) throw new Error((await r.text()) || `HTTP ${r.status}`);
  return r.json();
}
