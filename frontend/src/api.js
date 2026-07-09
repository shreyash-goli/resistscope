// Thin client for the ResistScope FastAPI backend.
// In dev, Vite proxies /api/* -> http://localhost:8000 (see vite.config.js).
const BASE = "/api";

async function getJSON(path) {
  const r = await fetch(`${BASE}${path}`);
  if (!r.ok) throw new Error((await r.text()) || `HTTP ${r.status}`);
  return r.json();
}

export const health = () => getJSON("/health");
export const listDrugs = () => getJSON("/drugs");
export const getDrug = (abbrev) => getJSON(`/drug/${abbrev}`);
export const getBenchmark = () => getJSON("/benchmark");
export const validationPlotUrl = () => `${BASE}/validation/plot`;

// Live triage of a custom SMILES (slow: minutes of docking).
export async function triage(smiles, subset = "primary") {
  const r = await fetch(`${BASE}/triage`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ smiles, target: "HIV1_PR", subset }),
  });
  if (!r.ok) throw new Error((await r.text()) || `HTTP ${r.status}`);
  return r.json();
}
