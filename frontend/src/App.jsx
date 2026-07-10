import { useEffect, useState } from "react";
import { listTargets, listDrugs, getDrug, triage, health } from "./api";
import ScoreCard from "./components/ScoreCard";
import MutationTable from "./components/MutationTable";
import ValidationTab from "./components/ValidationTab";
import AddTargetModal from "./components/AddTargetModal";

export default function App() {
  const [targets, setTargets] = useState([]);
  const [target, setTarget] = useState("HIV1_PR");
  const [drugs, setDrugs] = useState([]);
  const [selected, setSelected] = useState("");
  const [smiles, setSmiles] = useState("");
  const [result, setResult] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [tab, setTab] = useState("triage");
  const [live, setLive] = useState(null); // {live_docking, docking_backend}
  const [showAdd, setShowAdd] = useState(false);

  async function onTargetSaved(name) {
    setShowAdd(false);
    try {
      const d = await listTargets();
      setTargets(d.targets);
    } catch { /* keep old list */ }
    setTarget(name);
  }

  useEffect(() => {
    listTargets().then((d) => setTargets(d.targets)).catch(() => {});
    health().then(setLive).catch(() => {});
  }, []);

  // Load the drug list whenever the target changes; reset any prior result.
  useEffect(() => {
    setDrugs([]); setSelected(""); setSmiles(""); setResult(null); setError(null);
    listDrugs(target).then((d) => setDrugs(d.drugs)).catch(() => {});
  }, [target]);

  const activeTarget = targets.find((t) => t.name === target);
  const dockingReady = activeTarget ? activeTarget.docking_ready : true;

  // Selecting a benchmark drug fills SMILES and loads instant precomputed results.
  async function onSelectDrug(abbrev) {
    setSelected(abbrev);
    setError(null);
    const drug = drugs.find((d) => d.abbrev === abbrev);
    setSmiles(drug?.smiles || "");
    if (!abbrev) return;
    setLoading(true);
    try {
      setResult(await getDrug(abbrev, target));
    } catch (e) {
      setError(e.status === 404
        ? `${activeTarget?.label || target}: docking not run yet — see the runbook to dock on GPU.`
        : String(e));
    } finally {
      setLoading(false);
    }
  }

  async function onTriage() {
    if (!smiles.trim()) return;
    setError(null);
    setLoading(true);
    setResult(null);
    try {
      setResult(selected ? await getDrug(selected, target) : await triage(smiles.trim(), target));
    } catch (e) {
      setError(String(e));
    } finally {
      setLoading(false);
    }
  }

  const isBenchmark = !!selected && result?.precomputed;

  return (
    <div className="min-h-screen">
      <header className="border-b border-slate-200">
        <div className="max-w-5xl mx-auto px-6 py-4 flex items-baseline justify-between gap-4">
          <div>
            <h1 className="text-xl font-semibold text-slate-800">
              Resist<span className="text-teal-600">Scope</span>
            </h1>
            <p className="text-sm text-slate-500">
              Resistance-aware triage · {activeTarget?.label || "HIV-1 protease"}
            </p>
          </div>
          <div className="flex items-center gap-3">
            {targets.length > 0 && (
              <div className="flex items-center gap-2">
                <label className="text-xs text-slate-500 flex items-center gap-1.5">
                  Target
                  <select
                    value={target}
                    onChange={(e) => setTarget(e.target.value)}
                    className="border border-slate-300 rounded px-2 py-1 bg-white text-slate-700"
                  >
                    {targets.map((t) => (
                      <option key={t.name} value={t.name}>
                        {t.label}{t.user ? " ★" : ""}{t.docking_ready ? "" : " (docking pending)"}
                      </option>
                    ))}
                  </select>
                </label>
                <button
                  onClick={() => setShowAdd(true)}
                  className="text-xs text-teal-600 hover:text-teal-700 border border-teal-200 hover:border-teal-300 rounded px-2 py-1"
                  title="Bring your own target: upload a PDB, Claude finds the drugs + mutations"
                >
                  + Add
                </button>
              </div>
            )}
            <nav className="flex gap-1 text-sm">
              {["triage", "validation"].map((t) => (
                <button
                  key={t}
                  onClick={() => setTab(t)}
                  className={`px-3 py-1.5 rounded ${
                    tab === t ? "bg-slate-100 text-slate-800 font-medium" : "text-slate-500 hover:text-slate-700"
                  }`}
                >
                  {t === "triage" ? "Triage" : "Validation"}
                </button>
              ))}
            </nav>
          </div>
        </div>
      </header>

      <main className="max-w-5xl mx-auto px-6 py-6">
        {tab === "validation" ? (
          <ValidationTab target={target} targetLabel={activeTarget?.label} dockingReady={dockingReady} />
        ) : (
          <div className="space-y-6">
            {!dockingReady && (
              <div className="rounded-lg border border-amber-200 bg-amber-50 text-amber-800 text-sm px-4 py-3">
                <strong>{activeTarget?.label}</strong> is wired end-to-end (dataset, panels,
                docking box, and literature-grounded mechanisms) but its docking run is
                pending — it runs on GPU (Uni-Dock / A100). Drugs and mechanisms below are
                ready; precomputed ΔΔG results appear once that run completes.
              </div>
            )}
            {/* Input panel */}
            <div className="rounded-lg border border-slate-200 p-5 space-y-4">
              <div className="grid grid-cols-1 sm:grid-cols-3 gap-4">
                <label className="text-sm">
                  <span className="text-slate-500">Benchmark drug (instant)</span>
                  <select
                    value={selected}
                    onChange={(e) => onSelectDrug(e.target.value)}
                    className="mt-1 w-full border border-slate-300 rounded px-2 py-1.5 bg-white"
                  >
                    <option value="">— select —</option>
                    {drugs.map((d) => (
                      <option key={d.abbrev} value={d.abbrev}>
                        {d.name} ({d.abbrev})
                      </option>
                    ))}
                  </select>
                </label>
                <label className="text-sm sm:col-span-2">
                  <span className="text-slate-500">or paste a SMILES (custom compound)</span>
                  <textarea
                    value={smiles}
                    onChange={(e) => {
                      setSmiles(e.target.value);
                      setSelected("");
                    }}
                    rows={2}
                    placeholder="CC(C)CN(C[C@H](...))..."
                    className="mt-1 w-full border border-slate-300 rounded px-2 py-1.5 font-mono text-xs resize-none"
                  />
                </label>
              </div>
              <div className="flex items-center gap-3">
                <button
                  onClick={onTriage}
                  disabled={loading || !smiles.trim()}
                  className="px-4 py-1.5 rounded bg-teal-600 text-white text-sm font-medium disabled:bg-slate-300 hover:bg-teal-700"
                >
                  {selected ? "Show results" : "Triage compound"}
                </button>
                <span className="text-xs text-slate-400">
                  Benchmark drugs load instantly. A custom SMILES runs live docking.
                </span>
              </div>
              {/* Live-docking backend status so the user knows a custom SMILES will work */}
              {live && (
                <div className="flex items-center gap-2 text-xs">
                  <span
                    className={`inline-block w-1.5 h-1.5 rounded-full ${
                      live.live_docking ? "bg-emerald-500" : "bg-amber-500"
                    }`}
                  />
                  {live.live_docking ? (
                    <span className="text-slate-500">
                      Live docking ready via <span className="font-mono">{live.docking_backend}</span>
                      {live.docking_backend?.startsWith("remote") ? " (GPU worker)" : ""}
                    </span>
                  ) : (
                    <span className="text-slate-500">
                      Live docking not configured — benchmark drugs work now; to dock custom SMILES,
                      run the GPU worker and set <span className="font-mono">RESISTSCOPE_DOCKING_URL</span> (see README).
                    </span>
                  )}
                </div>
              )}
            </div>

            {error && (
              <div className="rounded border border-red-200 bg-red-50 text-red-700 text-sm px-4 py-3">
                {error}
              </div>
            )}

            {loading && (
              <div className="rounded-lg border border-slate-200 p-8 text-center text-slate-500">
                <div className="animate-pulse">
                  {selected ? "Loading precomputed results…" : "Docking against the resistance panel — this can take a few minutes…"}
                </div>
              </div>
            )}

            {result && !loading && (
              <div className="space-y-5">
                {isBenchmark && (
                  <div className="text-xs text-slate-400">
                    Showing precomputed benchmark results (GPU Uni-Dock run).
                  </div>
                )}
                <ScoreCard result={result} />
                <MutationTable mutations={result.mutations} target={target} />
              </div>
            )}
          </div>
        )}
      </main>

      <footer className="max-w-5xl mx-auto px-6 py-6 text-xs text-slate-400 border-t border-slate-200 mt-8">
        ResistScope · Built with Claude: Life Sciences · docking = AutoDock Vina / Uni-Dock,
        explanations = Claude (Haiku 4.5)
      </footer>

      {showAdd && <AddTargetModal onClose={() => setShowAdd(false)} onSaved={onTargetSaved} />}
    </div>
  );
}
