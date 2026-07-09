import { useEffect, useState } from "react";
import { listDrugs, getDrug, triage } from "./api";
import ScoreCard from "./components/ScoreCard";
import MutationTable from "./components/MutationTable";
import ValidationTab from "./components/ValidationTab";

export default function App() {
  const [drugs, setDrugs] = useState([]);
  const [selected, setSelected] = useState("");
  const [smiles, setSmiles] = useState("");
  const [result, setResult] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [tab, setTab] = useState("triage");

  useEffect(() => {
    listDrugs().then((d) => setDrugs(d.drugs)).catch(() => {});
  }, []);

  // Selecting a benchmark drug fills SMILES and loads instant precomputed results.
  async function onSelectDrug(abbrev) {
    setSelected(abbrev);
    setError(null);
    const drug = drugs.find((d) => d.abbrev === abbrev);
    setSmiles(drug?.smiles || "");
    if (!abbrev) return;
    setLoading(true);
    try {
      setResult(await getDrug(abbrev));
    } catch (e) {
      setError(String(e));
    } finally {
      setLoading(false);
    }
  }

  // Triage the current input: instant if a benchmark drug is selected,
  // otherwise a live dock of the custom SMILES (minutes).
  async function onTriage() {
    if (!smiles.trim()) return;
    setError(null);
    setLoading(true);
    setResult(null);
    try {
      setResult(selected ? await getDrug(selected) : await triage(smiles.trim()));
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
        <div className="max-w-5xl mx-auto px-6 py-4 flex items-baseline justify-between">
          <div>
            <h1 className="text-xl font-semibold text-slate-800">
              Resist<span className="text-teal-600">Scope</span>
            </h1>
            <p className="text-sm text-slate-500">Resistance-aware triage for HIV-1 protease inhibitors</p>
          </div>
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
      </header>

      <main className="max-w-5xl mx-auto px-6 py-6">
        {tab === "validation" ? (
          <ValidationTab />
        ) : (
          <div className="space-y-6">
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
                  Benchmark drugs load instantly. A custom SMILES runs live AutoDock Vina
                  (~2–5 min on CPU).
                </span>
              </div>
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
                <MutationTable mutations={result.mutations} />
              </div>
            )}
          </div>
        )}
      </main>

      <footer className="max-w-5xl mx-auto px-6 py-6 text-xs text-slate-400 border-t border-slate-200 mt-8">
        ResistScope · Built with Claude: Life Sciences · docking = AutoDock Vina / Uni-Dock,
        explanations = Claude (Haiku 4.5)
      </footer>
    </div>
  );
}
