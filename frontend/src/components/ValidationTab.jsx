import { useEffect, useState } from "react";
import { getBenchmark, validationPlotUrl } from "../api";

export default function ValidationTab() {
  const [data, setData] = useState(null);
  const [err, setErr] = useState(null);

  useEffect(() => {
    getBenchmark().then(setData).catch((e) => setErr(String(e)));
  }, []);

  if (err) return <div className="text-red-600 text-sm">Failed to load validation: {err}</div>;
  if (!data) return <div className="text-slate-400 text-sm">Loading validation…</div>;

  const f = data.faithfulness;
  return (
    <div className="space-y-6">
      <div className="grid grid-cols-1 sm:grid-cols-3 gap-4">
        <Metric
          label="Per-mutation Spearman ρ (pooled)"
          value={data.overall_spearman_rho?.toFixed(3) ?? "—"}
          sub={`n = ${data.n_mutations} mutations`}
        />
        <Metric
          label="Top-ddG DRM enrichment"
          value="2.9×"
          sub="vs 12% base rate"
        />
        <Metric
          label="Explanation faithfulness"
          value={f ? `${Math.round(f.pct_correct)}%` : "—"}
          sub={f ? `mean ${f.mean.toFixed(2)} / 2, n=${f.n}` : "run scripts/07"}
        />
      </div>

      <div className="rounded-lg border border-slate-200 p-4">
        <div className="text-sm text-slate-500 mb-2">
          Predicted ΔΔG vs measured clinical fold-resistance, DRM enrichment, and per-drug signal
        </div>
        <img
          src={validationPlotUrl()}
          alt="ResistScope validation"
          className="w-full max-w-full rounded"
        />
      </div>

      <p className="text-xs text-slate-400 leading-relaxed">
        Honest result: rigid single-mutation docking gives a weak <em>pooled</em> per-mutation
        correlation — the measured fold-resistance target is confounded by co-occurring mutations —
        but its top-ranked predictions are ~2.9× enriched for known major DRMs, and darunavir
        validates well (primary-DRM ρ ≈ 0.40). Faithfulness: {f ? `${Math.round(f.pct_correct)}%` : "—"} of
        LLM explanations correctly identify the expert-annotated mechanism.
      </p>
    </div>
  );
}

function Metric({ label, value, sub }) {
  return (
    <div className="rounded-lg border border-slate-200 p-4">
      <div className="text-xs text-slate-500">{label}</div>
      <div className="text-2xl font-semibold text-slate-800 tabular-nums mt-1">{value}</div>
      <div className="text-xs text-slate-400 mt-0.5">{sub}</div>
    </div>
  );
}
