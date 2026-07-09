function scoreColor(s) {
  if (s == null) return "text-slate-400";
  if (s >= 70) return "text-emerald-600";
  if (s >= 40) return "text-amber-500";
  return "text-red-600";
}
function scoreLabel(s) {
  if (s == null) return "—";
  if (s >= 70) return "Robust";
  if (s >= 40) return "Moderate";
  return "Vulnerable";
}
const fmt = (v, d = 2) => (v == null ? "—" : v.toFixed(d));

export default function ScoreCard({ result }) {
  const s = result.robustness_score;
  const sm = result.scoring_methods || {};
  return (
    <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
      <div className="md:col-span-1 rounded-lg border border-slate-200 p-5 flex flex-col items-center justify-center">
        <div className={`text-6xl font-semibold tabular-nums ${scoreColor(s)}`}>
          {s == null ? "—" : Math.round(s)}
        </div>
        <div className="text-sm text-slate-500 mt-1">robustness / 100</div>
        <div className={`mt-2 text-sm font-medium ${scoreColor(s)}`}>{scoreLabel(s)}</div>
      </div>

      <div className="md:col-span-2 rounded-lg border border-slate-200 p-5">
        <div className="grid grid-cols-2 gap-y-3 gap-x-6 text-sm">
          <Stat label="Wildtype binding" value={`${fmt(result.wildtype_binding)} kcal/mol`} />
          <Stat label="Mutations scored" value={result.n_mutations_scored ?? "—"} mono />
          <Stat label="Mean ΔΔG" value={`${fmt(sm.simple_mean)} kcal/mol`} />
          <Stat label="Prevalence-weighted ΔΔG" value={`${fmt(sm.weighted)} kcal/mol`} />
          <Stat label="Worst-case ΔΔG" value={`${fmt(sm.worst_case)} kcal/mol`} />
          <Stat label="Failed docks" value={result.n_mutations_failed ?? 0} mono />
        </div>
        <p className="text-xs text-slate-400 mt-4">
          ΔΔG = mutant − wildtype binding. Positive = weaker binding = predicted resistance.
        </p>
      </div>
    </div>
  );
}

function Stat({ label, value, mono }) {
  return (
    <div>
      <div className="text-slate-500">{label}</div>
      <div className={`text-slate-800 font-medium ${mono ? "font-mono" : ""}`}>{value}</div>
    </div>
  );
}
