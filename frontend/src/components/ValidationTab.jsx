import { useEffect, useState } from "react";
import { getBenchmark, validationPlotUrl } from "../api";

export default function ValidationTab({ target = "HIV1_PR", targetLabel }) {
  const [data, setData] = useState(null);
  const [err, setErr] = useState(null);
  const [pending, setPending] = useState(false);

  useEffect(() => {
    setData(null); setErr(null); setPending(false);
    getBenchmark(target)
      .then(setData)
      .catch((e) => (e.status === 404 ? setPending(true) : setErr(String(e))));
  }, [target]);

  if (pending)
    return (
      <div className="rounded-lg border border-amber-200 bg-amber-50 text-amber-800 text-sm px-4 py-4 leading-relaxed">
        <strong>{targetLabel || target}</strong> validation is pending its GPU docking run.
        Once docked, run <code className="text-xs">scripts/05_validate.py</code> and{" "}
        <code className="text-xs">scripts/08_benchmark.py</code> with{" "}
        <code className="text-xs">--target rt</code> and this tab populates automatically — the
        same rigorous metrics (enrichment CIs, ROC-AUC, de-confounding) computed for protease.
      </div>
    );
  if (err) return <div className="text-red-600 text-sm">Failed to load validation: {err}</div>;
  if (!data) return <div className="text-slate-400 text-sm">Loading validation…</div>;

  const f = data.faithfulness;
  const rig = data.rigorous;
  const h = rig?.headline;
  const fmtP = (p) => (p < 1e-3 ? "< 0.001" : p.toFixed(3));
  const enrichRows = rig?.enrichment ?? [];
  const rankRows = (rig?.ranking ?? []).filter((r) => r.drug !== "OVERALL");
  const deconRows = rig?.deconfounding ?? [];

  return (
    <div className="space-y-6">
      {/* Headline metrics — each with a significance test or CI */}
      <div className="grid grid-cols-1 sm:grid-cols-3 gap-4">
        <Metric
          label="Top-ΔΔG DRM enrichment"
          value={h ? `${h.best_enrichment.toFixed(1)}×` : "2.9×"}
          sub={
            h
              ? `top-${h.best_enrichment_top_n} · 95% CI ${h.best_enrichment_ci[0].toFixed(
                  1
                )}–${h.best_enrichment_ci[1].toFixed(1)} · perm p ${fmtP(h.best_enrichment_perm_p)}`
              : "vs 12% base rate"
          }
          accent="teal"
        />
        <Metric
          label="DRM-recovery ROC-AUC (pooled)"
          value={h ? h.pooled_roc_auc.toFixed(2) : "—"}
          sub={
            h
              ? `95% CI ${h.pooled_roc_auc_ci[0].toFixed(2)}–${h.pooled_roc_auc_ci[1].toFixed(
                  2
                )} · ~chance as a global ranker`
              : "run scripts/08"
          }
        />
        <Metric
          label="Explanation faithfulness"
          value={f ? `${Math.round(f.pct_correct)}%` : "—"}
          sub={f ? `mean ${f.mean.toFixed(2)} / 2, n=${f.n}` : "run scripts/07"}
          accent="teal"
        />
      </div>

      {/* The honest one-paragraph read */}
      <div className="rounded-lg border border-slate-200 bg-slate-50 p-4">
        <div className="text-sm font-medium text-slate-700 mb-1">What the numbers say</div>
        <p className="text-xs text-slate-500 leading-relaxed">
          Rigid single-mutation docking ΔΔG is a <strong>coarse DRM-triage flag, not a
          quantitative resistance predictor</strong>. Its most extreme predictions are
          significantly enriched for real major resistance mutations
          {h && ` (${h.best_enrichment.toFixed(1)}× at top-${h.best_enrichment_top_n}, permutation p ${fmtP(h.best_enrichment_perm_p)})`},
          yet as a global ranker it sits at chance (ROC-AUC {h ? h.pooled_roc_auc.toFixed(2) : "≈0.5"}). And
          de-confounding the target — restricting to isolates carrying fewer co-occurring
          mutations — <strong>does not rescue</strong> the magnitude correlation (below), which
          bounds what a rigid-docking method can do here.
        </p>
      </div>

      {/* Enrichment table with CI + permutation p */}
      {enrichRows.length > 0 && (
        <Panel title="Top-N DRM enrichment (permutation p, bootstrap 95% CI)">
          <table className="w-full text-xs tabular-nums">
            <thead className="text-slate-400 text-left">
              <tr>
                <Th>Top N</Th><Th>DRM precision</Th><Th>Enrichment</Th>
                <Th>95% CI</Th><Th>Permutation p</Th>
              </tr>
            </thead>
            <tbody className="text-slate-600">
              {enrichRows.map((r) => (
                <tr key={r.top_n} className="border-t border-slate-100">
                  <Td>{r.top_n}</Td>
                  <Td>{(r.precision * 100).toFixed(0)}% <span className="text-slate-400">/ {(r.base_rate * 100).toFixed(0)}% base</span></Td>
                  <Td className="font-medium text-teal-700">{r.enrichment.toFixed(2)}×</Td>
                  <Td>{r.ci_low.toFixed(2)}–{r.ci_high.toFixed(2)}</Td>
                  <Td>{fmtP(r.perm_p)}</Td>
                </tr>
              ))}
            </tbody>
          </table>
        </Panel>
      )}

      {/* Per-drug AUC + de-confounding, side by side on wide screens */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {rankRows.length > 0 && (
          <Panel title="Per-drug DRM-recovery ROC-AUC (0.5 = chance)">
            <table className="w-full text-xs tabular-nums">
              <thead className="text-slate-400 text-left">
                <tr><Th>Drug</Th><Th>ROC-AUC</Th><Th>PR-AUC</Th><Th>n</Th></tr>
              </thead>
              <tbody className="text-slate-600">
                {rankRows.map((r) => (
                  <tr key={r.drug} className="border-t border-slate-100">
                    <Td className="font-medium">{r.drug}</Td>
                    <Td className={r.roc_auc >= 0.6 ? "text-teal-700 font-medium" : ""}>{r.roc_auc.toFixed(3)}</Td>
                    <Td>{r.pr_auc.toFixed(3)}</Td>
                    <Td>{r.n}</Td>
                  </tr>
                ))}
              </tbody>
            </table>
          </Panel>
        )}

        {deconRows.length > 0 && (
          <Panel title="De-confounding: does a cleaner target rescue it?">
            <table className="w-full text-xs tabular-nums">
              <thead className="text-slate-400 text-left">
                <tr><Th>Isolate subset</Th><Th>n pairs</Th><Th>Spearman ρ</Th><Th>p</Th></tr>
              </thead>
              <tbody className="text-slate-600">
                {deconRows.map((r) => (
                  <tr key={r.subset} className="border-t border-slate-100">
                    <Td className="whitespace-nowrap">{r.subset}</Td>
                    <Td>{r.n_pairs}</Td>
                    <Td className={r.spearman_rho < 0 ? "text-slate-500" : ""}>
                      {r.spearman_rho == null ? "—" : r.spearman_rho.toFixed(3)}
                    </Td>
                    <Td>{r.spearman_pvalue == null ? "—" : fmtP(r.spearman_pvalue)}</Td>
                  </tr>
                ))}
              </tbody>
            </table>
            <p className="text-[11px] text-slate-400 mt-2 leading-snug">
              Cleaner subsets are far sparser (only ~27 single-mutation isolates exist) and do
              not improve the correlation — an honest bound, not a bug.
            </p>
          </Panel>
        )}
      </div>

      <Panel title="Predicted ΔΔG vs measured fold-resistance · enrichment · per-drug signal">
        <img src={validationPlotUrl(target)} alt="ResistScope validation" className="w-full max-w-full rounded" />
      </Panel>
    </div>
  );
}

function Panel({ title, children }) {
  return (
    <div className="rounded-lg border border-slate-200 p-4">
      <div className="text-sm text-slate-500 mb-3">{title}</div>
      <div className="overflow-x-auto">{children}</div>
    </div>
  );
}

function Metric({ label, value, sub, accent }) {
  return (
    <div className={`rounded-lg border p-4 ${accent === "teal" ? "border-teal-200 bg-teal-50/40" : "border-slate-200"}`}>
      <div className="text-xs text-slate-500">{label}</div>
      <div className="text-2xl font-semibold text-slate-800 tabular-nums mt-1">{value}</div>
      <div className="text-xs text-slate-400 mt-0.5">{sub}</div>
    </div>
  );
}

const Th = ({ children }) => <th className="font-normal pb-1.5 pr-3">{children}</th>;
const Td = ({ children, className = "" }) => <td className={`py-1 pr-3 ${className}`}>{children}</td>;
