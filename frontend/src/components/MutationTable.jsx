import { useState, useMemo, Fragment, lazy, Suspense } from "react";

// NGL is ~1.3 MB; only load it when a user actually expands a mutation row.
const StructureViewer = lazy(() => import("./StructureViewer"));

const SEV = {
  high: "bg-red-100 text-red-700",
  medium: "bg-amber-100 text-amber-700",
  low: "bg-emerald-100 text-emerald-700",
  unknown: "bg-slate-100 text-slate-500",
};

export default function MutationTable({ mutations, target = "HIV1_PR" }) {
  const [sortKey, setSortKey] = useState("delta_delta_g");
  const [asc, setAsc] = useState(false);
  const [open, setOpen] = useState(null);
  const [primaryOnly, setPrimaryOnly] = useState(false);

  const maxAbs = useMemo(
    () => Math.max(0.01, ...mutations.map((m) => Math.abs(m.delta_delta_g || 0))),
    [mutations]
  );

  const rows = useMemo(() => {
    let r = primaryOnly ? mutations.filter((m) => m.is_primary) : mutations.slice();
    r.sort((a, b) => {
      const av = a[sortKey], bv = b[sortKey];
      const c = typeof av === "string" ? av.localeCompare(bv) : (av ?? 0) - (bv ?? 0);
      return asc ? c : -c;
    });
    return r;
  }, [mutations, sortKey, asc, primaryOnly]);

  const toggle = (k) => (k === sortKey ? setAsc(!asc) : (setSortKey(k), setAsc(false)));
  const arrow = (k) => (k === sortKey ? (asc ? " ↑" : " ↓") : "");

  return (
    <div className="rounded-lg border border-slate-200 overflow-hidden">
      <div className="flex items-center justify-between px-4 py-2 bg-slate-50 border-b border-slate-200">
        <span className="text-sm text-slate-500">{rows.length} mutations</span>
        <label className="text-sm text-slate-600 flex items-center gap-2 cursor-pointer">
          <input type="checkbox" checked={primaryOnly} onChange={(e) => setPrimaryOnly(e.target.checked)} />
          major DRMs only
        </label>
      </div>
      <table className="w-full text-sm">
        <thead>
          <tr className="text-left text-slate-500 border-b border-slate-200">
            <Th onClick={() => toggle("mutation")}>Mutation{arrow("mutation")}</Th>
            <Th onClick={() => toggle("delta_g")} className="text-right">ΔG{arrow("delta_g")}</Th>
            <Th onClick={() => toggle("delta_delta_g")} className="w-1/2">ΔΔG (resistance){arrow("delta_delta_g")}</Th>
            <Th onClick={() => toggle("severity")}>Severity{arrow("severity")}</Th>
          </tr>
        </thead>
        <tbody>
          {rows.map((m) => {
            const isOpen = open === m.mutation;
            const w = (Math.abs(m.delta_delta_g || 0) / maxAbs) * 50; // % of half-width
            const pos = (m.delta_delta_g || 0) >= 0;
            return (
              <Fragment key={m.mutation}>
                <tr
                  onClick={() => setOpen(isOpen ? null : m.mutation)}
                  className="border-b border-slate-100 hover:bg-slate-50 cursor-pointer"
                >
                  <td className="px-4 py-2 font-mono font-medium text-slate-800">
                    {m.mutation}
                    {m.is_primary && (
                      <span className="ml-2 text-[10px] uppercase tracking-wide text-teal-600 bg-teal-50 rounded px-1 py-0.5">
                        major
                      </span>
                    )}
                  </td>
                  <td className="px-4 py-2 text-right font-mono tabular-nums text-slate-600">
                    {m.delta_g == null ? "—" : m.delta_g.toFixed(2)}
                  </td>
                  <td className="px-4 py-2">
                    <div className="relative h-4 flex items-center">
                      <div className="absolute left-1/2 top-0 bottom-0 w-px bg-slate-200" />
                      <div
                        className={`absolute h-3 rounded-sm ${pos ? "bg-red-400" : "bg-emerald-400"}`}
                        style={{ left: pos ? "50%" : `${50 - w}%`, width: `${w}%` }}
                      />
                      <span className="absolute right-0 font-mono tabular-nums text-xs text-slate-500">
                        {(m.delta_delta_g >= 0 ? "+" : "") + m.delta_delta_g.toFixed(2)}
                      </span>
                    </div>
                  </td>
                  <td className="px-4 py-2">
                    <span className={`text-xs rounded px-2 py-0.5 ${SEV[m.severity] || SEV.unknown}`}>
                      {m.severity}
                    </span>
                  </td>
                </tr>
                {isOpen && (
                  <tr className="bg-slate-50">
                    <td colSpan={4} className="px-4 py-3">
                      <div className="grid grid-cols-1 lg:grid-cols-5 gap-4 items-start">
                        <div className="lg:col-span-2">
                          <Suspense
                            fallback={
                              <div className="h-[300px] rounded border border-slate-200 flex items-center justify-center text-xs text-slate-400">
                                loading 3D viewer…
                              </div>
                            }
                          >
                            <StructureViewer
                              mutation={m.mutation}
                              position={parseInt(m.mutation.slice(1, -1), 10)}
                              severity={m.severity}
                              target={target}
                            />
                          </Suspense>
                        </div>
                        <div className="lg:col-span-3 text-slate-600 text-sm leading-relaxed">
                          {m.explanation ? (
                            <>
                              <span className="text-xs uppercase tracking-wide text-slate-400">
                                Mechanistic hypothesis (Claude)
                              </span>
                              <p className="mt-1">{m.explanation}</p>
                              {m.citations && m.citations.length > 0 && (
                                <div className="mt-3">
                                  <span className="text-xs uppercase tracking-wide text-slate-400">
                                    Literature (PubMed)
                                  </span>
                                  <ul className="mt-1 space-y-0.5">
                                    {m.citations.map((c) => (
                                      <li key={c.pmid} className="text-xs">
                                        <a
                                          href={c.url}
                                          target="_blank"
                                          rel="noreferrer"
                                          className="text-teal-600 hover:underline"
                                        >
                                          {c.journal} {c.year}
                                        </a>
                                        <span className="text-slate-500"> — {c.title}</span>
                                      </li>
                                    ))}
                                  </ul>
                                </div>
                              )}
                            </>
                          ) : (
                            <span className="text-slate-400 italic">
                              No cached explanation for this mutation. The 3D view still shows where
                              this residue sits relative to the inhibitor pocket.
                            </span>
                          )}
                        </div>
                      </div>
                    </td>
                  </tr>
                )}
              </Fragment>
            );
          })}
        </tbody>
      </table>
    </div>
  );
}

function Th({ children, onClick, className = "" }) {
  return (
    <th onClick={onClick} className={`px-4 py-2 font-medium select-none cursor-pointer hover:text-slate-700 ${className}`}>
      {children}
    </th>
  );
}
