import { useEffect, useRef, useState } from "react";
import { Stage } from "ngl";

// Match the severity palette used across the app (red / amber / emerald / slate).
const SEV_COLOR = { high: "#dc2626", medium: "#d97706", low: "#059669", unknown: "#64748b" };

/**
 * 3D structure of the resistance mutation in context: the wildtype receptor
 * (cartoon), the co-crystal inhibitor filling the binding pocket (ball+stick),
 * and the mutated residue highlighted + colored by ΔΔG severity. This visually
 * grounds the numbers — you can see whether the mutation directly lines the
 * pocket (near the ligand) or acts allosterically.
 */
export default function StructureViewer({ mutation, position, severity = "unknown", target = "HIV1_PR", height = 300 }) {
  const containerRef = useRef(null);
  const [status, setStatus] = useState("loading");
  const [msg, setMsg] = useState("");

  useEffect(() => {
    const el = containerRef.current;
    if (!el) return;
    let disposed = false;
    const color = SEV_COLOR[severity] || SEV_COLOR.unknown;
    const sel = String(position); // resno in every kept chain (homodimer → both)
    const stage = new Stage(el, { backgroundColor: "white" });

    // NGL grabs the canvas size at construction; inside a freshly-expanded grid
    // cell that can be 0 → a blank canvas. Re-sync the WebGL viewport whenever
    // the container's box changes (and once, right after layout settles).
    const resize = () => { try { stage.handleResize(); } catch { /* pre-init */ } };
    const ro = new ResizeObserver(resize);
    ro.observe(el);

    async function loadPDB(url, name) {
      const r = await fetch(url);
      if (!r.ok) throw new Error(`${name}: HTTP ${r.status}`);
      const text = await r.text();
      return stage.loadFile(new Blob([text], { type: "text/plain" }), { ext: "pdb", name });
    }

    (async () => {
      try {
        const receptor = await loadPDB(`/api/structure/receptor?target=${target}`, "receptor");
        if (disposed) return;
        receptor.addRepresentation("cartoon", { color: "#c3cbd6", opacity: 0.85, side: "front" });
        receptor.addRepresentation("ball+stick", { sele: sel, color, aspectRatio: 1.5, radiusScale: 1.4 });
        receptor.addRepresentation("spacefill", { sele: sel, color, opacity: 0.18 });

        try {
          const lig = await loadPDB(`/api/structure/ligand?target=${target}`, "ligand");
          if (!disposed) lig.addRepresentation("ball+stick", { color: "element", multipleBond: "symmetric" });
        } catch { /* co-crystal ligand is optional context */ }

        if (disposed) return;
        resize();
        stage.autoView(300);
        setStatus("ready");
        // One deferred kick after the expand animation / fonts settle — this is
        // what reliably fills the canvas when the row was 0-sized at mount.
        setTimeout(() => { if (!disposed) { resize(); stage.autoView(0); } }, 120);
      } catch (e) {
        if (!disposed) { setStatus("error"); setMsg(String(e.message || e)); }
      }
    })();

    return () => {
      disposed = true;
      ro.disconnect();
      try { stage.dispose(); } catch { /* already gone */ }
    };
  }, [mutation, position, severity, target]);

  return (
    <div className="relative rounded border border-slate-200 overflow-hidden bg-white" style={{ height }}>
      <div ref={containerRef} className="w-full h-full" />
      {status !== "ready" && (
        <div className="absolute inset-0 flex items-center justify-center text-xs text-slate-400 pointer-events-none px-4 text-center">
          {status === "error" ? `structure unavailable — ${msg}` : "loading 3D structure…"}
        </div>
      )}
      <div className="absolute top-1.5 left-2 text-[10px] text-slate-400 pointer-events-none">
        <span className="inline-block w-2 h-2 rounded-full align-middle mr-1" style={{ background: SEV_COLOR[severity] || SEV_COLOR.unknown }} />
        mutated residue · inhibitor in pocket
      </div>
      <div className="absolute bottom-1 right-2 text-[10px] text-slate-400 pointer-events-none">
        drag · scroll to zoom
      </div>
    </div>
  );
}
