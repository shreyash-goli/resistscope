import { useState } from "react";
import { intakeTarget, assembleTarget, saveTarget } from "../api";

/**
 * Bring-your-own-target flow: give a PDB (id or uploaded file) + a protein name,
 * the API parses the structure (chains + candidate pocket ligands) and a Claude
 * agent proposes the inhibitors + resistance mutations. The user confirms the
 * pocket and edits the mutation list, then saves — the target appears in the
 * selector (triage-only until its receptors are built on the GPU worker).
 */
export default function AddTargetModal({ onClose, onSaved }) {
  const [step, setStep] = useState("input"); // input | review
  const [busy, setBusy] = useState("");
  const [error, setError] = useState(null);

  const [pdbId, setPdbId] = useState("");
  const [pdbText, setPdbText] = useState(null);
  const [protein, setProtein] = useState("");

  const [intake, setIntake] = useState(null);
  const [pocketIdx, setPocketIdx] = useState(0);
  const [spec, setSpec] = useState(null);
  const [citations, setCitations] = useState([]);

  const [targetId, setTargetId] = useState("");
  const [label, setLabel] = useState("");
  const [drugsText, setDrugsText] = useState("");
  const [mutText, setMutText] = useState("");

  function onFile(e) {
    const f = e.target.files?.[0];
    if (!f) return;
    const reader = new FileReader();
    reader.onload = () => { setPdbText(String(reader.result)); setPdbId(f.name.replace(/\.(pdb|cif)$/i, "")); };
    reader.readAsText(f);
  }

  async function analyze() {
    setError(null);
    if (!protein.trim() || (!pdbId.trim() && !pdbText)) {
      setError("Enter a PDB id (or upload a .pdb) and a protein name.");
      return;
    }
    try {
      setBusy("Parsing structure…");
      const ik = (await intakeTarget(pdbText ? { pdb_text: pdbText } : { pdb_id: pdbId.trim() })).intake;
      setIntake(ik);
      setPocketIdx(0);

      setBusy("Claude is researching drugs + resistance mutations (PubMed)…");
      const a = await assembleTarget({ protein: protein.trim(), hint: ik.title || "" });
      setSpec(a.spec);
      setCitations(a.citations || []);
      setTargetId((a.spec.target_id || pdbId || "USER_TARGET").toUpperCase());
      setLabel(a.spec.label || protein.trim());
      setDrugsText((a.spec.drugs || []).map((d) => `${d.abbrev} = ${d.name}`).join("\n"));
      setMutText((a.spec.resistance_mutations || []).map((m) => m.mutation).join(", "));
      setStep("review");
    } catch (e) {
      setError(e.status === 502 ? "The research agent failed — try a more specific protein name." : String(e));
    } finally {
      setBusy("");
    }
  }

  async function save() {
    setError(null);
    const lig = intake?.ligands?.[pocketIdx];
    const drugs = drugsText.split("\n").map((l) => {
      const [ab, ...rest] = l.split(/[=:]/);
      return ab.trim() ? { abbrev: ab.trim(), name: (rest.join("=").trim() || ab.trim()) } : null;
    }).filter(Boolean);
    const mutations = mutText.split(/[,\s]+/).map((s) => s.trim().toUpperCase()).filter(Boolean);
    if (!drugs.length || !mutations.length) {
      setError("Need at least one drug and one mutation.");
      return;
    }
    try {
      setBusy("Saving target…");
      const body = {
        target_id: targetId.trim(),
        label: label.trim() || targetId.trim(),
        pdb_id: pdbText ? null : pdbId.trim(),
        pdb_text: pdbText || null,
        ligand_hetcode: lig?.hetcode || null,
        docking_center: lig?.centroid || intake?.suggested?.docking_center || [0, 0, 0],
        chains: intake?.suggested?.chains || ["A"],
        mutate_chains: intake?.suggested?.mutate_chains || ["A"],
        drugs,
        mutations,
        binding_site: spec?.binding_site || "",
        validatable: !!spec?.dataset?.exists,
        citations,
      };
      const res = await saveTarget(body);
      onSaved(res.saved);
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy("");
    }
  }

  return (
    <div className="fixed inset-0 bg-black/30 flex items-start justify-center z-50 p-4 overflow-y-auto" onClick={onClose}>
      <div className="bg-white rounded-lg shadow-xl w-full max-w-2xl my-8" onClick={(e) => e.stopPropagation()}>
        <div className="flex items-center justify-between px-5 py-3 border-b border-slate-200">
          <h2 className="text-base font-semibold text-slate-800">Add a target</h2>
          <button onClick={onClose} className="text-slate-400 hover:text-slate-600 text-lg leading-none">×</button>
        </div>

        <div className="p-5 space-y-4">
          {error && <div className="rounded border border-red-200 bg-red-50 text-red-700 text-sm px-3 py-2">{error}</div>}

          {step === "input" && (
            <>
              <p className="text-sm text-slate-500">
                Give a receptor structure and a protein name. We parse the pocket from the structure,
                and a Claude agent researches its inhibitors + resistance mutations from the literature.
              </p>
              <label className="block text-sm">
                <span className="text-slate-500">PDB id (with an inhibitor bound)</span>
                <input value={pdbId} onChange={(e) => { setPdbId(e.target.value); setPdbText(null); }}
                  placeholder="e.g. 2HU4" className="mt-1 w-full border border-slate-300 rounded px-2 py-1.5 font-mono text-sm" />
              </label>
              <div className="text-xs text-slate-400">or upload a .pdb / .cif file: <input type="file" accept=".pdb,.cif,.mmcif" onChange={onFile} className="text-xs" /></div>
              <label className="block text-sm">
                <span className="text-slate-500">Protein / target name</span>
                <input value={protein} onChange={(e) => setProtein(e.target.value)}
                  placeholder="e.g. influenza A neuraminidase" className="mt-1 w-full border border-slate-300 rounded px-2 py-1.5 text-sm" />
              </label>
              <button onClick={analyze} disabled={!!busy}
                className="px-4 py-1.5 rounded bg-teal-600 text-white text-sm font-medium disabled:bg-slate-300 hover:bg-teal-700">
                {busy || "Analyze structure + research"}
              </button>
            </>
          )}

          {step === "review" && intake && spec && (
            <>
              <div className="grid grid-cols-2 gap-3">
                <label className="text-sm"><span className="text-slate-500">Target id</span>
                  <input value={targetId} onChange={(e) => setTargetId(e.target.value.toUpperCase())} className="mt-1 w-full border border-slate-300 rounded px-2 py-1 font-mono text-sm" /></label>
                <label className="text-sm"><span className="text-slate-500">Label</span>
                  <input value={label} onChange={(e) => setLabel(e.target.value)} className="mt-1 w-full border border-slate-300 rounded px-2 py-1 text-sm" /></label>
              </div>

              <label className="block text-sm">
                <span className="text-slate-500">Pocket ligand (confirm this is the inhibitor site, not a glycan/ion)</span>
                <select value={pocketIdx} onChange={(e) => setPocketIdx(Number(e.target.value))}
                  className="mt-1 w-full border border-slate-300 rounded px-2 py-1.5 text-sm font-mono">
                  {intake.ligands.map((l, i) => (
                    <option key={i} value={i}>{l.hetcode} · chain {l.chain} · {l.n_atoms} atoms · center ({l.centroid.map((v) => v.toFixed(1)).join(", ")})</option>
                  ))}
                </select>
                <span className="text-[11px] text-slate-400">Structure: {intake.title || "—"} · chains {intake.suggested.chains.join("/")}</span>
              </label>

              <label className="block text-sm">
                <span className="text-slate-500">Inhibitors (one <span className="font-mono">ABBR = name</span> per line)</span>
                <textarea value={drugsText} onChange={(e) => setDrugsText(e.target.value)} rows={3}
                  className="mt-1 w-full border border-slate-300 rounded px-2 py-1.5 font-mono text-xs" />
              </label>

              <label className="block text-sm">
                <span className="text-slate-500">Resistance mutations (comma-separated, e.g. H275Y, E119V)</span>
                <textarea value={mutText} onChange={(e) => setMutText(e.target.value)} rows={2}
                  className="mt-1 w-full border border-slate-300 rounded px-2 py-1.5 font-mono text-xs" />
              </label>

              <div className="text-xs text-slate-500">
                Dataset: {spec.dataset?.exists
                  ? <span className="text-emerald-600">a public resistance dataset exists ({spec.dataset.name || "found"}) — validatable</span>
                  : <span>none found — this target is triage-only</span>}
                {citations.length > 0 && <> · {citations.length} citations</>}
              </div>

              <div className="rounded border border-amber-200 bg-amber-50 text-amber-800 text-xs px-3 py-2">
                On save this target joins the selector immediately. Live ΔΔG docking needs its mutant
                receptors built on the GPU worker (structure prep + scripts/03) — until then it shows
                "docking pending", exactly like the RT target.
              </div>

              <div className="flex gap-2">
                <button onClick={save} disabled={!!busy}
                  className="px-4 py-1.5 rounded bg-teal-600 text-white text-sm font-medium disabled:bg-slate-300 hover:bg-teal-700">
                  {busy || "Save target"}
                </button>
                <button onClick={() => setStep("input")} disabled={!!busy} className="px-3 py-1.5 rounded text-sm text-slate-500 hover:text-slate-700">Back</button>
              </div>
            </>
          )}
        </div>
      </div>
    </div>
  );
}
