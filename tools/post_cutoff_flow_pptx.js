// Editable (PowerPoint) version of the post-cutoff approval-set flow chart, in Nature artwork style:
// 183 mm double-column width, Arial, 5 to 7 pt text, 8 pt bold lowercase panel labels, 0.5 pt lines,
// greyscale with one muted accent. Every number comes from flow_counts.json, which
// tools/make_post_cutoff_flow_pptx.py re-derives from the shipped tables before calling this script.
//
//   node tools/post_cutoff_flow_pptx.js <counts.json> <out.pptx>
const fs = require("fs");
const pptxgen = require("pptxgenjs");

const [, , countsPath, outPath] = process.argv;
const c = JSON.parse(fs.readFileSync(countsPath, "utf8"));

const W = 7.2, H = 8.45;                       // 183 x 215 mm
const INK = "000000", GREY = "595959", RULE = "7F7F7F", FILL = "F2F2F2", ACCENT_FILL = "E4ECEA", ACCENT = "1F4E46", DROP = "8C1D18";
const FONT = "Arial";
const fmt = n => n.toLocaleString("en-US");

const pres = new pptxgen();
pres.defineLayout({ name: "NATURE_2COL", width: W, height: H });
pres.layout = "NATURE_2COL";
const s = pres.addSlide();

function box(x, y, w, h, runs, opts = {}) {
  s.addText(runs, {
    shape: pres.shapes.RECTANGLE, x, y, w, h,
    fill: { color: opts.fill || "FFFFFF" },
    line: { color: opts.lineColor || INK, width: opts.lineWidth || 0.5 },
    fontFace: FONT, fontSize: opts.fontSize || 6.5, color: INK,
    margin: opts.margin === undefined ? [3, 4, 3, 4] : opts.margin,
    align: opts.align || "left", valign: opts.valign || "middle",
    isTextBox: true, paraSpaceAfter: opts.paraSpaceAfter || 0,
  });
}
function label(x, y, w, h, runs, opts = {}) {
  s.addText(runs, { x, y, w, h, fontFace: FONT, fontSize: opts.fontSize || 6, color: opts.color || GREY,
    margin: 0, align: opts.align || "left", valign: opts.valign || "top", isTextBox: true });
}
function vline(x, y0, y1, arrow = true) {
  s.addShape(pres.shapes.LINE, { x, y: y0, w: 0, h: y1 - y0, line: { color: INK, width: 0.5, endArrowType: arrow ? "triangle" : undefined } });
}
function hline(x0, x1, y, arrow = true) {
  s.addShape(pres.shapes.LINE, { x: x0, y, w: x1 - x0, h: 0, line: { color: INK, width: 0.5, endArrowType: arrow ? "triangle" : undefined } });
}
const B = t => ({ text: t, options: { bold: true } });
const N = t => ({ text: t });
const line = (parts) => parts.map((p, i) => i === parts.length - 1 ? { ...p, options: { ...(p.options || {}), breakLine: true } } : p);

// ---------------------------------------------------------------- panel a: regulator sources
label(0.0, 0.0, 0.3, 0.16, [B("a")], { fontSize: 8, color: INK });
const REG = [
  ["FDA", "FDA (CDER)", "efficacy supplements, sub-type New Indication"],
  ["EMA", "EMA (CHMP)", "positive opinions on extensions of indication"],
  ["HC", "Health Canada", "SNDS whose reason text names an indication"],
  ["PMDA", "PMDA", "approvals whose notes say new indication"],
  ["TGA", "TGA", "extension-of-indications registrations"],
];
const rw = 1.32, rg = (W - 5 * rw) / 4, ry = 0.17, rh = 0.66;
REG.forEach(([k, name, what], i) => {
  const x = i * (rw + rg);
  box(x, ry, rw, rh, [
    { text: `n = ${fmt(c.rows[k])}`, options: { bold: true, breakLine: true } },
    { text: name, options: { breakLine: true } },
    { text: what, options: { fontSize: 5.5, color: GREY } },
  ], { fill: FILL, valign: "top", margin: [3, 4, 2, 4] });
  vline(x + rw / 2, ry + rh, 0.98);
});

// ---------------------------------------------------------------- spine
const SX = 0.0, SW = 4.05, EX = 4.45, EW = W - EX;
// S0
box(SX, 0.98, W, 0.44, [
  B(`n = ${fmt(c.n0)}`), N("  candidate rows, one per regulator decision (S0)"), { text: "", options: { breakLine: true } },
  { text: "2023-01-01 to 2026-09-17; every row cites the regulator's document", options: { fontSize: 5.5, color: GREY } },
], { fill: FILL });
vline(SX + 0.5, 1.42, 1.60);
// annotation step
box(SX, 1.60, W, 0.74, [
  B("Indication text, identifiers, review"), { text: "", options: { breakLine: true } },
  { text: `Indication text from the regulator record: FDA approval letter, “provides for” sentence (${c.text.FDA} of ${c.rows.FDA}); EMA EPAR procedural-steps scope matched to the opinion date (${c.text.EMA} of ${c.rows.EMA}); Health Canada, PMDA and TGA text as published. Drug to DrugBank and disease to MONDO with the NCATS Name Resolver 1.5.2 and Node Normalizer 2.4.1 (responses cached) and the MONDO SSSOM. Class, disease and drug reviewed in three passes: 2,258 recorded decisions, each with a rationale.`, options: { fontSize: 5.5, color: GREY } },
], { valign: "top" });
// S0 -> S1 exclusions
const ex1y = 2.52, ex1h = 1.62;
vline(SX + 0.5, 2.34, 4.32);
hline(SX + 0.5, EX, ex1y + 0.35);
const nn = c.not_new, nc = c.not_change;
const notNew = nn.line_extension + nn.expanded_population + nn.ingredient_prior_indication + nn.formulation_or_route + nn.biosimilar_alignment + nn.first_indication;
const notChange = nc.wording_only + nc.unclear + nc.out_of_scope_veterinary + nc.commitment_fulfilment + nc.excluded + nc.indication_removed;
const ex1total = c.no_text + notNew + notChange + c.no_target + c.pending;
box(EX, ex1y, EW, ex1h, [
  ...line([B(`Excluded before S1, n = ${fmt(ex1total)}`)]),
  ...line([B(`Not a new indication for the ingredient, n = ${notNew}`)]),
  ...line([{ text: `line extension ${nn.line_extension}; population extension ${nn.expanded_population}; indication already held by the ingredient in another product or approved before 2023 elsewhere ${nn.ingredient_prior_indication}; new formulation ${nn.formulation_or_route}; biosimilar ${nn.biosimilar_alignment}; first approval of a new substance ${nn.first_indication}`, options: { fontSize: 5.5, color: GREY } }]),
  ...line([B(`No indication text in the source, n = ${c.no_text}`)]),
  ...line([{ text: "129 CHMP cards without a matched EPAR entry, 10 FDA supplements without a letter, NOC reasons that say only “New indication”", options: { fontSize: 5.5, color: GREY } }]),
  ...line([B(`Not an indication event, n = ${notChange}`)]),
  ...line([{ text: `wording only ${nc.wording_only}; unclear ${nc.unclear}; veterinary ${nc.out_of_scope_veterinary}; commitment fulfilment ${nc.commitment_fulfilment}; withdrawn or excluded ${nc.excluded + nc.indication_removed}`, options: { fontSize: 5.5, color: GREY } }]),
  ...line([B(`No disease target, n = ${c.no_target}`), { text: "  diagnostic, procedural, cosmetic, contraceptive", options: { fontSize: 5.5, color: GREY } }]),
  B(`Identifier pending a live lookup, n = ${c.pending}`),
], { valign: "top", paraSpaceAfter: 1.5, lineColor: RULE });
// S1
box(SX, 4.32, SW, 0.44, [
  B(`n = ${fmt(c.n1)}`), N("  reviewed new-indication rows, both identifiers resolved (S1)"),
], { fill: FILL });
// S1 -> S2
vline(SX + 0.5, 4.76, 5.36);
hline(SX + 0.5, EX, 5.08);
const s2 = c.s2, s2total = Object.values(s2).reduce((a, b) => a + b, 0);
box(EX, 4.84, EW, 0.5, [
  ...line([B(`Not mappable to a graph node, n = ${s2total}`)]),
  { text: `vaccines ${s2.vaccine || 0}; no DrugBank accession ${s2["no DrugBank accession"] || 0}; tumour-agnostic approvals whose target is the root “cancer” term ${s2["generic disease target"] || 0}; diagnostic agents ${s2["diagnostic agent"] || 0}`, options: { fontSize: 5.5, color: GREY } },
], { valign: "top", lineColor: RULE });
// S2
box(SX, 5.36, SW, 0.44, [
  B(`n = ${fmt(c.n2)}`), N("  KG-mappable rows (S2): single agent with a DrugBank accession, specific MONDO term"),
], { fill: FILL });
// S2 -> S3
vline(SX + 0.5, 5.80, 6.42);
hline(SX + 0.5, EX, 6.12);
box(EX, 5.86, EW, 0.52, [
  ...line([B(`Collapsed across regulators, ${c.events - c.pairs} events`)]),
  { text: `${c.n2} rows name ${c.events} (row, disease) events; one (DrugBank id, MONDO id) pair each, earliest decision date kept, every approving regulator listed; ${c.multi} pairs approved by two or more regulators`, options: { fontSize: 5.5, color: GREY } },
], { valign: "top", lineColor: RULE });
// S3
box(SX, 6.42, SW, 0.52, [
  B(`n = ${c.pairs}`), N(`  unique drug-disease pairs (S3): ${c.drugs} drugs, ${c.diseases} diseases`), { text: "", options: { breakLine: true } },
  { text: "the post-cutoff approval set; per-graph checks in panel c", options: { fontSize: 5.5, color: GREY } },
], { fill: ACCENT_FILL, lineWidth: 1, lineColor: ACCENT });

// ---------------------------------------------------------------- panel b: composition by regulator (native table)
label(0.0, 7.12, 0.3, 0.16, [B("b")], { fontSize: 8, color: INK });
const order = ["FDA", "EMA", "HC", "PMDA", "TGA"];
const names = { FDA: "FDA", EMA: "EMA", HC: "Health Canada", PMDA: "PMDA", TGA: "TGA" };
const hdr = (t, a) => ({ text: t, options: { bold: true, fontSize: 6, fontFace: FONT, color: INK, align: a || "right", border: [{ type: "none" }, { type: "none" }, { type: "solid", color: INK, pt: 0.5 }, { type: "none" }], margin: [2, 3, 2, 3] } });
const cell = (t, a) => ({ text: String(t), options: { fontSize: 6, fontFace: FONT, color: INK, align: a || "right", border: [{ type: "none" }, { type: "none" }, { type: "none" }, { type: "none" }], margin: [1.5, 3, 1.5, 3] } });
const rows = [[hdr("Regulator", "left"), hdr("Rows in"), hdr("Rows in a pair"), hdr("Pairs approved"), hdr("Pairs, first approver")]];
order.forEach(k => rows.push([cell(names[k], "left"), cell(fmt(c.rows[k])), cell(c.rows_by_reg_s3[k] || 0), cell(c.any[k] || 0), cell(c.first[k] || 0)]));
rows.push([cell("All", "left"), cell(fmt(c.n0)), cell(c.rows_s3), cell(c.pairs), cell(c.pairs)]);
rows[rows.length - 1].forEach(x => { x.options.border[0] = { type: "solid", color: INK, pt: 0.5 }; x.options.bold = true; });
s.addTable(rows, { x: 0.15, y: 7.12, w: SW - 0.15, colW: [0.95, 0.6, 0.75, 0.7, 0.9], rowH: 0.13, fontFace: FONT, fontSize: 6 });

// ---------------------------------------------------------------- panel c: per-graph checks (native table)
label(EX - 0.25, 7.12, 0.3, 0.16, [B("c")], { fontSize: 8, color: INK });
const kgs = Object.keys(c.kg);
const rows2 = [[hdr(`Check, applied to the ${c.pairs} pairs`, "left"), ...kgs.map(k => hdr(k === "biokg" ? "BioKG" : k === "hetionet" ? "Hetionet" : k))]];
const checkRows = [["C2 covered (both nodes, drug degree ≥ 5)", "C2"], ["C1 prospective (after snapshot, no direct edge)", "C1"], ["C3 repurposing (> 2 is_a steps from prior)", "C3"]];
checkRows.forEach(([t, st]) => rows2.push([cell(t, "left"), ...kgs.map(k => cell(c.kg[k][st]))]));
rows2.push([cell("Common support (evaluable in every graph)", "left"), cell(c.common.C1), cell("")]);
s.addTable(rows2, { x: EX - 0.1, y: 7.12, w: EW + 0.1, colW: [1.9, 0.5, 0.45], rowH: 0.13, fontFace: FONT, fontSize: 6 });

label(0.15, 8.08, SW, 0.12, [N(`${c.multi} of the ${c.pairs} pairs were approved by two or more regulators.`)], { fontSize: 5.5 });
label(0.0, H - 0.14, W, 0.12, [N(`repurposing-dataset ${c.version}; every count re-derived from the shipped tables by tools/make_post_cutoff_flow_pptx.py`)], { fontSize: 5, color: GREY });

pres.writeFile({ fileName: outPath }).then(() => console.log("wrote", outPath));
