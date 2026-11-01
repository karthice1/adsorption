// Build Elsevier-style Word files from the Markdown sources.
//   node build_docx.js                                   -> CASSIA_manuscript.docx, CASSIA_highlights.docx
//   node build_docx.js plain IN.md OUT.docx              -> plain document (response letter, supplement)
// Footer line-number suppression is applied afterwards by postprocess_docx.py.
// Conventions: A4, 2.54 cm margins, Times New Roman 12 pt, double spacing,
// continuous line numbers, page numbers, title page, highlights/abstract on
// their own pages, three-rule tables (caption above), figures (caption below),
// Harvard references with hanging indent.
const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, ImageRun, Table, TableRow, TableCell,
  AlignmentType, BorderStyle, WidthType, Footer, PageNumber, LevelFormat,
  LineNumberRestartFormat, HeadingLevel, ShadingType, LineRuleType,
} = require("docx");

const DIR = __dirname;
const MODE = process.argv[2] === "plain" ? "plain" : "manuscript";
const IN = MODE === "plain" ? process.argv[3] : (process.env.CASSIA_IN || "CASSIA_manuscript.md");
const OUT = MODE === "plain" ? process.argv[4] : "CASSIA_manuscript.docx";
const md = fs.readFileSync(path.join(DIR, IN), "utf8").split("\n");
const FONT = "Times New Roman";
const PAGE_W = 11906, MARGIN = 1440, TEXT_W = PAGE_W - 2 * MARGIN; // A4, DXA
const DOUBLE = { line: 480, lineRule: LineRuleType.AUTO }, SINGLE = { line: 240, lineRule: LineRuleType.AUTO };

// ---------- inline markdown -> runs ----------
// Handles **bold**, *italic*, `code`, and X_sub subscripts (outside code).
function runs(text, base = {}) {
  const out = [];
  const re = /(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`)/g;
  let last = 0, m;
  const plain = (t0, fmt) => {
    // explicit _{sub} / ^{sup} (from latex2uni.py), then the implicit X_sub rule
    const parts = t0.split(/(_\{[^{}]*\}|\^\{[^{}]*\})/);
    parts.forEach((t) => {
      if (/^_\{.*\}$/.test(t)) { out.push(new TextRun({ text: t.slice(2, -1), subScript: true, ...fmt })); return; }
      if (/^\^\{.*\}$/.test(t)) { out.push(new TextRun({ text: t.slice(2, -1), superScript: true, ...fmt })); return; }
      plain1(t, fmt);
    });
  };
  const plain1 = (t, fmt) => {
    // subscripts: q_max, C_e,VI, q_app,VI, pH_pzc, E_a, R_s
    const sre = /([A-Za-zτσ)])_([A-Za-z0-9]+(?:,[A-Za-z]+)?)/g;
    let l = 0, s;
    while ((s = sre.exec(t))) {
      if (s.index + 1 > l) out.push(new TextRun({ text: t.slice(l, s.index + 1), ...fmt }));
      out.push(new TextRun({ text: s[2], subScript: true, ...fmt }));
      l = s.index + s[0].length;
    }
    if (l < t.length) out.push(new TextRun({ text: t.slice(l), ...fmt }));
  };
  while ((m = re.exec(text))) {
    if (m.index > last) plain(text.slice(last, m.index), base);
    const tok = m[0];
    if (tok.startsWith("**")) plain(tok.slice(2, -2), { ...base, bold: true });
    else if (tok.startsWith("`")) out.push(new TextRun({ text: tok.slice(1, -1), font: "Courier New", size: 20, ...base }));
    else plain(tok.slice(1, -1), { ...base, italics: true });
    last = m.index + tok.length;
  }
  if (last < text.length) plain(text.slice(last), base);
  return out;
}

const para = (text, opts = {}) => new Paragraph({ spacing: DOUBLE, children: runs(text, opts.run), ...opts.p });

// ---------- blocks ----------
function heading(level, text, breakBefore = false) {
  const isRefs = /^(References|Highlights|Abstract)$/.test(text);
  const fmt = level === 2 ? { bold: true } : level === 3 ? { bold: true, italics: true } : { italics: true };
  return new Paragraph({
    heading: [HeadingLevel.HEADING_1, HeadingLevel.HEADING_1, HeadingLevel.HEADING_2, HeadingLevel.HEADING_3][level - 1],
    spacing: { ...DOUBLE, before: level === 2 ? 240 : 120 },
    keepNext: true, pageBreakBefore: breakBefore,
    children: runs(text, { ...fmt, font: FONT, size: 24, color: "000000" }),
    ...(isRefs ? {} : {}),
  });
}

function table(lines) {
  const rows = lines.filter((l) => !/^\|\s*-/.test(l)).map((l) =>
    l.replace(/^\||\|$/g, "").split("|").map((c) => c.trim()));
  const n = rows[0].length;
  // width proportional to max content length, min share
  const lens = Array.from({ length: n }, (_, j) => Math.max(...rows.map((r) => (r[j] || "").length)));
  // clamp so short columns stay readable and long prose columns do not dominate
  const wt = lens.map((x) => Math.min(Math.max(x, 13), 40));
  const tot = wt.reduce((a, b) => a + b, 0);
  const widths = wt.map((x) => Math.floor((TEXT_W * x) / tot));
  widths[n - 1] += TEXT_W - widths.reduce((a, b) => a + b, 0);
  const none = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };
  const rule = { style: BorderStyle.SINGLE, size: 8, color: "000000" };
  return new Table({
    width: { size: TEXT_W, type: WidthType.DXA },
    columnWidths: widths,
    borders: { top: rule, bottom: rule, left: none, right: none, insideHorizontal: none, insideVertical: none },
    rows: rows.map((r, i) => new TableRow({
      tableHeader: i === 0,
      children: r.map((c, j) => new TableCell({
        width: { size: widths[j], type: WidthType.DXA },
        borders: i === 0 ? { bottom: rule, top: rule } : i === rows.length - 1 ? { bottom: rule } : {},
        margins: { top: 40, bottom: 40, left: 80, right: 80 },
        children: [new Paragraph({ spacing: SINGLE, children: runs(c, { size: 20, bold: i === 0 }) })],
      })),
    })),
  });
}

function figure(src) {
  const file = path.join(DIR, src);
  const buf = fs.readFileSync(file);
  const w = buf.readUInt32BE(16), h = buf.readUInt32BE(20);
  const maxW = (w / h < 1.6 ? 4.3 : 6.2) * 96; // px at 96 dpi; near-square plots narrower
  const width = Math.round(Math.min(maxW, w)), height = Math.round((width * h) / w);
  return new Paragraph({
    alignment: AlignmentType.CENTER, spacing: { ...SINGLE, before: 240, after: 120 }, keepNext: true,
    children: [new ImageRun({ type: "png", data: buf, transformation: { width, height },
      altText: { title: path.basename(src), description: path.basename(src), name: path.basename(src) } })],
  });
}

// ---------- title page ----------
const titleLine = md.find((l) => l.startsWith("# ")).slice(2);
const titlePage = [
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { ...DOUBLE, after: 240 },
    children: [new TextRun({ text: titleLine, bold: true, size: 32 })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: DOUBLE, children: [
    new TextRun("Karthik [Surname]"), new TextRun({ text: "a,", superScript: true }), new TextRun({ text: "*", superScript: true })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { ...DOUBLE, after: 360 }, children: [
    new TextRun({ text: "a", superScript: true }),
    new TextRun({ text: " [Department, Institution, Street address, City, Postcode, Country]", italics: true })] }),
  para("*Corresponding author.* E-mail address: karthikresearchworks@gmail.com"),
  para("ORCID: [to be added]"),
  para("Article type: Research article"),
  para("Declarations of interest: none."),
];

// ---------- body ----------
const body = [];
let i = MODE === "plain" ? md.findIndex((l) => l.startsWith("# ")) + 1 : md.findIndex((l) => l.startsWith("## Highlights"));
const highlights = [];
let inHighlights = false, inRefs = false;
const pageBreakBefore = MODE === "plain" ? new Set() : new Set(["Highlights", "Abstract", "1. Introduction", "References"]);
let afterTable = false;

while (i < md.length) {
  const line = md[i];
  if (/^\s*$/.test(line) || /^---\s*$/.test(line)) { i++; continue; }
  const h = line.match(/^(#{2,4})\s+(.*)$/);
  if (h) {
    const text = h[2].trim();
    inHighlights = text === "Highlights";
    inRefs = text === "References";
    body.push(heading(h[1].length, text, pageBreakBefore.has(text)));
    afterTable = false;
    i++; continue;
  }
  if (line.startsWith("|")) {
    const block = [];
    while (i < md.length && md[i].startsWith("|")) block.push(md[i++]);
    body.push(table(block));
    afterTable = true;
    continue;
  }
  if (line.startsWith("```")) {
    i++;
    while (i < md.length && !md[i].startsWith("```")) {
      body.push(new Paragraph({ spacing: SINGLE, indent: { left: 567 }, keepLines: true,
        children: [new TextRun({ text: md[i], font: "Courier New", size: 18 })] }));
      i++;
    }
    i++; continue;
  }
  if (line.startsWith("<center>")) {
    body.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { ...DOUBLE, before: 120, after: 120 }, children: runs(line.slice(8), { italics: true }) }));
    i++; continue;
  }
  const img = line.match(/^!\[[^\]]*\]\(([^)]+)\)/);
  if (img) { body.push(figure(img[1])); i++; continue; }
  const bullet = line.match(/^(\s*)- (.*)$/);
  if (bullet) {
    const lvl = bullet[1].length >= 2 ? 1 : 0;
    if (inHighlights) highlights.push(bullet[2]);
    body.push(new Paragraph({ numbering: { reference: "bullets", level: lvl }, spacing: DOUBLE, children: runs(bullet[2]) }));
    i++; continue;
  }
  const num = line.match(/^(\d+)\.\s+(.*)$/);
  if (num) {
    body.push(new Paragraph({ spacing: DOUBLE, indent: { left: 567, hanging: 567 },
      children: [new TextRun(`${num[1]}.\t`), ...runs(num[2])], tabStops: [{ type: "left", position: 567 }] }));
    i++; continue;
  }
  if (inRefs) {
    body.push(new Paragraph({ spacing: DOUBLE, indent: { left: 567, hanging: 567 }, children: runs(line) }));
    i++; continue;
  }
  // figure / table captions: keep with neighbouring object
  const isTableCap = /^\*\*Table \d+\.\*\*/.test(line);
  const isFigCap = /^\*\*Fig\. \d+\.\*\*/.test(line);
  const sp = isFigCap ? { ...SINGLE, after: 240 } : isTableCap ? { ...SINGLE, after: 120 } : { ...DOUBLE };
  if (afterTable) { sp.before = 200; afterTable = false; }
  body.push(new Paragraph({
    spacing: sp,
    keepNext: isTableCap, alignment: AlignmentType.JUSTIFIED,
    children: runs(line, isFigCap || isTableCap ? { size: 22 } : {}),
  }));
  i++;
}

const common = {
  styles: {
    default: { document: { run: { font: FONT, size: 24 }, paragraph: { spacing: DOUBLE } } },
    paragraphStyles: [1, 2, 3].map((n) => ({
      id: `Heading${n}`, name: `Heading ${n}`, basedOn: "Normal", next: "Normal", quickFormat: true,
      run: { font: FONT, size: 24, bold: n < 3, italics: n > 1, color: "000000" },
      paragraph: { spacing: { before: 240, after: 120 }, outlineLevel: n - 1 },
    })),
  },
  numbering: { config: [{ reference: "bullets", levels: [
    { level: 0, format: LevelFormat.BULLET, text: "\u2022", alignment: AlignmentType.LEFT,
      style: { paragraph: { indent: { left: 567, hanging: 340 } } } },
    { level: 1, format: LevelFormat.BULLET, text: "\u2013", alignment: AlignmentType.LEFT,
      style: { paragraph: { indent: { left: 1134, hanging: 340 } } } }] }] },
};
const pageProps = {
  page: { size: { width: PAGE_W, height: 16838 }, margin: { top: MARGIN, bottom: MARGIN, left: MARGIN, right: MARGIN } },
};
const footer = { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER,
  children: [new TextRun({ children: [PageNumber.CURRENT] })] })] }) };

const plainTitle = [new Paragraph({ spacing: { ...DOUBLE, after: 240 },
  children: [new TextRun({ text: titleLine, bold: true, size: 28 })] })];
const doc = new Document({
  ...common,
  creator: "Karthik", title: titleLine,
  description: "Prepared to Elsevier Guide for Authors conventions (review format)",
  sections: [{
    properties: MODE === "plain" ? pageProps
      : { ...pageProps, lineNumbers: { countBy: 1, restart: LineNumberRestartFormat.CONTINUOUS, distance: 360 } },
    footers: footer,
    children: MODE === "plain" ? [...plainTitle, ...body] : [...titlePage, ...body],
  }],
});

const hlDoc = new Document({
  ...common,
  sections: [{ properties: pageProps, children: [
    new Paragraph({ spacing: DOUBLE, children: [new TextRun({ text: "Highlights", bold: true })] }),
    ...highlights.map((h) => new Paragraph({ numbering: { reference: "bullets", level: 0 }, spacing: DOUBLE, children: runs(h) })),
  ] }],
});

(async () => {
  fs.writeFileSync(path.join(DIR, OUT), await Packer.toBuffer(doc));
  if (MODE === "manuscript") fs.writeFileSync(path.join(DIR, "CASSIA_highlights.docx"), await Packer.toBuffer(hlDoc));
  console.log("written", OUT, "highlights:", highlights.length);
})();
