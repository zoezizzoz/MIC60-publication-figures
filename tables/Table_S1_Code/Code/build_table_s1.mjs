import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { Workbook, SpreadsheetFile } from "@oai/artifact-tool";

const scriptDir = path.dirname(fileURLToPath(import.meta.url));
const tableDir = path.resolve(scriptDir, "..");
const inputDir = process.env.TABLE_S1_INPUT_DIR
  ? path.resolve(process.env.TABLE_S1_INPUT_DIR)
  : path.join(tableDir, "Original_Data");
const femalePath = path.join(inputDir, "DEG_CSF_vs_WRF.csv");
const malePath = path.join(inputDir, "DEG_CSM_vs_WRM.csv");
const normalizedPath = path.join(inputDir, "normalized_counts.csv");
const outputDir = path.join(tableDir, "Rebuilt_Output");
const outputPath = path.join(outputDir, "Table_S1.xlsx");

function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"' && text[i + 1] === '"') {
        field += '"';
        i += 1;
      } else if (ch === '"') {
        quoted = false;
      } else {
        field += ch;
      }
    } else if (ch === '"') {
      quoted = true;
    } else if (ch === ",") {
      row.push(field);
      field = "";
    } else if (ch === "\n") {
      row.push(field.replace(/\r$/, ""));
      rows.push(row);
      row = [];
      field = "";
    } else {
      field += ch;
    }
  }
  if (field.length || row.length) {
    row.push(field.replace(/\r$/, ""));
    rows.push(row);
  }
  return rows;
}

function numeric(value) {
  if (value === undefined || value === null || value === "" || value === "NA" || value === "NaN") return null;
  const result = Number(value);
  return Number.isFinite(result) ? result : null;
}

function loadDeTable(text) {
  const rows = parseCsv(text);
  const header = rows[0];
  const idx = Object.fromEntries(header.map((name, i) => [name, i]));
  return rows.slice(1).filter((r) => r.length > 1).map((r) => [
    r[idx.gene],
    numeric(r[idx.baseMean]),
    numeric(r[idx.log2FoldChange]),
    numeric(r[idx.lfcSE]),
    numeric(r[idx.stat]),
    numeric(r[idx.pvalue]),
    numeric(r[idx.padj]),
  ]);
}

function loadNormalizedCounts(text) {
  const rows = parseCsv(text);
  const headers = ["gene", ...rows[0].slice(1)];
  const values = rows.slice(1).filter((r) => r.length > 1).map((r) => [r[0], ...r.slice(1).map(numeric)]);
  return { headers, values };
}

const [femaleText, maleText, normalizedText] = await Promise.all([
  fs.readFile(femalePath, "utf8"),
  fs.readFile(malePath, "utf8"),
  fs.readFile(normalizedPath, "utf8"),
]);
const femaleRows = loadDeTable(femaleText);
const maleRows = loadDeTable(maleText);
const normalized = loadNormalizedCounts(normalizedText);

if (femaleRows.length !== 17173 || maleRows.length !== 17173 || normalized.values.length !== 17173) {
  throw new Error(`Unexpected row count: female=${femaleRows.length}, male=${maleRows.length}, normalized=${normalized.values.length}`);
}

const workbook = Workbook.create();
const readme = workbook.worksheets.add("README");
const female = workbook.worksheets.add("Female DE");
const male = workbook.worksheets.add("Male DE");
const counts = workbook.worksheets.add("Normalized counts");

const navy = "#17365D";
const paleBlue = "#EEF5FB";
const gray = "#E7E6E6";
const green = "#E2F0D9";
const greenText = "#375623";
const red = "#FCE4D6";
const redText = "#9C0006";
const amber = "#FFF2CC";
const amberText = "#7F6000";
const line = "#D9E2F3";

readme.showGridLines = false;
readme.mergeCells("A1:H1");
readme.getRange("A1").values = [["Supplementary Table S1. RNA-seq data"]];
readme.getRange("A1:H1").format = {
  fill: navy,
  font: { name: "Aptos Display", size: 16, bold: true, color: "#FFFFFF" },
  verticalAlignment: "center",
};
readme.getRange("A1:H1").format.rowHeight = 30;

readme.getRange("A3:B15").values = [
  ["Dataset", "Drosophila bulk RNA-seq"],
  ["Comparison", "dMIC60-CS versus dMIC60-WT, analyzed separately by sex"],
  ["Biological replication", "n = 3 pooled libraries per genotype and sex"],
  ["Library input", "Five 96-h pupae of the same genotype and sex per library"],
  ["Groups", "CSF, WTF, CSM, and WTM"],
  ["Model", "Female and male libraries fitted separately with design ~ genotype; WT was the reference level"],
  ["Analysis strategy", "Female and male libraries analyzed separately; no manual count prefilter was applied"],
  ["Gene universe", "17,173 genes entered into each sex-specific DESeq2 analysis"],
  ["Genes with non-missing padj", "Female: 7,456; male: 9,695"],
  ["Statistical test", "DESeq2 Wald test"],
  ["Multiple testing", "Genome-wide Benjamini–Hochberg adjustment"],
  ["Independent filtering", "DESeq2 results() default setting; significance was classified at adjusted P < 0.05"],
  ["Missing values", "Blank statistical cells are missing/not estimable, not zero"],
];
readme.getRange("A3:A15").format = { fill: paleBlue, font: { bold: true, color: navy } };
readme.getRange("A3:B15").format.borders = { preset: "insideHorizontal", style: "thin", color: line };
readme.getRange("A3:B15").format.wrapText = true;

readme.getRange("D3:H7").values = [
  ["Sex", "Up in CS", "Down in CS", "Not significant", "Not testable"],
  ["Female", null, null, null, null],
  ["Male", null, null, null, null],
  ["Total", null, null, null, null],
  ["Direction rule", "padj < 0.05", "padj < 0.05", "Does not meet both cutoffs", "padj missing"],
];
readme.getRange("D3:H3").format = { fill: navy, font: { bold: true, color: "#FFFFFF" }, horizontalAlignment: "center" };
readme.getRange("D4:D7").format = { fill: paleBlue, font: { bold: true, color: navy } };
readme.getRange("E4:H6").format.numberFormat = "#,##0";
for (const [row, sheetName, n] of [[4, "Female DE", femaleRows.length], [5, "Male DE", maleRows.length]]) {
  readme.getRange(`E${row}`).formulas = [[`=COUNTIF('${sheetName}'!$H$2:$H$${n + 1},"Up in CS")`]];
  readme.getRange(`F${row}`).formulas = [[`=COUNTIF('${sheetName}'!$H$2:$H$${n + 1},"Down in CS")`]];
  readme.getRange(`G${row}`).formulas = [[`=COUNTIF('${sheetName}'!$H$2:$H$${n + 1},"Not significant")`]];
  readme.getRange(`H${row}`).formulas = [[`=COUNTIF('${sheetName}'!$H$2:$H$${n + 1},"Not testable")`]];
}
for (const col of ["E", "F", "G", "H"]) readme.getRange(`${col}6`).formulas = [[`=SUM(${col}4:${col}5)`]];
readme.getRange("D6:H6").format = { fill: gray, font: { bold: true, color: "#444444" } };
readme.getRange("D7:H7").format = { fill: amber, font: { color: amberText }, wrapText: true };

readme.getRange("A19:B21").values = [
  ["Classification threshold", "Value"],
  ["Adjusted P value (FDR)", 0.05],
  ["Absolute log2 fold change", 0.58],
];
readme.getRange("A19:B19").format = { fill: navy, font: { bold: true, color: "#FFFFFF" } };
readme.getRange("B20").format.numberFormat = "0.00";
readme.getRange("B21").format.numberFormat = "0.00";

readme.getRange("A24:C32").values = [
  ["DE column", "Definition", "Interpretation"],
  ["gene", "Drosophila gene symbol or identifier", "Text"],
  ["baseMean", "Mean DESeq2 normalized count across the six libraries for that sex", "Expression abundance"],
  ["log2FoldChange", "Estimated log2 fold change for CS versus WT", "Positive = higher in CS"],
  ["lfcSE", "Standard error of the log2-fold-change estimate", "log2 scale"],
  ["stat", "DESeq2 Wald statistic", "Signed test statistic"],
  ["pvalue", "Raw Wald-test P value", "Blank where not estimable"],
  ["padj", "Benjamini–Hochberg-adjusted P value", "Blank where not estimable"],
  ["direction", "Formula-driven classification using the thresholds above", "Up / down / not significant / not testable"],
];
readme.getRange("A24:C24").format = { fill: navy, font: { bold: true, color: "#FFFFFF" } };
readme.getRange("A25:A32").format = { fill: paleBlue, font: { bold: true, color: navy } };
readme.getRange("A24:C32").format.wrapText = true;

readme.getRange("E10:H15").values = [
  ["Sample ID", "Sex", "Genotype", "Replicate"],
  ["CSF1–CSF3", "Female", "dMIC60-CS", "1–3"],
  ["WTF1–WTF3", "Female", "dMIC60-WT", "1–3"],
  ["CSM1–CSM3", "Male", "dMIC60-CS", "1–3"],
  ["WTM1–WTM3", "Male", "dMIC60-WT", "1–3"],
  ["Normalized counts", "All", "All groups", "12 libraries"],
];
readme.getRange("E10:H10").format = { fill: navy, font: { bold: true, color: "#FFFFFF" }, horizontalAlignment: "center" };
readme.getRange("E11:E15").format = { fill: paleBlue, font: { bold: true, color: navy } };

readme.getRange("A35:B38").values = [
  ["Data component", "Analysis source"],
  ["Female differential expression", "DEG_CSF_vs_WRF.csv"],
  ["Male differential expression", "DEG_CSM_vs_WRM.csv"],
  ["Normalized counts", "DESeq2 size-factor-normalized counts; all 12 libraries retained"],
];
readme.getRange("A35:B35").format = { fill: navy, font: { bold: true, color: "#FFFFFF" } };
readme.getRange("A36:A38").format = { fill: paleBlue, font: { bold: true, color: navy } };
readme.getRange("A35:B38").format.wrapText = true;

readme.getRange("A3:H38").format.font = { name: "Aptos", size: 10, color: "#222222" };
for (const range of ["D3:H3", "A19:B19", "A24:C24", "E10:H10", "A35:B35"]) {
  readme.getRange(range).format.font = { name: "Aptos", size: 10, bold: true, color: "#FFFFFF" };
}
readme.getRange("A3:A38").format.columnWidth = 24;
readme.getRange("B3:B38").format.columnWidth = 68;
readme.getRange("C3:C38").format.columnWidth = 34;
readme.getRange("D3:D38").format.columnWidth = 17;
readme.getRange("E3:E38").format.columnWidth = 19;
readme.getRange("F3:H38").format.columnWidth = 18;
readme.getRange("A3:H38").format.verticalAlignment = "center";
readme.freezePanes.freezeRows(1);

function populateDataSheet(sheet, rows, tableName) {
  sheet.showGridLines = false;
  const lastRow = rows.length + 1;
  sheet.getRange("A1:H1").values = [["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "direction"]];
  sheet.getRangeByIndexes(1, 0, rows.length, 7).values = rows;
  sheet.getRange("H2").formulas = [["=IF(ISBLANK(G2),\"Not testable\",IF(AND(G2<'README'!$B$20,C2>='README'!$B$21),\"Up in CS\",IF(AND(G2<'README'!$B$20,C2<=-'README'!$B$21),\"Down in CS\",\"Not significant\")))"]];
  sheet.getRange(`H2:H${lastRow}`).fillDown();
  const table = sheet.tables.add(`A1:H${lastRow}`, true, tableName);
  table.style = "TableStyleMedium2";
  table.showFilterButton = true;
  table.showBandedRows = true;
  sheet.getRange(`A1:H${lastRow}`).format.font = { name: "Aptos", size: 9, color: "#222222" };
  sheet.getRange("A1:H1").format = { fill: navy, font: { name: "Aptos", size: 10, bold: true, color: "#FFFFFF" }, verticalAlignment: "center", wrapText: true };
  sheet.getRange("A1:H1").format.rowHeight = 30;
  sheet.getRange(`B2:B${lastRow}`).format.numberFormat = "#,##0.000";
  sheet.getRange(`C2:E${lastRow}`).format.numberFormat = "0.0000";
  sheet.getRange(`F2:G${lastRow}`).format.numberFormat = "0.00E+00";
  sheet.getRange(`A1:A${lastRow}`).format.columnWidth = 20;
  sheet.getRange(`B1:B${lastRow}`).format.columnWidth = 14;
  sheet.getRange(`C1:E${lastRow}`).format.columnWidth = 17;
  sheet.getRange(`F1:G${lastRow}`).format.columnWidth = 14;
  sheet.getRange(`H1:H${lastRow}`).format.columnWidth = 18;
  sheet.getRange(`H2:H${lastRow}`).conditionalFormats.add("containsText", { text: "Up in CS", format: { fill: green, font: { color: greenText } } });
  sheet.getRange(`H2:H${lastRow}`).conditionalFormats.add("containsText", { text: "Down in CS", format: { fill: red, font: { color: redText } } });
  sheet.getRange(`H2:H${lastRow}`).conditionalFormats.add("containsText", { text: "Not testable", format: { fill: amber, font: { color: amberText } } });
  sheet.freezePanes.freezeRows(1);
  sheet.freezePanes.freezeColumns(1);
}

populateDataSheet(female, femaleRows, "FemaleDETable");
populateDataSheet(male, maleRows, "MaleDETable");

function populateCountsSheet(sheet, normalizedData) {
  sheet.showGridLines = false;
  const rows = normalizedData.values;
  const headers = normalizedData.headers;
  const lastRow = rows.length + 1;
  const lastCol = String.fromCharCode(64 + headers.length);
  sheet.getRange(`A1:${lastCol}1`).values = [headers];
  sheet.getRangeByIndexes(1, 0, rows.length, headers.length).values = rows;
  const table = sheet.tables.add(`A1:${lastCol}${lastRow}`, true, "NormalizedCountsTable");
  table.style = "TableStyleMedium2";
  table.showFilterButton = true;
  table.showBandedRows = true;
  sheet.getRange(`A1:${lastCol}${lastRow}`).format.font = { name: "Aptos", size: 9, color: "#222222" };
  sheet.getRange(`A1:${lastCol}1`).format = { fill: navy, font: { name: "Aptos", size: 10, bold: true, color: "#FFFFFF" }, horizontalAlignment: "center", verticalAlignment: "center" };
  sheet.getRange(`A1:${lastCol}1`).format.rowHeight = 28;
  sheet.getRange(`B2:${lastCol}${lastRow}`).format.numberFormat = "#,##0.000";
  sheet.getRange(`A1:A${lastRow}`).format.columnWidth = 20;
  sheet.getRange(`B1:${lastCol}${lastRow}`).format.columnWidth = 13;
  sheet.freezePanes.freezeRows(1);
  sheet.freezePanes.freezeColumns(1);
}

populateCountsSheet(counts, normalized);

await fs.mkdir(outputDir, { recursive: true });
const exported = await SpreadsheetFile.exportXlsx(workbook);
await exported.save(outputPath);

const checks = [];
for (const [sheetName, range] of [
  ["README", "A1:H38"],
  ["Female DE", "A1:H8"],
  ["Male DE", "A1:H8"],
  ["Normalized counts", "A1:M8"],
]) {
  checks.push((await workbook.inspect({ kind: "table", range: `${sheetName}!${range.split("!").pop()}`, include: "values,formulas", tableMaxRows: 40, tableMaxCols: 14, maxChars: 14000 })).ndjson);
}
checks.push((await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 200 },
  summary: "final formula error scan",
  maxChars: 4000,
})).ndjson);
await fs.writeFile(path.join(outputDir, "Table_S1_verification.ndjson"), checks.join("\n---\n"), "utf8");

for (const [sheetName, range, filename] of [
  ["README", "A1:H38", "preview_readme.png"],
  ["Female DE", "A1:H24", "preview_female.png"],
  ["Male DE", "A1:H24", "preview_male.png"],
  ["Normalized counts", "A1:M18", "preview_normalized_counts.png"],
]) {
  const preview = await workbook.render({ sheetName, range, scale: 1.5, format: "png" });
  await fs.writeFile(path.join(outputDir, filename), new Uint8Array(await preview.arrayBuffer()));
}

console.log(JSON.stringify({ outputPath, femaleRows: femaleRows.length, maleRows: maleRows.length, normalizedRows: normalized.values.length }));
