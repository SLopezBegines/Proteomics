# Running the Proteomics Pipeline

This directory contains the main RMarkdown notebooks that orchestrate the full analysis.
Each notebook sources the modular scripts in `../code/` and produces a self-contained
HTML report alongside all tables and figures under `output/`.

---

## Files

| File | Description |
|------|-------------|
| `IP-CLN3_PXD031582.Rmd` | Example analysis: CLN3 lysosomal interactome, dataset PXD031582. Use this as the template for new analyses. |
| `rbioapi/` | STRING network images downloaded programmatically during the last analysis run. Regenerated automatically on each run. |

---

## Prerequisites

### 1. R environment

Restore the exact package environment before running:

```r
# From the project root (where renv.lock is located)
install.packages("renv")
renv::restore()
```

This installs all CRAN and Bioconductor packages at the versions recorded in `renv.lock`.
Expected R version: **≥ 4.3**. See `renv.lock` for exact tested versions.

### 2. Input data

Place the MaxQuant output file in the `rawdata/` folder at the project root:

```
Proteomics/
└── rawdata/
    └── ProteinGroups_{your_experiment}.xlsx   ← MaxQuant ProteinGroups export
```

MaxQuant must be run with **LFQ quantification enabled**. The pipeline expects a standard
`proteinGroups.txt` converted to `.xlsx`, or the direct Excel export from Perseus.
Required columns: `Protein IDs`, `Gene names`, `LFQ intensity *`, `Reverse`, `Potential contaminant`, `Only identified by site`.

---

## Adapting to a New Dataset

Open the `.Rmd` in RStudio and modify **only the setup chunk** (lines ~30–120). All other
sections run automatically from this configuration.

### Step 1 — Organism parameters

```r
# KEGG species code:
#   "hsa" = Homo sapiens
#   "mmu" = Mus musculus
#   "dre" = Danio rerio
kegg_organism <- "hsa"

# NCBI taxonomy ID:
#   9606  = Human
#   10090 = Mouse
#   7955  = Zebrafish
species <- 9606

# Bioconductor OrgDb annotation package:
#   "org.Hs.eg.db" = Human
#   "org.Mm.eg.db" = Mouse
#   "org.Dr.eg.db" = Zebrafish
organism <- "org.Hs.eg.db"
```

Install the OrgDb package if not already present:
```r
BiocManager::install("org.Hs.eg.db")  # or org.Mm.eg.db, etc.
```

### Step 2 — Output path

```r
# Output will be created at: ./output/YourExperimentName/
output_path <- "./output/YourExperimentName/"
image_number <- 1  # Reset to 1 for each new analysis
```

### Step 3 — Pairwise comparisons

```r
# Each string defines a limma contrast: "ConditionA_vs_ConditionB"
# ConditionA and ConditionB must match values in the 'condition' vector below.
comparisons <- c("KO_vs_WT", "KO_vs_CTRL", "WT_vs_CTRL")
```

### Step 4 — Sample mapping

```r
# label: exact LFQ column names from the MaxQuant ProteinGroups file
label <- c(
  "LFQ intensity Sample_A_01", "LFQ intensity Sample_A_02",
  "LFQ intensity Sample_B_01", "LFQ intensity Sample_B_02"
)

# columns_to_rename: short names for internal use (no spaces)
columns_to_rename <- c("A_1", "A_2", "B_1", "B_2")

# condition: group assignment (must align with label order)
condition <- as.factor(c("A", "A", "B", "B"))

replicate <- c(1, 2, 1, 2)
experiment <- rep("MyExperiment", length(replicate))
```

### Step 5 — Input data path

```r
prot_data <- readxl::read_xlsx(path = "../rawdata/ProteinGroups_YourExperiment.xlsx")
```

---

## Running the Analysis

### Option A — Knit from RStudio (recommended)

1. Open `IP-CLN3_PXD031582.Rmd` (or your adapted copy) in RStudio
2. Click **Knit** → **Knit to HTML**
3. The HTML report and all output files are generated automatically

### Option B — Command line

```bash
cd mains/
Rscript -e "rmarkdown::render('IP-CLN3_PXD031582.Rmd', output_format='html_document')"
```

### Option C — Run chunks interactively

Open the `.Rmd` in RStudio and run chunks sequentially with **Ctrl+Shift+Enter**.
Useful for debugging or inspecting intermediate objects.

---

## Expected Runtime

| Dataset size | Approximate runtime |
|--------------|---------------------|
| ~1,000 proteins, 2–3 comparisons | 5–15 min |
| ~3,000 proteins, 4–6 comparisons | 20–40 min |
| Large datasets with API queries (STRING, PANTHER) | 30–60 min |

Runtime is dominated by functional enrichment (clusterProfiler) and API calls (rbioapi).
Internet connection required for STRING and PANTHER queries.

---

## Output

All results are written to `output/{analysis_name}/`. See [`../docs/readme.md`](../docs/readme.md)
for a full description of every table, figure, and RData file produced.

---

## Troubleshooting

**`Error in read_xlsx: cannot open the connection`**
→ Check that the path in `read_xlsx()` points to an existing file. Use `file.exists(path)` to verify.

**`Error: object 'kegg_organism' not found`**
→ The setup chunk must be run first. In interactive mode, run chunks in order from the top.

**`clusterProfiler: no enrichment found`**
→ No terms passed the adjusted p-value threshold. Try relaxing `pvalueCutoff` in `06_GO.R`,
or verify that the gene IDs are correctly mapped (check `keyType` in `global_variables.R`).

**`rbioapi: Error in curl`**
→ Network timeout during STRING/PANTHER API call. Re-run the `rbioapi-queries` chunk.
The API has rate limits; wait 30 seconds between retries.

**`renv::restore() fails on Bioconductor packages`**
→ Run `BiocManager::install(version = "3.18")` first to set the correct Bioconductor release,
then retry `renv::restore()`.

**Figures are blank or show wrong groups**
→ Verify that `condition` and `label` vectors have matching lengths and that comparison strings
(e.g., `"KO_vs_WT"`) exactly match factor levels in `condition`.
