# Integrative multi-omics identifies coordinated molecular and microbial signatures associated with dietary lipid sources in pigs

## Reproducibility code — author-reviewed methods, code-only repository

This repository organizes the analysis code for a DIABLO (*mixOmics*) integrative
multi-omics study involving **44 immunocastrated male pigs**, with dietary groups
**Soy (18), Canola (14), Fish (12)**. The code is derived from the study's original
R scripts; the fitted model and input matrices were verified separately.

**Scope.** This is a **code-only archive**. Study matrices, individual-level
phenotype records, and the fitted `.rds` are **not included** and should be shared
only after approval. As a result, full end-to-end execution requires separately
obtained input files and a suitable R package environment. The cleaned scripts
have not been executed end-to-end in the target R environment; this is an
important distinction from verification of the **archived model**.

## Repository layout

All R scripts and documentation files are stored in the repository root (top-level folder).

## Script index

| File | Purpose |
| --- | --- |
| `01_candidate_gene_lists.R` | Deduplicate candidate-gene lists from the 3 tissue spreadsheets |
| `02_QC_gene_TPM_Top1000.R` | Filter genes (TPM > 1 in ≥50% of animals), compare across tissues, rank by expression SD and map Ensembl IDs to symbols |
| `03_microbiota_fezes.R` | Fecal ASV prevalence, genus aggregation, CZM and CLR |
| `04_microbiota_ceco.R` | Cecal ASV prevalence, genus aggregation, CZM and CLR |
| `05_fatty_acid_total_sum.R` | Generalized total-sum percentage normalization for muscle, liver and brain |
| `06_phenotype_normality_QC.R` | Shapiro–Wilk screening (p > 0.05) of phenotypes; final selected variables used without log/sqrt transformation |
| `01_fit_DIABLO_FUTURE_RERUN.R` | Model fitting, perf() and tune.block.splsda(); **prospective rerun only** |
| `02_figures_from_archived_model.R` | Plots, loading tables and visualizations using a fitted model |
| `03_Figure5_from_archived_model.R` | Corrected Figure 5 heatmap without fitting |
| `04_environment_versions.R` | Current software version report (not a historical version record) |

## Candidate-gene lists and transcriptome selection

Candidate-gene lists were provided as tissue-specific Excel files. The scripts
above document deduplication, expression-based TPM filtering and subsequent
selection by expression variability. A historical script refers to a DAVID-related
workbook; the exact derivation of the final tissue-specific candidate lists
requires author confirmation. The code does not establish the annotation criteria
used to construct those lists.

Historical source code for the liver list referred to `USARESSE_LIVER.xlsx`, while
the workbook supplied for this package is `Genes_Liver.xlsx`. Their identity
has not been independently established. The preprocessed final model inputs
were verified against the archived model, but candidate-list-to-TPM derivation
was not rerun in this package.

## Microbiota

The two supplied scripts start from **existing ASV count, taxonomy and metadata
tables**, not raw FASTQ files. They apply a 10% ASV prevalence filter, aggregate
at genus level, apply **CZM zero replacement** (`zCompositions`) and **CLR**
(`compositions`). Upstream read processing with **DADA2** is documented in the
related microbiota methodology; DADA2 package source code is not redistributed.
The final 44-sample subset is the one in the archived DIABLO input matrices.

## Fatty acids and phenotypes

- The author confirmed total-sum proportional normalization across **all three
  tissues**. `05_fatty_acid_total_sum.R` expresses the same operation as the
  historical muscle script, generalized to brain/liver/muscle. The 52-feature
  FA block also contains summary quantities/ratios and requires the original
  prepared matrices; the generalized script does not recreate all derived traits.
- The authors selected phenotypic variables screened by **Shapiro–Wilk p > 0.05**;
  the variables retained in the final DIABLO input were **not** log- or
  square-root-transformed. Their Portuguese names were **translated to English**
  for the analysis. The historical QC script also explored transformations,
  but they should not be described as applied to the final PHENO block.
- A non-significant normality test is not evidence of proven normality; this
  screening rule is reported as an author-selected preprocessing criterion,
  not a requirement of DIABLO.

## Archived model parameters (verified)

- **7 blocks**: Genes_MU (616), Genes_LI (731), Genes_BR (608), ASVs_CE (151),
  ASVs_FE (138), FA (52), PHENO (13), each with 44 samples.
- **5 fitted latent components**; Components 1 and 2 used for biological
  interpretation and key plots.
- **Design matrix**: off-diagonal weights `0.7` among the biological blocks;
  `0` on the diagonal; model response (`Y`) connected with weight `1`.
- Featurewise **mean-centering and unit-variance scaling** (`scale=TRUE`).
- Number-of-component performance evaluation: 10-fold repeated M-fold CV,
  100 repetitions (`perf`). Feature-number tuning: 10-fold repeated M-fold CV,
  10 repetitions (`tune.block.splsda`, `max.dist`). The number of components
  was selected from `choice.ncomp$WeightedVote["Overall.BER", "max.dist"]`
  (five components). Final model performance was assessed separately using
  10-fold CV, 20 repetitions and `max.dist`.
- The original script's `circosPlot()` and `network()` use correlation cutoff
  **0.80**, which is distinct from the 0.7 DIABLO design weight.
- The original input matrices were checked against the RDS: exactly matching
  sample and feature names/order and values after z-scoring (within floating-
  point precision, maximum deviations near 1e-15).
- Exact final `keepX` values: `ARCHIVED_MODEL_PARAMETERS.md` and Supplementary
  Table S4. Supplementary Tables S3, S5 and S6 provide feature correlations,
  within-block explained variance and loading values, respectively.
- Block-to-block correlations (Figure 3) are Pearson correlations between
  block-specific latent scores, calculated separately for Components 1 and 2.
  Circos links instead use the DIABLO similarity matrix (absolute cutoff 0.80).
  Both are descriptive; no independent significance testing was performed.
- No additional regression-based adjustment for biological or technical
  covariates was performed on the input features before DIABLO integration.

## Reproducibility and software

The archived RDS reports **R 4.5.1** as the serialization writer version. An
original script comment mentions `mixOmics 6.28.0` used in November 2025;
that note **does not independently establish** the package version at the
August 2026 final-model run. Figures 1–5 were regenerated in the available R 4.5.1 environment using
`mixOmics` 6.32.0; the revised manuscript also reports `igraph` 2.2.1,
`cowplot` 1.2.0, `viridis` 0.6.5, `ggplot2` 4.0.1 and `BiocParallel` 1.42.2
as the available computational environment. `04_environment_versions.R` reports
packages installed in a current environment; it must not be misrepresented
as historical `sessionInfo()`.

**Original cross-validation seed and fold assignments were not recovered.**
A *new* seed appears in the prospective-fitting script for reproducible new
runs, not as a claim to reproduce the original tuning partitions exactly.
The archived model parameters and inputs remain verifiable without rerunning
training; regenerated tuning could select somewhat different values.

## Data availability and execution

Expected input matrices: `musculo.txt`, `figado.txt`, `cerebro.txt`, `ceco.txt`,
`fezes.txt`, `fatty_acids.txt`, `phenos.txt`, `groups.txt` (samples and columns
formatted as in source code). They are not included here. Their access
mechanism should be stated in the accompanying manuscript's data-availability
section; do not upload raw/derived animal-level data or archived model without
appropriate permissions.

The script `01_fit_DIABLO_FUTURE_RERUN.R` expects inputs in `data/` (or
`DIABLO_INPUT_DIR`) and writes **new** results to `results/` (or
`DIABLO_OUTPUT_DIR`); it is not needed to recreate images from an archived RDS.
Run image generation only with the necessary mixOmics packages available.

## Regenerating manuscript figures without refitting

```bash
Rscript 02_figures_from_archived_model.R path/to/final_diablo_model.rds results/figures
```

An optional third argument is the original `all_loadings.txt` table from the
archived run. Its selected feature names and signed loading values are checked
against the model before plotting. Without that table, the script obtains
loading/group-contribution information using `plotLoadings()`. The revised
Figure 1 uses separate panels, Figure 4A excludes the response block Y, and
Figures 1–5 are exported at 600 dpi. Figure 5 can also be generated separately
with `03_Figure5_from_archived_model.R`.

Supplementary Figure S1 is the original 300 dpi image of the initial
100-repeat cross-validation, not a regenerated result. The historical execution
log confirms its origin; the corresponding performance object was not recovered.
No new CV was performed during figure preparation.

Sample metadata for the 44 pigs are provided in Supplementary Table S8 of
the accompanying manuscript. These are not the original model input matrices
and are not included in this code-only repository.

## Documentation status

The source model and its original inputs were verified. During manuscript
revision, the original Component 1 and 2 loading table was checked against the
archived model, block-score correlations were recalculated directly from that
model, and the main figures were re-exported without model refitting. The
consolidated public-facing scripts have **not** been tested end-to-end with all
upstream source tables; specific method
steps (e.g. phenotype name mapping and derived FA calculations) are described
and not presented as independently regenerated. See `READINESS_NOTE.md`.
