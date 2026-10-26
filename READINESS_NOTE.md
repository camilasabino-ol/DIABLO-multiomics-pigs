# Publication and execution note

This folder contains the **public-facing source-code selection only**: no raw
animal data, RDS, images, or personal directory structures. The scripts are
organized and ready for author review/upload as **code**, but **full
computational re-execution remains unverified**, largely because the R
installation available during code organization lacked `mixOmics`.

Authors have confirmed:
- DAVID functional annotation was used for candidate genes.
- Identical total-sum FA normalization across three tissues.
- PHENO feature selection used Shapiro–Wilk screening and did not use
  transformations; phenotype labels were translated to English.

Do not represent these files as a completely tested raw-data-to-figure pipeline.
The DADA2 raw-sequence upstream stage and DAVID web output are referenced as
external upstream procedures, not bundled software sources. Output packages
and data availability can be supplied by author-approved means.

Original CV random seed/folds are unknown. RDS writer R version is 4.5.1,
while mixOmics 6.28.0 appears in an older source comment and is unconfirmed
for the archived August 2026 fit. Preserve the historical model and datasets.
