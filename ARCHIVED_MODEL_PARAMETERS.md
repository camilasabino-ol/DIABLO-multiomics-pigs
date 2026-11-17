# Archived model parameters (from final_diablo_model.rds)

- Serialization writer R version: **4.5.1**; does **not** establish mixOmics historical package version.
- Model class: `block.splsda`; `scale = TRUE`.
- 44 pigs: Soy 18, Fish 12, Canola 14.
- 7 blocks: Genes_MU (616), Genes_LI (731), Genes_BR (608), ASVs_CE (151), ASVs_FE (138), FA (52), PHENO (13).
- 5 fitted components; interpreted components 1 and 2.
- Off-diagonal biological design weights: 0.7; relationships to response Y: 1.
- Historical fitting code: `perf()` Mfold, 10 folds, 100 repeats; `tune.block.splsda()` Mfold, 10 folds, 10 repeats, `dist='max.dist'`, 5 workers.
- The original seed/fold assignments are not available. The execution log confirms 5 components selected using weighted-vote Overall.BER with max.dist; final-model performance used 10-fold CV, 20 repeats, max.dist.
- **Figure 2 Circos and network code:** `cutoff = 0.80`, **not** the design weight 0.7.

| Block | KeepX comp1 | comp2 | comp3 | comp4 | comp5 |
|---|---:|---:|---:|---:|---:|
| Genes_MU | 18 | 12 | 12 | 8 | 12 |
| Genes_LI | 8 | 8 | 20 | 20 | 6 |
| Genes_BR | 12 | 20 | 12 | 8 | 20 |
| ASVs_CE | 13 | 4 | 4 | 13 | 4 |
| ASVs_FE | 4 | 10 | 10 | 4 | 10 |
| FA | 8 | 6 | 8 | 10 | 10 |
| PHENO | 10 | 8 | 8 | 4 | 10 |

**Verification completed:** all original seven input matrices match the archived RDS X matrices after per-feature z-scoring using sample standard deviation (max difference ~1e-15), with identical feature and sample labels and order. Does not prove full rerun/tuning reproducibility.

## Figure verification during revision

- Components 1 and 2: 73 and 68 nonzero loadings, respectively, verified against the archived original loading table.
- Pearson correlations between cecal and fecal block scores: 0.305 (Component 1), 0.495 (Component 2), rounded to 0.31 and 0.50 in the manuscript.
- Strongest pair in Component 1: liver transcriptome and fecal microbiota (approximately 0.715); in Component 2: muscle transcriptome and fecal microbiota (approximately 0.689).
- Figure 3 was regenerated from the archived model and manuscript correlation values were updated accordingly.
- Supplementary Figure S1 retains the original initial-model CV plot (300 dpi).
