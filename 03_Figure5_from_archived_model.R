############################################################
# Figure 5 — DIABLO heatmap
# Recreates Figure 5 from the previously fitted final model.
# Does not rerun or modify the DIABLO model.
#
# Usage (from a terminal):
#   Rscript 03_Figure5_from_archived_model.R path/to/final_diablo_model.rds path/to/Figure5_heatmap_corrected.png
#
# Both arguments are optional if the RDS is in the current folder.
############################################################

library(mixOmics)
library(viridis)

args <- commandArgs(trailingOnly = TRUE)
model_path <- if (length(args) >= 1) args[[1]] else "final_diablo_model.rds"
figure_path <- if (length(args) >= 2) args[[2]] else "Figure5_heatmap_corrected.png"

if (!file.exists(model_path)) {
  stop("Model file not found: ", model_path,
       ". Supply the path to the fitted final_diablo_model.rds as the first argument.")
}

final.diablo.model <- readRDS(model_path)

# Named treatment colors, independent of the order chosen by cimDiablo.
treatment_colors <- c(Soy = "purple", Fish = "gold", Canola = "green")
heatmap_group_order <- sort(unique(as.character(final.diablo.model$Y)))

if (!all(heatmap_group_order %in% names(treatment_colors))) {
  stop("Unrecognized treatment(s): ",
       paste(setdiff(heatmap_group_order, names(treatment_colors)), collapse = ", "))
}

heatmap_treatment_colors <- unname(treatment_colors[heatmap_group_order])
print(data.frame(Group = heatmap_group_order,
                 Color = heatmap_treatment_colors))

dir.create(dirname(figure_path), showWarnings = FALSE, recursive = TRUE)

png(filename = figure_path, width = 12, height = 15,
    units = "in", bg = "white", res = 300)

cimDiablo(
  final.diablo.model,
  color.blocks = c(
    "#228B22",     # Genes_MU
    "#32CD32",     # Genes_LI
    "springgreen", # Genes_BR
    "#4682B4",     # ASVs_FE
    "#8B4513",     # ASVs_CE
    "orange",      # FA
    "red3"         # PHENO
  ),
  comp = c(1, 2),
  legend.position = "topright",
  color.Y = heatmap_treatment_colors,
  transpose = TRUE,
  row.names = FALSE,
  size.legend = 0.8,
  color = viridis(n = 100, option = "magma"),
  trim = FALSE
)

dev.off()
cat("Figure 5 generated: ", figure_path, "\n", sep = "")
