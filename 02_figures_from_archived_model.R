# DIABLO figures and loading tables from archived model (DRAFT)
# Derived from the original submitted visualization code, preserving 0.80
# correlation cutoffs and the final model's first two components.
#
# This is not a fresh DIABLO fit and does NOT modify any original files.
# Run with: Rscript scripts/02_figures_from_archived_model.R \
#   data/final_diablo_model.rds results/figures
# Only run in an R environment with mixOmics and plot dependencies installed.
# Visual outputs need comparison with the manuscript before publication.

suppressPackageStartupMessages({
  library(mixOmics)
  library(ggplot2)
  library(igraph)
  library(cowplot)
  library(viridis)
})

args <- commandArgs(trailingOnly=TRUE)
model_path <- if (length(args) >= 1) args[[1]] else 'data/final_diablo_model.rds'
output_dir <- if (length(args) >= 2) args[[2]] else 'results/figures'
if (!file.exists(model_path)) stop('Missing RDS file: ', model_path)
final.diablo.model <- readRDS(model_path)
if (!inherits(final.diablo.model, 'block.splsda')) stop('Unexpected RDS object class')
Y <- final.diablo.model$Y
if (!all(c('Soy','Fish','Canola') %in% unique(as.character(Y)))) {
  stop('Expected Soy, Fish, Canola groups')
}
data <- final.diablo.model$X
# The original code uses this named treatment color convention.
treatment_colors <- c('Soy'='purple','Fish'='gold','Canola'='green')
# All image files are written to the chosen output folder, not to original data.
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
setwd(output_dir)

#plots
#component1
png("omics_correlation_comp1.png", width = 10, height = 8,bg = "white",units = "in",res = 300)

plotDiablo(final.diablo.model, ncomp = 1,col.per.group = c("purple","gold","green") ) 

dev.off()

#component2
png("omics_correlation_comp2.png", width = 10, height = 8,bg = "white",units = "in",res = 300)
plotDiablo(final.diablo.model, ncomp = 2,col.per.group = c("purple","gold","green"))
dev.off()


# Saving those figures as only one figure
# Getting the images from our directory
d1 <- ggdraw() + draw_image("omics_correlation_comp1.png")
d2 <- ggdraw() + draw_image("omics_correlation_comp2.png")

combined_plot <- plot_grid(d1, d2, ncol = 1, labels = c("A", "B"))

# Saving the Figure

ggsave("omics_correlation_combined_comp1_2.png", combined_plot, width = 10, height = 16, dpi = 300,bg = "white")


# The sample plot with the plotIndiv() function projects each sample into the space spanned by the components from each block, resulting in a series of graphs corresponding to each data set
explained.variance <- as.data.frame(final.diablo.model$prop_expl_var)
write.table(explained.variance,row.names = TRUE,file = "explained_variance_blocks.txt",sep = "\t")

############################################################
# FIGURE 4A
# Sample projection by omics block
# Rebuilt manually to ensure that points and ellipses use
# the same treatment colors in the archived plotting setup
############################################################

blocks_to_plot <- names(data)

plot_indiv_data <- do.call(
  rbind,
  lapply(blocks_to_plot, function(block_name) {
    
    block_scores <- final.diablo.model$variates[[block_name]]
    
    data.frame(
      Sample = rownames(block_scores),
      Component1 = block_scores[, 1],
      Component2 = block_scores[, 2],
      Block = block_name,
      Group = Y,
      stringsAsFactors = FALSE
    )
  })
)

plot_indiv_data$Group <- factor(
  plot_indiv_data$Group,
  levels = c("Soy", "Fish", "Canola")
)

plot_indiv_data$Block <- factor(
  plot_indiv_data$Block,
  levels = blocks_to_plot
)

e1 <- ggplot(
  plot_indiv_data,
  aes(
    x = Component1,
    y = Component2,
    color = Group
  )
) +
  geom_point(
    size = 1.5,
    alpha = 0.85
  ) +
  stat_ellipse(
    aes(group = Group),
    type = "norm",
    level = 0.95,
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ Block,
    scales = "free",
    ncol = 4
  ) +
  scale_color_manual(
    values = treatment_colors,
    breaks = c("Soy", "Fish", "Canola"),
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = NULL,
    color = "Groups"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(
      fill = "gray90",
      color = "gray70"
    ),
    strip.text = element_text(
      face = "bold",
      size = 9
    ),
    axis.text = element_text(size = 7),
    axis.title = element_blank(),
    legend.position = "right",
    legend.title = element_text(
      face = "bold",
      size = 9
    ),
    legend.text = element_text(size = 8),
    panel.grid.major = element_line(
      color = "gray90",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(
  filename = "block_plsda_blocked.png",
  plot = e1,
  width = 9,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

# In this next Figure I have returned to my final model and sum the total explained variance of each component and add here as x and y axis

############################################################
# FIGURE 4B
# Arrow plot with standardized treatment colors
############################################################

e2 <- plotArrow(
  final.diablo.model,
  group = Y,
  ind.names = FALSE,
  legend = TRUE,
  title = "",
  col.per.group = unname(
    treatment_colors[c("Soy", "Fish", "Canola")]
  ),
  legend.title = "Groups",
  style = "ggplot2"
)

e2 <- e2 +
  scale_color_manual(
    values = treatment_colors,
    breaks = c("Soy", "Fish", "Canola"),
    labels = c("Soy", "Fish", "Canola"),
    drop = FALSE
  )

ggsave(
  filename = "block_plsda_unblocked.png",
  plot = e2,
  width = 9,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)


#Joining the figures
e3 <- ggdraw() + draw_image("block_plsda_blocked.png")
e4 <- ggdraw() + draw_image("block_plsda_unblocked.png")

combined_plot1 <- plot_grid(e3, e4, ncol = 2, labels = c("A", "B"))


ggsave("multi-block_sPLS-DA_plots.png", combined_plot1, width = 18, height = 8, dpi = 300,bg = "white")

#Variable plots
png(filename = "correlation_circleplot.png",width = 8,height = 10,units = "in",bg = "white",res = 300)

f1_transc <- plotVar(final.diablo.model, var.names = TRUE,
                                  style = 'graphics', legend = TRUE,
                                  pch = rep(16, 7),
                                  cex = rep(1.0, 7),
                                  col = c('#228B22',   # Genes_MU
                                          '#32CD32',   # Genes_LI
                                          'springgreen', # Genes_BR
                                          '#4682B4',   # ASVs_FE
                                          '#8B4513',   # ASVs_CE
                                          'orange',    # FA
                                          'red3'),     # PHENO
                                  title = "",
                                  overlap = FALSE,
                                  blocks = c("Genes_MU","Genes_LI","Genes_BR", "ASVs_FE","ASVs_CE","FA","PHENO"),
                                  legend.title = "Blocks")
dev.off()


### CIRCOS PLOT COM 7 BLOCOS (incluindo FA e PHENO)
png(
  filename = "circosplot.png",
  width = 9, height = 7.5,
  units = "in", bg = "white", res = 300
)

f2 <- circosPlot(
  final.diablo.model,
  comp = c(1, 2),   # <- ADICIONEI ESTA LINHA
  cutoff = 0.80,
  line = TRUE,
  color.blocks = c(
    '#228B22',   # Genes_MU
    '#32CD32',   # Genes_LI
    'springgreen', # Genes_BR
    '#4682B4',   # ASVs_FE
    '#8B4513',   # ASVs_CE
    'orange',    # FA
    'red3'       # PHENO
  ),
  color.cor = c("lightblue", "orange"),
  size.labels = 1.0,
  color.Y = c("purple", "gold", "green"),
  size.variables = 0.45
)

dev.off()

# Exportando a matriz de correlação do circos
write.table(f2, "correlation_biomarkers.txt", sep = "\t")

library(viridis)
paleta_viridis <- viridis(n = 100, option = "A")

png(
  filename = "network1.png",
  width = 15, height = 9,
  units = "in", bg = "white", res = 400
)

f3 <- network(
  final.diablo.model,
  blocks = 1:7,  # agora temos 7 blocos em 'data'
  color.node = c(
    '#228B22',    # Genes_MU
    '#32CD32',    # Genes_LI
    'springgreen',# Genes_BR
    '#4682B4',    # ASVs_FE
    '#8B4513',    # ASVs_CE
    'orange',     # FA
    'red3'        # PHENO
  ),
  cutoff = 0.80,
  graph.scale = 0.09,
  cex.node.name = 0.6,
  layout.fun = layout_in_circle,
  shape.node = "rectangle",
  size.node = 0.8,
  lwd.edge = c(1, 1.5),
  color.edge = paleta_viridis
)

dev.off()

# Exportar rede para Cytoscape
myNetwork <- f3
write_graph(myNetwork$gR, file = "myNetwork.gml", format = "gml")

#Joining the figures
f4 <- ggdraw() + draw_image("correlation_circleplot.png")
f5 <- ggdraw() + draw_image("circosplot.png")

combined_plot2 <- plot_grid(f4,f5, ncol = 2, labels = c("A", "B"),rel_heights = c(1,2))


ggsave("correlation_circos_network.png", combined_plot2, width = 10, height = 6, dpi = 500,bg = "white")



png(filename = "biomarkers_candidates_1.png",width = 8,height = 12,units = "in",bg = "white",res = 300)
g1 <- plotLoadings(final.diablo.model, comp = 1, contrib = 'max', method = 'median',legend.color = c("purple", "gold", "green"),size.legend = 0.8,legend.title = "Groups",size.name = 0.55,title = "Component 1")

dev.off()


# Exporting tables

g1_table <- as.data.frame(g1$Genes_MU)
g1_table$Tissue <- "Genes_MU"
g1_table$Component <- "1"
g1_table$Variable <- rownames(g1_table)

write.table(g1_table,file = "comp1_variable_importance_genes_MU.txt",sep = "\t",row.names = TRUE)

g1_table2 <- as.data.frame(g1$Genes_LI)
g1_table2$Tissue <- "Genes_LI"
g1_table2$Component <- "1"
g1_table2$Variable <- rownames(g1_table2)
write.table(g1_table2,file = "comp1_variable_importance_genes_LI.txt",sep = "\t",row.names = TRUE)

g1_table3 <- as.data.frame(g1$Genes_BR)
g1_table3$Tissue <- "Genes_BR"
g1_table3$Component <- "1"
g1_table3$Variable <- rownames(g1_table3)
write.table(g1_table3,file = "comp1_variable_importance_genes_BR.txt",sep = "\t")

g1_table4 <- as.data.frame(g1$ASVs_FE)
g1_table4$Tissue <- "ASVs_FE"
g1_table4$Component <- "1"
g1_table4$Variable <- rownames(g1_table4)
write.table(g1_table4,file = "comp1_variable_importance_ASVs_FE.txt",sep = "\t")

g1_table5 <- as.data.frame(g1$ASVs_CE)
g1_table5$Tissue <- "ASVs_CE"
g1_table5$Component <- "1"
g1_table5$Variable <- rownames(g1_table5)
write.table(g1_table5,file = "comp1_variable_importance_ASVs_CE.txt",sep = "\t")

### NOVO: FA e PHENO na componente 1
g1_table6 <- as.data.frame(g1$FA)
g1_table6$Tissue <- "FA"
g1_table6$Component <- "1"
g1_table6$Variable <- rownames(g1_table6)
write.table(g1_table6, file = "comp1_variable_importance_FA.txt", sep = "\t")

g1_table7 <- as.data.frame(g1$PHENO)
g1_table7$Tissue <- "PHENO"
g1_table7$Component <- "1"
g1_table7$Variable <- rownames(g1_table7)
write.table(g1_table7, file = "comp1_variable_importance_PHENO.txt", sep = "\t")

# Second component

png(filename = "biomarkers_candidates_2.png",width = 8,height = 6,units = "in",bg = "white",res = 300)

g2 <- plotLoadings(final.diablo.model, comp = 2, contrib = 'max', method = 'median',legend.color = c("purple", "gold", "green"),size.legend = 0.8,legend.title = "Groups",size.name = 0.55,title = "Component 2")

dev.off()


# Exporting tables

g2_table <- as.data.frame(g2$Genes_MU)
g2_table$Tissue <- "Genes_MU"
g2_table$Component <- "2"
g2_table$Variable <- rownames(g2_table)
write.table(g2_table,file = "comp2_variable_importance_genes_MU.txt",sep = "\t",row.names = TRUE)

g2_table2 <- as.data.frame(g2$Genes_LI)
g2_table2$Tissue <- "Genes_LI"
g2_table2$Component <- "2"
g2_table2$Variable <- rownames(g2_table2)
write.table(g2_table2,file = "comp2_variable_importance_genes_LI.txt",sep = "\t",row.names = TRUE)

g2_table3 <- as.data.frame(g2$Genes_BR)
g2_table3$Tissue <- "Genes_BR"
g2_table3$Component <- "2"
g2_table3$Variable <- rownames(g2_table3)
write.table(g2_table3,file = "comp2_variable_importance_Genes_BR.txt",sep = "\t")

g2_table4 <- as.data.frame(g2$ASVs_FE)
g2_table4$Tissue <- "ASVs_FE"
g2_table4$Component <- "2"
g2_table4$Variable <- rownames(g2_table4)
write.table(g2_table4,file = "comp2_variable_importance_ASVs_FE.txt",sep = "\t")

g2_table5 <- as.data.frame(g2$ASVs_CE)
g2_table5$Tissue <- "ASVs_CE"
g2_table5$Component <- "2"
g2_table5$Variable <- rownames(g2_table5)
write.table(g2_table5,file = "comp2_variable_importance_ASVs_CE.txt",sep = "\t")

### NOVO: FA e PHENO na componente 2
g2_table6 <- as.data.frame(g2$FA)
g2_table6$Tissue <- "FA"
g2_table6$Component <- "2"
g2_table6$Variable <- rownames(g2_table6)
write.table(g2_table6, file = "comp2_variable_importance_FA.txt", sep = "\t")

g2_table7 <- as.data.frame(g2$PHENO)
g2_table7$Tissue <- "PHENO"
g2_table7$Component <- "2"
g2_table7$Variable <- rownames(g2_table7)
write.table(g2_table7, file = "comp2_variable_importance_PHENO.txt", sep = "\t")

# Joining all the loadings into one table
all_loadings <- rbind(
  g1_table, g1_table2, g1_table3, g1_table4, g1_table5, g1_table6, g1_table7,
  g2_table, g2_table2, g2_table3, g2_table4, g2_table5, g2_table6, g2_table7
)

all_loadings <- all_loadings[all_loadings$color != "white", ]
write.table(all_loadings, file = "all_loadings.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

# Making a better figure to the loadings

all_loadings1 <- all_loadings[all_loadings$Component == "1", ]

# Now the plot
comp1_loading <- ggplot(all_loadings1, aes(x = importance, y = reorder(Variable, importance), color = GroupContrib)) +
  
  # Lollipop sticks
  geom_segment(aes(x = 0, xend = importance, y = Variable, yend = Variable), 
               color = "darkgray", size = 1) +
  
  # Lollipop heads
  geom_point(aes(fill = color), size = 5, shape = 21) +
  
  # Variable labels
  geom_label(aes(label = Variable, x = -max(all_loadings1$importance) * 1.20), 
             size = 4.0, fill = "gray90", label.size = 0.5, color = "black") +
  
  # Custom colors
  scale_color_manual(values = c("Canola" = "green", "Soy" = "purple", "Fish" = "gold"),
                     name = "Group Contribution",
                     limits = c("Canola", "Soy", "Fish")) +
  scale_fill_identity() +
  
  theme_minimal() +
  
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right",
    panel.grid.major = element_line(color = "lightgray", size = 0.75),
    panel.grid.minor = element_line(color = "lightgray", size = 0.30),
    strip.background = element_rect(fill = "gray90", color = "gray70"),
    strip.text = element_text(face = "bold", color = "black", size = 14),
    
    # 🔑 Title and legend size adjustments
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"), # Bigger centered title
    legend.title = element_text(size = 17, face = "bold"),             # Bigger legend title
    legend.text  = element_text(size = 16)                             # Bigger legend labels
  ) +
  
  labs(
    title = "Component 1",
    x = "Loading Score",
    y = NULL,
    color = "Group Contribution"
  ) +
  
  facet_grid(rows = vars(Tissue), scales = "free_y", space = "free")

ggsave(filename = "loading_plot_comp1.png",plot = comp1_loading,device = "png",dpi = 300,width = 30,height = 18,bg = "white")

###### For component 2
# Making a better figure to the loadings

all_loadings2 <- all_loadings[all_loadings$Component == "2", ]

# Now, create the plot
# Now the plot
comp2_loading <- ggplot(all_loadings2, aes(x = importance, y = reorder(Variable, importance), color = GroupContrib)) +
  
  # Lollipop sticks
  geom_segment(aes(x = 0, xend = importance, y = Variable, yend = Variable), 
               color = "darkgray", size = 1) +
  
  # Lollipop heads
  geom_point(aes(fill = color), size = 5, shape = 21) +
  
  # Variable labels
  geom_label(aes(label = Variable, x = -max(all_loadings1$importance) * 1.32), 
             size = 4.0, fill = "gray90", label.size = 0.5, color = "black") +
  
  # Custom colors
  scale_color_manual(values = c("Canola" = "green", "Soy" = "purple", "Fish" = "gold"),
                     name = "Group Contribution",
                     limits = c("Canola", "Soy", "Fish")) +
  scale_fill_identity() +
  
  theme_minimal() +
  
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right",
    panel.grid.major = element_line(color = "lightgray", size = 0.75),
    panel.grid.minor = element_line(color = "lightgray", size = 0.30),
    strip.background = element_rect(fill = "gray90", color = "gray70"),
    strip.text = element_text(face = "bold", color = "black", size = 14),
    
    # 🔑 Title and legend size adjustments
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"), # Bigger centered title
    legend.title = element_text(size = 17, face = "bold"),             # Bigger legend title
    legend.text  = element_text(size = 16)                             # Bigger legend labels
  ) +
  
  labs(
    title = "Component 2",
    x = "Loading Score",
    y = NULL,
    color = "Group Contribution"
  ) +
  
  facet_grid(rows = vars(Tissue), scales = "free_y", space = "free")

ggsave(filename = "loading_plot_comp2.png",plot = comp2_loading,device = "png",dpi = 300,width = 24,height = 16,bg = "white")


# Joining the Figures

#g3 <- ggdraw() + draw_image("loading_plot_comp1.png")
#g4 <- ggdraw() + draw_image("loading_plot_comp2.png")


#combined_plot4 <- plot_grid(g3, g4, ncol = 2, labels = c("A", "B"))


#ggsave("biomarkers_candidates_joined_top10.png", combined_plot4, width = 12, height = 7, dpi = 300,bg = "white")




#windows(width = 14, height = 10) 

cat('Figure and loading-table generation commands finished.\n')
cat('Compare visually to the published versions before public release.\n')
