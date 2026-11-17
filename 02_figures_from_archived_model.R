# Recreate revised manuscript Figures 1–5 from the archived fitted DIABLO model.
# No model fitting, tuning or cross-validation is performed.
# Usage: Rscript 02_figures_from_archived_model.R final_diablo_model.rds results/figures [all_loadings.txt]
# The optional third argument is the verified historical loading/group-contribution table.
# Consolidated from the figure-specific commands run during the manuscript revision.
# Test the consolidated script in the target R environment before publication.
suppressPackageStartupMessages({
  library(mixOmics)
  library(ggplot2)
  library(cowplot)
  library(viridis)
})
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 2) stop('Supply the archived RDS and a new output directory.')
model <- readRDS(a[1])
out <- a[2]
if (dir.exists(out) && length(list.files(out, all.files=TRUE, no..=TRUE)))
  stop('Output directory is not empty; use a new directory.')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
f <- function(x) file.path(out, x)
blocks <- c('Genes_MU','Genes_LI','Genes_BR','ASVs_CE','ASVs_FE','FA','PHENO')
colors <- c(Soy='purple', Fish='gold', Canola='green')
# Preserve the actual block colors in the revised figures (CE blue, FE brown).
bcolors <- c(Genes_MU='#228B22', Genes_LI='#32CD32', Genes_BR='springgreen',
             ASVs_CE='#4682B4', ASVs_FE='#8B4513', FA='orange', PHENO='red3')
stopifnot(all(blocks %in% names(model$loadings)),
          all(blocks %in% names(model$variates)),
          identical(as.integer(table(factor(model$Y, levels=names(colors)))), c(18L,12L,14L)),
          all(vapply(model$variates[blocks], ncol, integer(1)) == 5L))
# Check loading names and signed values; do not use Y as a biological block.
validate_loadings <- function(tab) {
  for (i in 1:2) for (b in blocks) {
    L <- model$loadings[[b]]
    x <- setNames(L[,i], rownames(L)); x <- x[x != 0]
    t <- subset(tab, Component == i & Tissue == b)
    stopifnot(!anyDuplicated(t$Variable), setequal(names(x),t$Variable),
      isTRUE(all.equal(as.numeric(x[t$Variable]), as.numeric(t$importance), tolerance=1e-7)))
  }
  stopifnot(all(tab$GroupContrib %in% names(colors)))
}
if (length(a) >= 3) {
  tab <- read.delim(a[3], check.names=FALSE)
  tab <- subset(tab, Component %in% 1:2 & Tissue %in% blocks)
} else {
  # The plotLoadings output provides contribution groups under the historical median/max rule.
  # This extraction must be checked against the historical table if available.
  tab <- (function() {
    pdf(NULL, width=12, height=16)
    on.exit(dev.off())
    pieces <- list()
    for (i in 1:2) {
      o <- plotLoadings(model, comp=i, contrib='max', method='median',
                        legend.color=unname(colors), size.legend=.8,
                        legend.title='Groups', size.name=.55,
                        title=paste('Component', i))
      for (b in blocks) {
        d <- as.data.frame(o[[b]])
        d$Tissue <- b; d$Component <- i; d$Variable <- rownames(d)
        pieces[[paste(i,b)]] <- d
      }
    }
    d <- do.call(rbind,pieces)
    subset(d, color != 'white')
  })()
}
validate_loadings(tab)
write.table(tab, f('all_loadings_verified.txt'), sep='\t', row.names=FALSE, quote=FALSE)
for (i in 1:2) {
  d <- subset(tab, Component == i)
  d$key <- reorder(paste(d$Tissue, d$Variable, sep='|'), d$importance)
  p <- ggplot(d, aes(importance,key,color=GroupContrib)) +
    geom_segment(aes(x=0,xend=importance,yend=key),color='darkgray',linewidth=.6) +
    geom_point(size=3) +
    scale_y_discrete(labels=function(x) sub('^[^|]*\\|','',x)) +
    scale_color_manual(values=colors, breaks=c('Canola','Soy','Fish')) +
    facet_grid(Tissue ~ ., scales='free_y',space='free_y') +
    labs(title=paste('Component',i),x='Loading Score',y=NULL,color='Group Contribution') +
    theme_minimal(base_size=12) + theme(axis.text.y=element_text(size=9,color='black'),
      strip.background=element_rect(fill='gray90',color='gray70'),
      strip.text=element_text(face='bold'),plot.title=element_text(hjust=.5,face='bold'),
      legend.position='right')
  ggsave(f(paste0('Figure1',LETTERS[i],'_final.png')),p,width=9,height=13,units='in',dpi=600,bg='white')
}
# Figure 2: same parameters as the revised server export.
png(f('Figure2.png'),width=9,height=7.5,units='in',res=600,bg='white')
correlations <- circosPlot(model,comp=c(1,2),cutoff=.8,line=TRUE,
  color.blocks=unname(bcolors[blocks]),color.cor=c('lightblue','orange'),
  size.labels=1,color.Y=unname(colors),size.variables=.45)
dev.off()
write.table(correlations,f('circos_similarity_matrix.txt'),sep='\t')
# Figure 3: panel renders at 600 dpi before composition; export score correlations.
for (i in 1:2) {
  png(f(paste0('Figure3_comp',i,'.png')),width=10,height=8,units='in',res=600,bg='white')
  plotDiablo(model,ncomp=i,col.per.group=unname(colors))
  dev.off()
  scores <- sapply(model$variates[blocks],function(x)x[,i])
  write.table(cor(scores),f(paste0('block_correlations_comp',i,'.txt')),sep='\t')
}
figure3 <- plot_grid(ggdraw()+draw_image(f('Figure3_comp1.png')),
                     ggdraw()+draw_image(f('Figure3_comp2.png')),ncol=1,labels=c('A','B'))
ggsave(f('Figure3.png'),figure3,width=10,height=16,dpi=600,bg='white')
# Figure 4A: only the seven biological blocks, matching the manuscript order.
order4 <- c('ASVs_CE','ASVs_FE','FA','Genes_BR','Genes_LI','Genes_MU','PHENO')
scores <- do.call(rbind,lapply(order4,function(b) data.frame(
  Dimension1=model$variates[[b]][,1],Dimension2=model$variates[[b]][,2],
  Group=factor(model$Y,levels=names(colors)),Block=b)))
scores$Block <- factor(scores$Block,levels=order4)
figure4A <- ggplot(scores,aes(Dimension1,Dimension2,color=Group)) +
  geom_point(shape=1,size=2.2,stroke=.8) + stat_ellipse(type='norm',level=.95,linewidth=.9) +
  facet_wrap(~Block,ncol=3,scales='fixed',labeller=label_both) +
  scale_color_manual(values=colors) + coord_cartesian(xlim=c(-12,10),ylim=c(-8,10)) +
  labs(x=NULL,y=NULL,color='Groups') + theme_minimal(base_size=11) +
  theme(strip.background=element_rect(fill='gray85',color='gray40'),
        strip.text=element_text(face='bold'),legend.position='right')
ggsave(f('Figure4A_final.png'),figure4A,width=9,height=8,dpi=600,bg='white')
figure4B <- plotArrow(model,group=model$Y,ind.names=FALSE,legend=TRUE,title='',
  col.per.group=unname(colors),legend.title='Groups',style='ggplot2')
if (!inherits(figure4B,'ggplot')) figure4B <- figure4B$graph
stopifnot(inherits(figure4B,'ggplot'))
figure4B <- figure4B + scale_color_manual(values=colors,breaks=names(colors),drop=FALSE)
ggsave(f('Figure4B.png'),figure4B,width=9,height=8,dpi=600,bg='white')
figure4 <- plot_grid(figure4A,ggdraw()+draw_image(f('Figure4B.png')),ncol=2,labels=c('A','B'))
ggsave(f('Figure4_final.png'),figure4,width=18,height=8,dpi=600,bg='white')
# Figure 5: preserve heatmap color mapping and sorted group labels.
png(f('Figure5.png'),width=12,height=15,units='in',res=600,bg='white')
cimDiablo(model,color.blocks=unname(bcolors[blocks]),comp=c(1,2),legend.position='topright',
  color.Y=unname(colors[sort(unique(as.character(model$Y)))]),transpose=TRUE,row.names=FALSE,
  size.legend=.8,color=viridis(n=100,option='magma'),trim=FALSE)
dev.off()
# S1 is the historical initial-model CV figure; it is intentionally not regenerated here.
cat('Figure exports completed. Compare them to the final manuscript before use.\n')
