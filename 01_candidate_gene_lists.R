# 01_candidate_gene_lists.R — candidate gene lists per tissue (DRAFT)
# This reproduces the de-duplication logic of Salvando_ListadeGenes.R.
# Candidate genes were functionally annotated using DAVID (confirmed by authors).
# Assumes columns 'Genes' contain one or more IDs separated by commas.
# DAVID functional annotation was an upstream step, not reproduced by this script.
# This script begins with the three author-supplied candidate-gene spreadsheets.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
})
input_dir <- Sys.getenv('GENE_LIST_INPUT_DIR', unset='data/candidate_genes')
output_dir <- Sys.getenv('GENE_LIST_OUTPUT_DIR', unset='results/preprocessing/gene_lists')
files <- c(Brain='Genes_Brain.xlsx', Liver='Genes_Liver.xlsx', Muscle='Genes_Muscle.xlsx')

# The historical script references USARESSE_LIVER.xlsx for liver.
# The provided workbook is Genes_Liver.xlsx. This script uses the provided
# workbook; its identity to the historical alternate name was not verified.
if (dir.exists(output_dir) && length(list.files(output_dir)))
  stop('Output directory is not empty. Use a new output location.')
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

for (tissue in names(files)) {
  path <- file.path(input_dir,files[[tissue]])
  if (!file.exists(path)) stop('Missing candidate list: ',path)
  df <- readxl::read_excel(path)
  if (!'Genes' %in% names(df)) stop('Missing Genes column in ',path)
  genes <- df |>
    dplyr::select(Genes) |>
    tidyr::separate_rows(Genes, sep=',') |>
    dplyr::mutate(Genes=trimws(Genes)) |>
    dplyr::filter(!is.na(Genes), Genes!='') |>
    dplyr::distinct(Genes)
  write.table(genes$Genes,
              file=file.path(output_dir,paste0('genesunicos_',toupper(tissue),'_Multiomics.txt')),
              quote=FALSE,row.names=FALSE,col.names=FALSE)
  cat(tissue,': ',nrow(genes),' unique gene IDs\n',sep='')
}
