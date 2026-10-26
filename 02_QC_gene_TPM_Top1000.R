# 02_QC_gene_TPM_Top1000.R — structured from QC_Genes_MultiOmics_Top1000.R
# STATUS: unvalidated recreation, requires original TPM matrices and BioMart exports.
# Maintains original: TPM > 1 in >=50% samples, three-tissue intersection,
# top 1000 by SD within the shared and non-shared sets per tissue.
# "Non-shared" means not in the 3-way intersection (NOT exclusively one tissue).
# Historical sample exclusions and annotation mapping are retained.
# Compare output feature names/dimensions with archived RDS before claiming exact reproduction.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
input_dir <- Sys.getenv('GENE_TPM_INPUT_DIR',unset='data/TPM')
list_dir <- Sys.getenv('GENE_LIST_DIR',unset='results/preprocessing/gene_lists')
mart_dir <- Sys.getenv('GENE_MART_DIR',unset='data/Ensembl_GeneNames')
output_dir <- Sys.getenv('GENE_TPM_OUTPUT_DIR',unset='results/preprocessing/TPM_filtered')
if (dir.exists(output_dir) && length(list.files(output_dir))) stop('Use an empty output directory.')
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

paths <- c(Brain='TPM_brain_final.txt', Liver='tpm_liver_expression.txt', Muscle='TPM_Muscle.txt')
mart_files <- c(Brain='mart_export_brain_1000.txt', Liver='mart_export_liver_1000.txt', Muscle='mart_export_muscle_1000.txt')
load_expr <- function(path) {
  if (!file.exists(path)) stop('Missing TPM matrix: ',path)
  x <- read.delim(path,header=TRUE,sep='\t',stringsAsFactors=FALSE,check.names=FALSE)
  names(x)[1] <- 'gene_id'
  x
}
load_list <- function(tissue) {
  path <- file.path(list_dir,paste0('genesunicos_',toupper(tissue),'_Multiomics.txt'))
  if (!file.exists(path)) stop('Missing list: ',path)
  readLines(path,warn=FALSE)
}
filter_tissue <- function(expr, candidates, min_expr=1, fraction=0.5) {
  expr <- expr[expr$gene_id %in% candidates,,drop=FALSE]
  vals <- expr[,-1,drop=FALSE]
  retained <- apply(vals,1,function(x) mean(x>min_expr,na.rm=TRUE)>=fraction)
  expr[retained,,drop=FALSE]
}
top_by_sd <- function(expr,genes,n=1000L) {
  x <- expr[expr$gene_id %in% genes,,drop=FALSE]
  if(!nrow(x)) return(x)
  sds <- apply(x[,-1,drop=FALSE],1,sd,na.rm=TRUE)
  # 'sd' is used only to rank; original file briefly retained it by accident.
  x[order(sds,decreasing=TRUE,na.last=TRUE)[seq_len(min(nrow(x),n))],,drop=FALSE]
}
expr <- lapply(names(paths),function(t)filter_tissue(load_expr(file.path(input_dir,paths[[t]])),load_list(t)))
names(expr) <- names(paths)
common <- Reduce(intersect,lapply(expr,function(x)x$gene_id))
selected <- lapply(expr,function(x) {
  rbind(top_by_sd(x,common,1000),top_by_sd(x,setdiff(x$gene_id,common),1000))
})

# Exact historical exclusion list (one duplicate entry in source script).
animals_exclude <- unique(c('B1C','B21C','B5C','B31C','B32C','B37C',
'B45C','B47C','B67C','B62C','B71C','B51C','B16C','B34C','B50C','B6C',
'B13C','B14C','B28C','B33C','B35C','B41C','B43C','B49C','B4C','B56C',
'B66C','B7C'))
rename_sample <- function(nm) {
  if (nm=='gene_id') return(nm)
  digits <- gsub('[^0-9]','',sub('^X','',nm))
  if (nzchar(digits)) paste0('B',digits,'C') else nm
}
for (t in names(selected)) {
  df <- selected[[t]]
  names(df) <- vapply(names(df),rename_sample,character(1))
  df <- df[,!(names(df) %in% animals_exclude),drop=FALSE]
  mart_path <- file.path(mart_dir,mart_files[[t]])
  if (!file.exists(mart_path)) stop('Missing BioMart file: ',mart_path)
  ann <- readr::read_tsv(mart_path,show_col_types=FALSE) |>
    dplyr::select(`Gene stable ID`,`Gene name`) |>
    dplyr::rename(gene_id=`Gene stable ID`, gene_name=`Gene name`) |>
    dplyr::distinct(gene_id,.keep_all=TRUE)
  df <- dplyr::left_join(df,ann,by='gene_id') |>
    dplyr::mutate(gene_final=ifelse(!is.na(gene_name) & gene_name!='',gene_name,gene_id)) |>
    dplyr::select(gene_final, dplyr::everything(),-gene_id,-gene_name)
  names(df)[1] <- 'gene_id'
  out_file <- file.path(output_dir,paste0(tolower(t),'_final.txt'))
  write.table(df,file=out_file,sep='\t',quote=FALSE,row.names=FALSE)
  cat(t,': ',nrow(df),' genes, ',ncol(df)-1,' sample columns\n',sep='')
}
