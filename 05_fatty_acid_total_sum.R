# Fatty acid total-sum normalization — brain, liver, and muscle
# Procedure confirmed by authors to have been used for all three tissues.
# Original source script demonstrated this for muscle (AG_musculo).
# For each sample, columns 4 onward are fatty acid measurements.
# Columns 1-3 are metadata and carried over untouched.
# Derived FA sums/ratios in the 52-variable DIABLO block require independent
# documentation: do not compute ratios by normalizing precomputed ratios.
# The cleaned implementation is untested against historical tissue tables.
#
# Usage:
# Rscript preprocessing/05_fatty_acid_total_sum.R \
#   brain.csv liver.csv muscle.csv output_directory

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4L) {
  stop('Provide brain.csv liver.csv muscle.csv output_directory')
}
input_files <- setNames(args[1:3],c('brain','liver','muscle'))
out_dir <- args[[4]]
if (dir.exists(out_dir) && length(list.files(out_dir))) {
  stop('Choose an empty output directory; do not overwrite existing results.')
}
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

normalize_tissue <- function(path, tissue) {
  if (!file.exists(path)) stop('Missing input: ',path)
  df <- read.csv(path,check.names=FALSE,stringsAsFactors=FALSE)
  if (ncol(df) < 4L) stop(tissue,': expected sample ID, two metadata columns, fatty acids')
  vals <- df[,4:ncol(df),drop=FALSE]
  if (!all(vapply(vals,is.numeric,logical(1)))) {
    stop(tissue,': fatty-acid columns must be numeric')
  }
  sums <- rowSums(vals,na.rm=TRUE)
  if (any(!is.finite(sums) | sums<=0)) {
    stop(tissue,': missing or nonpositive total fatty-acid content')
  }
  normalized <- sweep(as.matrix(vals), 1L, sums, '/')*100
  result <- cbind(df[,1:3,drop=FALSE], as.data.frame(normalized,check.names=FALSE))
  result
}
for (tissue in names(input_files)) {
  data <- normalize_tissue(input_files[[tissue]], tissue)
  path <- file.path(out_dir,paste0('fatty_acids_normalized_',tissue,'.csv'))
  write.csv(data,path,row.names=FALSE)
  cat(tissue,': wrote ',nrow(data),' samples and ',ncol(data)-3L,' normalized FA columns\n',sep='')
}
