# Phenotypic QC: test normality and list variables eligible under the
# authors' selection criterion (Shapiro-Wilk p > 0.05).
# The final DIABLO PHENO matrix included author-selected variables retained
# without log/sqrt transformation; feature names were translated to English.
# This script DOES NOT reconstruct the translated final 'phenos.txt' because
# the variable-name mapping and original full phenotype table are not here.
# It documents the selection procedure and exports an inspection summary.
# The final 13-column matrix was independently verified against the RDS.
#
# Usage: Rscript preprocessing/06_phenotype_normality_QC.R raw_phenotypes.csv output_QC.csv
# Input CSV: one row per animal, optional metadata columns; numeric phenotypes.
# Edit metadata_fields to fit the source file.

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop('Provide raw_phenotypes.csv and output_QC.csv')
input_file <- args[[1]]; output_file <- args[[2]]
if (!file.exists(input_file)) stop('Missing phenotype input: ',input_file)
if (file.exists(output_file)) stop('Output exists; choose a fresh path.')
d <- read.csv(input_file,check.names=FALSE,stringsAsFactors=FALSE)
metadata_fields <- c('Amostra','tratamento','baia')
phenos <- d[,setdiff(names(d),metadata_fields),drop=FALSE]
if (!all(vapply(phenos,is.numeric,logical(1)))) {
  stop('Phenotype columns must be numeric; review metadata_fields.')
}
qc <- do.call(rbind,lapply(names(phenos),function(nm) {
  x <- phenos[[nm]]; x <- x[is.finite(x)]
  if (length(x)<3 || length(x)>5000 || length(unique(x))<3) {
    return(data.frame(Variavel=nm,W=NA_real_,p_value=NA_real_,Normal=NA_character_))
  }
  test <- shapiro.test(x)
  data.frame(Variavel=nm,W=unname(test$statistic),p_value=test$p.value,
             Normal=ifelse(test$p.value>0.05,'Sim','Não'))
}))
write.csv(qc,output_file,row.names=FALSE)
cat('Eligible by screening criterion (p > 0.05):\n')
print(qc$Variavel[!is.na(qc$p_value) & qc$p_value>0.05])
cat('Do not treat p > 0.05 as proof of normality.\n')
