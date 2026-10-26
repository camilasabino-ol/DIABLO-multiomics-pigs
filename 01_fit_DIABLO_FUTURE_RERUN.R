# DIABLO multi-omics study — documented fitting workflow (DRAFT)
# Derived from script_DIABLO_PhenoFA_0.7_NewPlots_ColorCorrect(2).R
# IMPORTANT: not yet executed/validated; do not replace the original analysis.
# The original random seed / CV splits are unknown. A new seed below makes FUTURE
# reruns controlled but does NOT reconstruct the original cross-validation.

suppressPackageStartupMessages(library(mixOmics))
suppressPackageStartupMessages(library(BiocParallel))

input_dir <- Sys.getenv('DIABLO_INPUT_DIR', unset = 'data')
output_dir <- Sys.getenv('DIABLO_OUTPUT_DIR', unset = 'results')
if (normalizePath(input_dir, mustWork=FALSE) == normalizePath(output_dir, mustWork=FALSE)) {
  stop('Input and output directories must differ')
}
if (file.exists(file.path(output_dir, 'diablo_model_NEW_RUN.rds'))) {
  stop('New-run output already exists; choose a fresh DIABLO_OUTPUT_DIR to avoid overwriting')
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

read_block <- function(filename) {
  path <- file.path(input_dir, filename)
  if (!file.exists(path)) stop('Missing input file: ', path)
  x <- read.table(path, sep = '\t', header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE)
  x <- x[, sort(colnames(x)), drop = FALSE]
  x <- t(as.matrix(x))
  storage.mode(x) <- 'double'
  x
}

# Files are expected as variables x animals, with animal IDs in column names.
data <- list(
  Genes_MU = read_block('musculo.txt'),
  Genes_LI = read_block('figado.txt'),
  Genes_BR = read_block('cerebro.txt'),
  ASVs_CE  = read_block('ceco.txt'),
  ASVs_FE  = read_block('fezes.txt'),
  FA       = read_block('fatty_acids.txt'),
  PHENO    = read_block('phenos.txt')
)

groups_path <- file.path(input_dir, 'groups.txt')
if (!file.exists(groups_path)) stop('Missing group file: ', groups_path)
group <- read.table(groups_path, sep = '\t', header = TRUE,
                    stringsAsFactors = FALSE)
if (!all(c('X', 'group') %in% names(group)))
  stop('groups.txt must contain columns X and group')
group <- group[order(group$X), , drop = FALSE]
raw_groups <- as.character(group$group)
if (!all(raw_groups %in% c('soja', 'peixe', 'canola')))
  stop('Unexpected dietary group label(s)')
Y <- factor(c(soja='Soy', peixe='Fish', canola='Canola')[raw_groups],
            levels = c('Soy', 'Fish', 'Canola'))

# Explicit checks before fitting; do not silently reorder samples beyond the
# original script's sorting. Compare sample IDs and RDS model before release.
ids <- rownames(data[[1]])
if (!all(vapply(data, function(x) identical(rownames(x), ids), logical(1))))
  stop('Sample identifiers/order differ across blocks')
if (!identical(ids, as.character(group$X)))
  stop('Group IDs and block sample IDs/order differ; inspect input files')
if (length(Y) != 44L || anyNA(Y)) stop('Expected 44 nonmissing groups')

expected_dimensions <- c(Genes_MU=616, Genes_LI=731, Genes_BR=608,
                         ASVs_CE=151, ASVs_FE=138, FA=52, PHENO=13)
actual_dimensions <- vapply(data, ncol, integer(1))
if (!identical(actual_dimensions, expected_dimensions))
  stop('Input dimensions differ from archived RDS. Observed: ',
       paste(names(actual_dimensions), actual_dimensions, collapse=', '))

# Design matrix: 0.7 between all biological blocks, 0 on the diagonal.
design <- matrix(0.7, nrow=length(data), ncol=length(data),
                 dimnames=list(names(data), names(data)))
diag(design) <- 0

# The original fitting script used 10-fold CV, 100 repeats for perf(),
# and 10-fold CV, 10 repeats for tune.block.splsda().
# This seed is NEW, only for prospective reproducible runs.
set.seed(20261022)
basic.diablo.model <- block.splsda(X=data, Y=Y, ncomp=5, design=design)
perf.diablo <- perf(basic.diablo.model, validation='Mfold', folds=10,
                    nrepeat=100)
saveRDS(perf.diablo, file.path(output_dir, 'perf_diablo.rds'))

# Original selection rule. Verify against the archived model: ncomp must be 5.
ncomp <- as.integer(perf.diablo$choice.ncomp$WeightedVote['Overall.BER','max.dist'])
if (is.na(ncomp)) stop('Could not extract ncomp from perf()')

# The CV split seed in the original historical run is unknown. This NEW seed
# documents a future reproducibility attempt; it is not a recovered original seed.
# Grid copied from the original script.
test.keepX <- list(
  Genes_MU=c(6,8,12,18,20), Genes_LI=c(6,8,12,18,20),
  Genes_BR=c(6,8,12,18,20), ASVs_CE=seq(4,15,3),
  ASVs_FE=seq(4,15,3), FA=c(4,6,8,10), PHENO=c(4,6,8,10)
)
set.seed(20261022)
tune.result <- tune.block.splsda(
  X=data, Y=Y, ncomp=ncomp, test.keepX=test.keepX, design=design,
  validation='Mfold', folds=10, nrepeat=10, dist='max.dist',
  progressBar=TRUE, BPPARAM=BiocParallel::SnowParam(workers=5)
)
saveRDS(tune.result, file.path(output_dir, 'tuning_diablo.rds'))
list.keepX <- tune.result$choice.keepX
final.diablo.model <- block.splsda(X=data, Y=Y, ncomp=ncomp,
                                   keepX=list.keepX, design=design)
saveRDS(final.diablo.model, file.path(output_dir, 'diablo_model_NEW_RUN.rds'))
writeLines(capture.output(sessionInfo()), file.path(output_dir,'sessionInfo_NEW_RUN.txt'))

# IMPORTANT: compare ncomp, keepX, and selected variables against the
# archived final_diablo_model.rds. Differences may arise from CV randomness.
