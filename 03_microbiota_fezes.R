# Archived workflow, minimally adapted for portable directories.
# The scientific operations and original file-reading conventions are preserved.
# Starts from pre-existing ASV count, taxonomy, metadata tables: NOT DADA2.
# Not rerun/validated in an R environment as part of the repository preparation.

### 🚀 PIPELINE MICROBIOMA — FEZES E CECO ###
# 👉 Filtragem por prevalência + Agregação taxonômica + CLR com tratamento de zeros
# 👑 Camila Sabino 

### 1. Preparação ----------------------------------------------------

# Definir diretório
# Input/output paths are configurable and do not change the session working directory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript <this_script.R> <input_directory> <output_directory>")
input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Carregar bibliotecas
library(dplyr)
library(compositions) # Pacote recomendado para CLR, lida com zeros
library(zCompositions)  # <- Pacote usado para imputar zeros

# Nome dos arquivos
asv_file <- "asv_table_original_fezes.txt"
metadata_file <- "metadata_fezes.txt"
taxonomy_file <- "tax_table_original_fezes.txt"

# Carregar os dados
df_asv <- read.delim(file.path(input_dir, asv_file), sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
df_meta <- read.delim(file.path(input_dir, metadata_file), sep = "\t", header = TRUE, check.names = FALSE)
df_tax <- read.delim(file.path(input_dir, taxonomy_file), sep = " ", header = TRUE, row.names = 1, check.names = FALSE)

### 2. Filtragem por amostras ----------------------------------------

# Verificar se a coluna da metadata é 'sampleNames'
sample_ids <- unique(df_meta$sampleNames)
amostras_para_manter <- intersect(colnames(df_asv), sample_ids)
df_asv_filtrada <- df_asv[, amostras_para_manter]

### 3. Filtragem por prevalência -------------------------------------

# Definir prevalência mínima (ex.: 10%)
prevalencia_minima <- 0.10
min_amostras <- ceiling(prevalencia_minima * length(amostras_para_manter))
cat("Filtrando ASVs presentes em pelo menos", min_amostras, "amostras\n")
prevalencia <- rowSums(df_asv_filtrada > 0)
df_asv_prevalente <- df_asv_filtrada[prevalencia >= min_amostras, ]
cat("Número de ASVs após filtro de prevalência:", nrow(df_asv_prevalente), "\n")

### 4. Agregar por nível taxonômico ----------------------------------

# Definir nível taxonômico desejado: "Genus", "Family", "Phylum"...
nivel_taxonomico <- "Genus"  # <-- AJUSTAR SE QUISER OUTRO NÍVEL

if (!nivel_taxonomico %in% colnames(df_tax)) {
  stop(paste("A coluna", nivel_taxonomico, "não foi encontrada na tabela de taxonomia"))
}

# Extrair apenas a coluna taxonômica e adicionar ID da ASV
df_tax_com_asv <- df_tax[, nivel_taxonomico, drop = FALSE]
df_tax_com_asv$ASV_ID <- rownames(df_tax_com_asv)
colnames(df_tax_com_asv)[1] <- "Taxon"

# Remover entradas não informativas
df_tax_com_asv <- df_tax_com_asv[!(is.na(df_tax_com_asv$Taxon) | df_tax_com_asv$Taxon == "" | df_tax_com_asv$Taxon == "uncharacterized"), ]

# Adicionar coluna de ID no df_asv_prevalente
df_asv_com_id <- df_asv_prevalente
df_asv_com_id$ASV_ID <- rownames(df_asv_com_id)

# Juntar taxonomia + abundância
df_merged <- merge(df_tax_com_asv, df_asv_com_id, by = "ASV_ID")

# Agrupar por táxon e somar abundâncias
df_taxa_sum <- aggregate(. ~ Taxon, data = df_merged[,-1], FUN = sum)
rownames(df_taxa_sum) <- df_taxa_sum$Taxon
df_taxa_sum <- df_taxa_sum[, -1]  # Remove coluna "Taxon" duplicada

cat("Tabela no nível taxonômico", nivel_taxonomico, "gerada! Dimensões:", dim(df_taxa_sum), "\n")


### 5 e 6. Tratamento de Zeros com Imputação + Transformação CLR -----

# Transpor a matriz: amostras nas linhas, táxons nas colunas
df_taxa_sum_t <- t(df_taxa_sum)

# ✅ Substituir zeros usando método CZM (Count Zero Multiplicative)
df_taxa_sum_t_repl <- cmultRepl(df_taxa_sum_t, method = "CZM", output = "p-counts")

# Converter para objeto composicional
df_taxa_sum_t_acomp <- acomp(df_taxa_sum_t_repl)

# Aplicar CLR
df_taxa_clr_t <- clr(df_taxa_sum_t_acomp)

# Transpor de volta para o formato original
df_taxa_clr <- as.data.frame(t(df_taxa_clr_t))

cat("Matriz CLR no nível de", nivel_taxonomico, "gerada! Dimensões:", dim(df_taxa_clr), "\n")

### 7. Filtrar Metadata ----------------------------------------------

amostras_finais <- colnames(df_taxa_clr)
df_meta_filtrada <- df_meta %>%
  filter(sampleNames %in% amostras_finais) %>%
  arrange(match(sampleNames, amostras_finais))

# Verificar se a ordem bate
if(!all(df_meta_filtrada$sampleNames == colnames(df_taxa_clr))){
  stop("ERRO: A ordem dos metadados e da tabela CLR não bate!")
}

cat("Metadata filtrada e ordenada! Dimensões:", dim(df_meta_filtrada), "\n")

### 8. Salvar tudo ---------------------------------------------------

write.table(df_taxa_sum, file.path(output_dir, paste0("Contagens_", nivel_taxonomico, "_SUM.txt")), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(df_taxa_clr, file.path(output_dir, paste0("Contagens_", nivel_taxonomico, "_CLR.txt")), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(df_meta_filtrada, file.path(output_dir, "Metadata_filtrada_ordenada.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("🚀 Arquivos salvos com sucesso!\n")

##Se quiser conferir quais taxons foram excluídos, pois tinham soma total zero (ou quase zero) após a agregação.
taxa_antes <- colnames(df_taxa_sum_t)
taxa_depois <- colnames(df_taxa_clr_t)
taxa_excluidos <- setdiff(taxa_antes, taxa_depois)
length(taxa_excluidos) # diagnostic only; result depends on actual input
taxa_excluidos          # veja quem são

# --- 8b. EXPORTAR DE FORMA "EXCEL-SAFE" ---

# Garanta que é matriz numérica (evita exportar como texto estranho)
df_taxa_clr <- as.data.frame(df_taxa_clr)
stopifnot(all(sapply(df_taxa_clr, is.numeric)))

# (Opcional) arredondar só para reduzir tamanho de arquivo
# CLR não precisa de 15 casas decimais; 6–8 já é mais que suficiente
df_taxa_clr_out <- round(df_taxa_clr, 8)

# 1) CSV (melhor para o Excel importar sem quebrar notação científica)
write.csv(
  df_taxa_clr_out,
  file = file.path(output_dir, paste0("Contagens_", nivel_taxonomico, "_CLR_Fezes.csv")),
  row.names = TRUE,
  quote = FALSE
)

# 2) TSV "puro" (também funciona, mas o CSV costuma importar melhor)
write.table(
  df_taxa_clr_out,
  file = file.path(output_dir, paste0("Contagens_", nivel_taxonomico, "_CLR_Fezes.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

cat("✅ Exportação segura concluída (CSV e TSV).\n")

