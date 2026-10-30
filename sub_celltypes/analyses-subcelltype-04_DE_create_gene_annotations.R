library(tidyverse)
library(openxlsx)

# Annotation of DEG

input_dir <- "metadata/"
output_dir <- "paper_output/metadata/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


#-----------------------------------------------------------------------------
## Load annotations

annotation_cellranger <- read.csv(
  file.path(input_dir, "features.tsv"),
  sep = "\t", header = FALSE
)
colnames(annotation_cellranger) <- c("GeneStableID", "GeneName", "Source")
annotation_cellranger <- annotation_cellranger[
  !duplicated(annotation_cellranger$GeneName),
]

annotation_ensembl <- read.csv(
  file.path(
    input_dir,
    "EnsGeneDescription.tsv"
  ),
  sep = "\t",
  header = TRUE
)
annotation_ensembl <- annotation_ensembl[
  !duplicated(annotation_ensembl$Gene.stable.ID),
]

annotation_ensembl_biotype <- read.csv(
  file.path(
    input_dir, "biotypes.csv"
  ),
  sep = ","
)

annotations_raw <- merge(
  annotation_cellranger,
  annotation_ensembl,
  by.x = "GeneStableID",
  by.y = "Gene.stable.ID",
  all.x = FALSE
)

annotations_raw <- merge(
  annotations_raw,
  annotation_ensembl_biotype,
  by.x = "GeneStableID",
  by.y = "Gene.stable.ID",
  all.x = FALSE
)


gene_annotations <- data.frame(
  row.names = annotations_raw$GeneName,
  gene_id = annotations_raw$GeneStableID,
  gene_description = annotations_raw$Gene.description,
  gene_type = annotations_raw$Gene.type
)
# Flag mitochondrial ribosomal
gene_annotations <- gene_annotations %>%
  mutate(gene_type = case_when(
    grepl(
      "^MT-",
      rownames(gene_annotations)
    ) ~ paste0("mitochondrial_", gene_type),
    grepl(
      "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",
      rownames(gene_annotations)
    ) &
      !grepl("lncRNA", gene_type) ~ paste0("ribosomal_", gene_type),
    TRUE ~ gene_type
  ))

write.csv(
  gene_annotations,
  file = paste0(output_dir, "gene_annotations.csv"),
  row.names = TRUE
)
#-----------------------------------------------------------------------------
