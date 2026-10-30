# Annotation of zoom in from Adenomyosis project
# Step 1: Subset the different objects
# Author: Jaime Llera Oyola

suppressPackageStartupMessages({
  library("Seurat")
  library("RColorBrewer")
  library("ggpubr")
  library("grid")
  library("ggplot2")
  library("ggrepel")
  library("dittoSeq")
  library("tidyverse")
  library("dplyr")
  library("EnhancedVolcano")
  library("SingleCellExperiment")
  library("future")
})

set.seed(47)

input_dir <- "input/"
output_dir <- "paper_output/zoom-in-mito/data/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

### RDSSeurat files to load
sc_multiple <- LoadSeuratRds(
  file.path(input_dir, "integration.bioinformatics_id.v5_annotated.rds")
)


DimPlot(
  sc_multiple,
  reduction = "umap.harmony",
  group.by = "celltype"
)
ggsave("paper_output/zoom-in-mito/general_umap_v5.png")

sc_multiple

sc_epi <- subset(
  x = sc_multiple,
  subset = celltype == "Ciliated" |
    celltype == "Epithelial"
)
sc_immune <- subset(
  x = sc_multiple,
  subset = celltype == "Immune_lymphoid" |
    celltype == "Immune_myeloid" |
    celltype == "Mast_cells"
)
sc_endo <- subset(
  x = sc_multiple,
  subset = celltype == "Endothelial" |
    celltype == "Lymphatic_endothelial"
)
sc_pv <- subset(x = sc_multiple, subset = celltype == "Perivascular")
sc_smc <- subset(x = sc_multiple, subset = celltype == "SMCs")
sc_stromal <- subset(x = sc_multiple, subset = celltype == "Stromal")

SaveSeuratRds(object = sc_epi, file = file.path(output_dir, "sc_epi.Rds"))
SaveSeuratRds(object = sc_pv, file = file.path(output_dir, "sc_pv.Rds"))
SaveSeuratRds(object = sc_immune, file = file.path(output_dir, "sc_immune.Rds"))
SaveSeuratRds(object = sc_endo, file = file.path(output_dir, "sc_endo.Rds"))
SaveSeuratRds(object = sc_smc, file = file.path(output_dir, "sc_smc.Rds"))
SaveSeuratRds(
  object = sc_stromal,
  file = file.path(output_dir, "sc_stromal.Rds")
)
