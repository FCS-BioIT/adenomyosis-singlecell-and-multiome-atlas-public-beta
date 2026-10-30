library(Seurat)
library(dittoSeq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(harmony)
library(tidyverse)

# Define input and output directories
input_dir <- "paper_output/zoom-in-mito/"
output_dir <- "paper_output/zoom-in-mito/"

# Ensure output folders exist --------------------------------------------------
celltypes <- c("endothelial", "immune", "pv", "stromal", "epithelial", "smc")
subdirs <- c("plots", "data")

for (ct in celltypes) {
  for (sd in subdirs) {
    dirs <- c(
      file.path(output_dir, sd, ct, "annotations")
    )
    for (d in dirs) {
      if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
  }
}

# Define plots in function -----------------------------------------------------
ann_plots <- function(sc_obj, genes, markers) {
  db <- dittoBarPlot(
    sc_obj,
    var = "patient",
    group.by = "sub_celltype",
    split.by = "group"
  ) + dittoBarPlot(
    sc_obj,
    var = "patient",
    group.by = "sub_celltype",
    scale = "count",
    split.by = "group"
  )
  ggsave(
    plot = db,
    filename = file.path(
      output_dir,
      "plots",
      celltype,
      "annotations/dittobarplot.png"
    ),
    device = "png",
    width = 10,
    height = 6,
    dpi = 400
  )

  dimpl <- DimPlot(
    sc_obj,
    group.by = "sub_celltype",
    split.by = "group",
    label = TRUE
  )
  ggsave(
    plot = dimpl,
    filename = file.path(
      output_dir,
      "plots",
      celltype,
      "annotations/dimplot-group.png"
    ),
    device = "png",
    width = 10,
    height = 6,
    dpi = 400
  )

  dimpl <- DimPlot(
    sc_obj,
    group.by = "sub_celltype",
    split.by = "protocol",
    label = TRUE
  )
  ggsave(
    plot = dimpl,
    filename = file.path(
      output_dir,
      "plots",
      celltype,
      "annotations/dimplot-pcol.png"
    ),
    device = "png",
    width = 14,
    height = 6,
    dpi = 400
  )

  dimpl <- DimPlot(
    sc_obj,
    group.by = "sub_celltype",
    label = TRUE
  )
  ggsave(
    plot = dimpl,
    filename = file.path(
      output_dir,
      "plots/",
      celltype,
      "/annotations/dimplot.png"
    ),
    device = "png",
    width = 8,
    height = 6,
    dpi = 400
  )

  # QC plots:
  dimpl <- FeaturePlot(
    sc_obj,
    features = c(
      "nCount_RNA",
      "mitoRatio",
      "scDblFinder_RNA.score"
    ),
    min.cutoff = "q5",
    max.cutoff = "q95",
    ncol = 3
  )
  ggsave(
    plot = dimpl,
    filename = file.path(
      output_dir,
      "plots/",
      celltype,
      "/annotations/qcplot.png"
    ),
    device = "png",
    width = 9,
    height = 3,
    dpi = 400,
    bg = "white"
  )
  # Define genes and extract expression data

  cell_types <- "sub_celltype" # Change to your actual cell type metadata

  # Extract expression data and format for ggplot
  expr_data <- FetchData(sc_obj, vars = genes) %>%
    rownames_to_column(var = "cell") %>%
    mutate(CellType = sc_obj@meta.data[cell, cell_types]) %>%
    pivot_longer(cols = genes, names_to = "Gene", values_to = "Expression")

  # Generate the faceted violin plot
  vlns <- ggplot(
    expr_data,
    aes(x = CellType, y = Expression, fill = CellType)
  ) +
    geom_violin(scale = "width", trim = TRUE) +
    facet_grid(rows = vars(Gene), cols = vars(CellType), scales = "free") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      y = "Expression Level",
      x = "Cell Type",
      title = "Violin Plot of Markers per Cell Type"
    ) +
    theme(panel.background = element_blank(), panel.grid = element_blank())

  ggsave(
    plot = vlns,
    filename = file.path(
      output_dir,
      "plots/",
      celltype,
      "/annotations/violins.png"
    ),
    device = "png",
    width = 8,
    height = length(genes) * 0.9,
    dpi = 400, bg = "white"
  )

  dp <- DotPlot(
    sc_obj,
    features = markers,
    cols = "RdYlBu"
  ) +
    coord_flip() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  ggsave(
    plot = dp,
    filename = file.path(
      output_dir,
      "plots",
      celltype,
      "annotations/dotplot.png"
    ),
    device = "png",
    dpi = 400,
    bg = "white",
    height = 0.4 * length(markers),
    width = 1.5 * length(unique(sc_obj$sub_celltype))
  )

  ggsave(
    filename = file.path(
      output_dir,
      "plots",
      celltype,
      "annotations/dotplot_umap_ct.png"
    ),
    plot = dimpl + dp,
    device = "png",
    width = 12,
    height = 7,
    bg = "white"
  )
}

## ENDOTHELIAL res 0.3
# 0: Venous CCL14, ACKR1
# 1: Venous NOSTRIN, LIFR, TLL1, LEF1
# 2: Venous: ACKR1 NOSTRIN
# 3: Arterial: NKAIN2, FN1, IGFBP3
# 4: Venous: LIFR, NOSTRIN
# 5: Pericytes: STEAP4, THBS1, ACTA2, MYH11 (less than 6)
# 6: Contractile capillary: TPM2, ACTA2, MYH11, MYL9
# 7: Lymphatic endothelial: CCL21, PDPN, COLEC12. Prev annotation

celltype <- "endothelial"
cell_nick <- "endo"
features <- c(
  "CCL14", "ACKR1",
  "NOSTRIN", "LIFR", "TLL1", "LEF1",
  "NKAIN2", "FN1", "IGFBP3",
  "STEAP4", "THBS1", "ACTA2", "MYH11",
  "TPM2", "MYL9",
  "CCL21", "PDPN", "COLEC12"
)
sc_obj <- LoadSeuratRds(
  paste0(
    input_dir,
    "/data/annotations/",
    celltype,
    "/sc_",
    cell_nick,
    "_clustered.Rds"
  )
)

sc_obj$sub_celltype <- sc_obj$celltype

sc_obj$sub_celltype[which(sc_obj$RNA_snn_res.0.2 == 0)] <- "Venous"
sc_obj$sub_celltype[which(sc_obj$RNA_snn_res.0.2 == 1)] <- "Venous"
sc_obj$sub_celltype[which(sc_obj$RNA_snn_res.0.2 == 2)] <- "Venous"
sc_obj$sub_celltype[which(sc_obj$RNA_snn_res.0.2 == 3)] <- "Capillary"
sc_obj$sub_celltype[which(sc_obj$RNA_snn_res.0.2 == 4)] <- "Arterial"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 5)
] <- "Lymphatic endothelial"



print("Compute endothelial markers")
# MARKER TABLES ----------------------------------------------------------------
sc_obj <- SetIdent(sc_obj, value = "sub_celltype")
# Now we will find the markers for the annotated celltypes
future::plan("multisession", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
DefaultAssay(object = sc_obj) <- "RNA"
sub_celltypes_markers <- FindAllMarkers(
  object = sc_obj,
  assay = "RNA",
  min.pct = 0.3
)
future::plan("sequential")

write.table(
  sub_celltypes_markers,
  file = file.path(
    output_dir,
    "data/annotations",
    celltype,
    paste0(celltype, "_sub_celltypes.csv")
  ),
  sep = ",",
  col.names = TRUE,
  row.names = FALSE
)

top5 <- sub_celltypes_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

ann_plots(sc_obj, features, markers = unique(top5$gene))

# SAVE OBJECT ------------------------------------------------------------------
SaveSeuratRds(
  sc_obj,
  file = file.path(
    output_dir,
    "data/annotations/",
    celltype,
    paste0(
      "/sc_",
      cell_nick,
      "_annoted.Rds"
    )
  )
)



## IMMUNE res 0.3 #TODO Put canonical markers here
# 0: Lymphoid (T/NK cells) CD3E NKG7 GZMB
# 1: Monocytes NKG7 LYZ FCGR3A
# 4: Granulocytes CSF3R MPO CEBPE
# 5: Unidentified (contd)
# 6: B Cells MS4A1 CD19 IGHM
# 7: Contaminant
# 8: Contaminant
# TODO: Put markers from immune
celltype <- "immune"
cell_nick <- "immune"
features <- c(
  "PTPRC", "NKG7",
  "CSF1R", "CD14",
  "KLRC1", "KIT",
  "IGKC", "IGHM", "MS4A1"
)
# features <- c(
#   "CD3E",   # T/NK Cells
#   "NKG7",   # T/NK Cells
#   "GZMB",   # T/NK Cells
#   "CD14",   # Monocytes
#   "LYZ",    # Monocytes
#   "FCGR3A", # Monocytes
#   "CSF3R",  # Granulocytes
#   "MPO",    # Granulocytes
#   "CEBPE",  # Granulocytes
#   "MS4A1",  # B Cells
#   "CD19",   # B Cells
#   "IGHM"    # B Cells
# )

sc_obj <- LoadSeuratRds(
  file.path(
    input_dir,
    "data/annotations/",
    celltype,
    paste0(
      "sc_",
      cell_nick,
      "_clustered.Rds"
    )
  )
)

sc_obj$sub_celltype <- sc_obj$celltype

sc_obj$sub_celltype[which(sc_obj$RNA_snn_res.0.2 == 0)] <- "T/NK Cells"
sc_obj$sub_celltype[which(sc_obj$RNA_snn_res.0.2 == 1)] <- "Macrophages"
sc_obj$sub_celltype[which(sc_obj$RNA_snn_res.0.2 == 2)] <- "T/NK Cells"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 3)
] <- "Unidentified Doublets (?)"
sc_obj$sub_celltype[which(sc_obj$RNA_snn_res.0.2 == 4)] <- "Granulocytes"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 5)
] <- "Unidentified Doublets (?)"
sc_obj$sub_celltype[which(sc_obj$RNA_snn_res.0.2 == 6)] <- "T/NK Cells"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 7)
] <- "Unidentified Non Conclusive"
sc_obj$sub_celltype[which(sc_obj$RNA_snn_res.0.2 == 8)] <- "B Cells"

DimPlot(sc_obj, group.by = "bioinformatics_id")

print("Compute immune markers")
# MARKER TABLES ----------------------------------------------------------------
sc_obj <- SetIdent(sc_obj, value = "sub_celltype")
# Now we will find the markers for the annotated celltypes
future::plan("multisession", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
DefaultAssay(object = sc_obj) <- "RNA"
sub_celltypes_markers <- FindAllMarkers(
  object = sc_obj,
  assay = "RNA",
  min.pct = 0.3
)
future::plan("sequential")

top5 <- sub_celltypes_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Plots
ann_plots(sc_obj, features, unique(top5$gene))

write.table(
  sub_celltypes_markers,
  file = file.path(
    output_dir,
    "data/annotations/",
    celltype,
    paste0(
      celltype,
      "_sub_celltypes.csv"
    )
  ),
  sep = ",",
  col.names = TRUE,
  row.names = FALSE
)

# SAVE OBJECT ------------------------------------------------------------------
SaveSeuratRds(
  sc_obj,
  file = file.path(
    output_dir,
    "data/annotations/",
    celltype,
    paste0(
      "sc_",
      cell_nick,
      "_annoted.Rds"
    )
  )
)



# PV res 0.3
# 0: ECM remodelling PV COL1A1 COL1A2 EGFLAM LAMB1 PLXDC1
# 1: Contractile PV ELN KCNAB1 CARML1 SYN3
# 2: Immune modulatory PV CCL2 IL1RAP C1R
# 3: Vascular homeostasis Regulatory PV
# (Vascular homeostasis KLF2 PHLDA RERGL) (MUSTN1 ADIRF mecanosensibility)
# 4: Vascular regulatory PV (endothelium homeostasis stress FOXC1 BCR)
# (ECM composition HS3ST2 GMDS)
# 5: Hormone response PV (PGR, FOXO1 GREB1)
# 6: Contractile PV DES ACTG2
# 7: Vascular PV? (RALGAPA2 EMCN FLT1)

celltype <- "pv"
cell_nick <- "pv"
features <- c(
  "COL1A1", "COL1A2", "EGFLAM", "LAMB1", "PLXDC1",
  "ELN", "CCL2", "IL1RAP", "C1R",
  "KLF2", "PGR", "FOXO1", "GREB1",
  "DES", "ACTG2",
  "RALGAPA2", "EMCN", "FLT1"
)

sc_obj <- LoadSeuratRds(
  file.path(
    input_dir,
    "data/annotations/",
    celltype,
    paste0(
      "sc_",
      cell_nick,
      "_clustered.Rds"
    )
  )
)

sc_obj$sub_celltype <- sc_obj$celltype

sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 0)
] <- "Pericytes / Contractile PV"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 1)
] <- "Fibroblast-like PV"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 2)
] <- "Immune associated PV"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 3)
] <- "Nerve-vascular modulated PV"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 4)
] <- "Contractile PV"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 5)
] <- "Endothelial-like PV"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 6)
] <- "Immune associated PV"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 7)
] <- "Stem-like PV (?)"

print("Compute PV markers")

# MARKER TABLES ----------------------------------------------------------------
sc_obj <- SetIdent(sc_obj, value = "sub_celltype")
# Now we will find the markers for the annotated celltypes
future::plan("multisession", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
DefaultAssay(object = sc_obj) <- "RNA"
sub_celltypes_markers <- FindAllMarkers(
  object = sc_obj,
  assay = "RNA",
  min.pct = 0.3
)
future::plan("sequential")

sub_celltypes_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) -> top5

# Plots
ann_plots(sc_obj, features, unique(top5$gene))


write.table(
  sub_celltypes_markers,
  file = file.path(
    output_dir,
    "data/annotations/",
    celltype,
    paste0(
      celltype,
      "_sub_celltypes.csv"
    )
  ),
  sep = ",",
  col.names = TRUE,
  row.names = FALSE
)

# SAVE OBJECT ----------------------------------------------------------------
SaveSeuratRds(
  sc_obj,
  file = file.path(
    output_dir,
    "data/annotations/",
    celltype,
    paste0(
      "sc_",
      cell_nick,
      "_annoted.Rds"
    )
  )
)



# STROMAL res 0.3
# 0: Active fibroblasts COL4A4, ITGA11, ADAMTS3
# 1: Responsive fibroblast: SDK1, EFNA5, IL15, SYNPO2, SIPA1L1
# 2: Responsive fibroblast: PRLR, RORB, ADAMTS9, TRPC4
# 3: Stressed fibroblasts: Ribosomal genes. RPS, RPL
# 4: Immune modulated fibroblasts: NFTC1,
#  NFKB1, ITGA5, ICAM1, ADAMTSL4-AS1, TGFB1
# 5: Myofibroblasts: CARMN
# 6: Myofibroblasts: MYH11


celltype <- "stromal"
cell_nick <- "stromal"
features <- c(
  "COL4A1", "ITGA11", "ADAMTS3",
  "SDK1", "EFNA5", "IL15", "SYNPO2", "SIPA1L1",
  "PRLR", "RORB", "ADAMTS9", "TRPC4",
  "RPS5", "RPL8",
  "NFATC1", "NFKB1", "ITGA5", "ICAM1", "ADAMTSL4-AS1", "TGFB1",
  "CARMN",
  "MYH11"
)

sc_obj <- LoadSeuratRds(
  file.path(
    input_dir,
    "data/annotations/",
    celltype,
    paste0(
      "sc_",
      cell_nick,
      "_clustered.Rds"
    )
  )
)

sc_obj$sub_celltype <- sc_obj$celltype

sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 0)
] <- "Response Fibroblasts"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 1)
] <- "Stressed Fibroblasts"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 2)
] <- "Active Fibroblasts"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 3)
] <- "WNT Modulated Fibroblasts"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 4)
] <- "Decidualization fibroblasts (?)"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 5)
] <- "Myofibroblasts"


print("Compute stromal markers")
# MARKER TABLES ----------------------------------------------------------------
sc_obj <- SetIdent(sc_obj, value = "sub_celltype")
# Now we will find the markers for the annotated celltypes
future::plan("multisession", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
DefaultAssay(object = sc_obj) <- "RNA"
sub_celltypes_markers <- FindAllMarkers(
  object = sc_obj,
  assay = "RNA",
  min.pct = 0.3
)
future::plan("sequential")

top5 <- sub_celltypes_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Plots
ann_plots(sc_obj, features, unique(top5$gene))


write.table(
  sub_celltypes_markers,
  file = file.path(
    output_dir,
    "data/annotations",
    celltype,
    paste0(
      celltype,
      "_sub_celltypes.csv"
    )
  ),
  sep = ",",
  col.names = TRUE,
  row.names = FALSE
)

# SAVE OBJECT ------------------------------------------------------------------
SaveSeuratRds(
  sc_obj,
  file = file.path(
    output_dir,
    "data/annotations/",
    celltype,
    paste0(
      "sc_",
      cell_nick,
      "_annoted.Rds"
    )
  )
)




# EPITHELIAL res 0.3
# 0: Glandular secretory PAEP DPP4 PTPRR
# 1: Glandular SCGB1D2 SCGB1D4
# 2: Luminal PTGS1
# 3: Luminal PDGFC
# 4: Pre Glandular secretory CEACAM SUFU PGR
# 5: Luminal PDGFC
# 6: Ciliated PIFO FOXJ1 TP73
# 7: OUT
# 8: Luminal cycling TOP2A MKI67 CDK1

celltype <- "epithelial"
cell_nick <- "epi"
features <- c(
  "PAEP", "DPP4", "PTPRR",
  "SCGB1D2", "SCGB1D4",
  "PTGS1", "PDGFC",
  "CEACAM1", "SUFU",
  "PIFO", "FOXJ1", "TP73",
  "TOP2A", "MKI67", "CDK1"
)
sc_obj <- LoadSeuratRds(
  file.path(
    input_dir,
    "data/annotations/",
    celltype,
    paste0(
      "sc_",
      cell_nick,
      "_clustered.Rds"
    )
  )
)

sc_obj$sub_celltype <- sc_obj$celltype

sc_obj$sub_celltype[which(sc_obj$RNA_snn_res.0.2 == 0)] <- "Glandular secretory"
sc_obj$sub_celltype[which(sc_obj$RNA_snn_res.0.2 == 1)] <- "Glandular"
sc_obj$sub_celltype[which(sc_obj$RNA_snn_res.0.2 == 2)] <- "Luminal"
sc_obj$sub_celltype[which(sc_obj$RNA_snn_res.0.2 == 3)] <- "Glandular"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 4)
] <- "Unidentified Doublets (?)" # CEACAM SUFU PGR
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 5)
] <- "Ciliated"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 6)
] <- "Unidentified CTNNA2+"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 7)
] <- "Cycling Cells"


DimPlot(sc_obj, group.by = "bioinformatics_id")

print("Compute epithelial markers")
# MARKER TABLES ----------------------------------------------------------------
sc_obj <- SetIdent(sc_obj, value = "sub_celltype")
# Now we will find the markers for the annotated celltypes
future::plan("multisession", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
DefaultAssay(object = sc_obj) <- "RNA"
sub_celltypes_markers <- FindAllMarkers(
  object = sc_obj,
  assay = "RNA",
  min.pct = 0.3
)
future::plan("sequential")


sub_celltypes_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) -> top5

# Plots
ann_plots(sc_obj, features, unique(top5$gene))

write.table(
  sub_celltypes_markers,
  file = file.path(
    output_dir,
    "data/annotations/",
    celltype,
    paste0(
      celltype,
      "_sub_celltypes.csv"
    )
  ),
  sep = ",",
  col.names = TRUE,
  row.names = FALSE
)

# SAVE OBJECT ------------------------------------------------------------------
SaveSeuratRds(
  sc_obj,
  file = file.path(
    output_dir,
    "data/annotations",
    celltype,
    paste0(
      "/sc_",
      cell_nick,
      "_annoted.Rds"
    )
  )
)



# SMC res 0.3
# 0: Contractility regulator SMC: RYR3, SGCD
# 1: Contractile SMC: CARMN, KCNMA1, KANSC1
# 2: ECM remodel SMC: FBN1, ADAMTS19, ESR1
# 3: Stimuli response SMC: HOMER1
# 4: Canonical SMC: TAGLN, S100A6
# 5: Canonical SMC: TAGLN, S100A6
# 6: Stimuli response SMC: HOMER1
# 7: Unidentified
# 8: OUT

celltype <- "smc"
cell_nick <- "smc"
features <- c(
  "RYR3", "SGCD",
  "CARMN", "KCNMA1", "KANSL1",
  "FBN1", "ADAMTS19", "ESR1",
  "HOMER1",
  "TAGLN", "S100A6"
)

sc_obj <- LoadSeuratRds(
  file.path(
    input_dir,
    "data/annotations/",
    celltype,
    paste0(
      "/sc_",
      cell_nick,
      "_clustered.Rds"
    )
  )
)

sc_obj$sub_celltype <- sc_obj$celltype

sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 0)
] <- "Vascular SMC"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 1)
] <- "Stimuli response SMC"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 2)
] <- "Stimuli response SMC"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 3)
] <- "Canonical SMC"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 4)
] <- "Canonical SMC"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 5)
] <- "Stimuli response SMC"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 6)
] <- "Canonical SMC"
sc_obj$sub_celltype[
  which(sc_obj$RNA_snn_res.0.2 == 7)
] <- "Unidentified Poor quality cells"

DimPlot(sc_obj, group.by = "bioinformatics_id")

print("Compute smc markers")

# MARKER TABLES ----------------------------------------------------------------
sc_obj <- SetIdent(sc_obj, value = "sub_celltype")
# Now we will find the markers for the annotated celltypes
future::plan("multisession", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
DefaultAssay(object = sc_obj) <- "RNA"
sub_celltypes_markers <- FindAllMarkers(
  object = sc_obj,
  assay = "RNA",
  min.pct = 0.3
)
future::plan("sequential")


top5 <- sub_celltypes_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Plots
ann_plots(sc_obj, features, unique(top5$gene))


write.table(
  sub_celltypes_markers,
  file = file.path(
    output_dir,
    "data/annotations/",
    celltype,
    paste0(
      celltype,
      "_sub_celltypes.csv"
    )
  ),
  sep = ",",
  col.names = TRUE,
  row.names = FALSE
)

# SAVE OBJECT ------------------------------------------------------------------
SaveSeuratRds(
  sc_obj,
  file = file.path(
    output_dir,
    "data/annotations",
    celltype,
    paste0(
      "sc_",
      cell_nick,
      "_annoted.Rds"
    )
  )
)
