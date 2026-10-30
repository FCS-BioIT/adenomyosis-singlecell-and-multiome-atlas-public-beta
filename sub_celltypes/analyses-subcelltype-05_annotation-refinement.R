library(Seurat)
library(tidyverse)
library(ggplot2)
library(dittoSeq)
library(writexl)
library(patchwork)
library(harmony)
library(dplyr)
library(tidyr)
library(tidyverse)

input_dir <- "paper_output/zoom-in-mito/data/"
output_dir <- "paper_output/zoom-in-mito/"


celltypes <- c("epithelial", "pv", "endothelial", "smc", "immune", "stromal")
for (ct in celltypes) {
  dir.create(
    file.path(output_dir, "plots", "refined", ct, "annotations"),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# Ensure all necessary directories exist
celltypes <- c("endothelial", "epithelial", "smc", "pv", "immune", "stromal")

# Create all required subdirectories recursively
for (cell in celltypes) {
  dirs <- c(
    file.path(output_dir, "plots/", cell, "/re_annotations"),
    file.path(output_dir, "data/annotations/", cell)
  )
  sapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
}



gene_annotations <- read.csv(
  file = "paper_output/metadata/gene_annotations.csv",
  row.names = 1
)

sc_endo <- LoadSeuratRds(
  file = file.path(input_dir, "annotations/endothelial/sc_endo_annoted.Rds")
)
sc_epi <- LoadSeuratRds(
  file = file.path(input_dir, "annotations/epithelial/sc_epi_annoted.Rds")
)
sc_smc <- LoadSeuratRds(
  file = file.path(input_dir, "annotations/smc/sc_smc_annoted.Rds")
)
sc_pv <- LoadSeuratRds(
  file = file.path(input_dir, "annotations/pv/sc_pv_annoted.Rds")
)
sc_immune <- LoadSeuratRds(
  file = file.path(input_dir, "annotations/immune/sc_immune_annoted.Rds")
)
sc_stromal <- LoadSeuratRds(
  file = file.path(input_dir, "annotations/stromal/sc_stromal_annoted.Rds")
)


markers_endo <- c()
markers_epi <- c()
markers_smc <- c()
markers_pv <- c()
markers_immune <- c()
markers_stromal <- c()


df_endo <- read.csv(
  file.path(input_dir, "annotations/endothelial/endothelial_sub_celltypes.csv")
)
df_epi <- read.csv(
  file.path(input_dir, "annotations/epithelial/epithelial_sub_celltypes.csv")
)
df_smc <- read.csv(
  file.path(input_dir, "annotations/smc/smc_sub_celltypes.csv")
)
df_pv <- read.csv(
  file.path(input_dir, "annotations/pv/pv_sub_celltypes.csv")
)
df_immune <- read.csv(
  file.path(input_dir, "annotations/immune/immune_sub_celltypes.csv")
)
df_stromal <- read.csv(
  file.path(input_dir, "annotations/stromal/stromal_sub_celltypes.csv")
)

top5_endo <- df_endo %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)
top5_epi <- df_epi %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)
top5_smc <- df_smc %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)
top5_pv <- df_pv %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)
top5_immune <- df_immune %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)
top5_stromal <- df_stromal %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)


DimPlot(sc_endo, group.by = "sub_celltype") +
  DimPlot(sc_epi, group.by = "sub_celltype") +
  DimPlot(sc_smc, group.by = "sub_celltype") +
  DimPlot(sc_pv, group.by = "sub_celltype") +
  DimPlot(sc_immune, group.by = "sub_celltype") +
  DimPlot(sc_stromal, group.by = "sub_celltype")



unique(sc_smc@meta.data$lib_type)

# Provisional annotation for immune


c("CD3D", "CD3E", "CD4", "CD8A", "IL7R")
FeaturePlot(
  sc_immune,
  features = "nCount_RNA",
  min.cutoff = "q5",
  max.cutoff = "q95"
)

VlnPlot(
  sc_immune,
  features = c(
    "nFeature_RNA",
    "scDblFinder_RNA.score"
  ),
  group.by = "sub_celltype",
  pt.size = 0
)
DimPlot(sc_immune, group.by = "sub_celltype")

DimPlot(sc_immune, group.by = "RNA_snn_res.0.2")

DotPlot(
  sc_immune,
  features = c(
    "UCA1",
    "AC024230.1",
    "WFDC2",
    "EPCAM",
    "KLF5",
    "SDC4",
    "GABRP",
    "PAX8",
    "TMEM101",
    "ASRGL1",
    "ACSL5",
    "CD3D",
    "CD3E",
    "CD4",
    "CD8A",
    "IL7R"
  ),
  group.by = "sub_celltype"
) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Annotation pv
FeaturePlot(
  sc_pv,
  features = c("LAMA2", "PECAM1", "ELN", "TOP2A"),
  min.cutoff = "q5",
  max.cutoff = "q95"
) +
  DimPlot(sc_pv,
    label = TRUE,
    group.by = "sub_celltype"
  ) +
  DimPlot(
    sc_pv,
    group.by = "RNA_snn_res.0.2"
  )

# Reannotation of pv to gene defining clusters
sc_pv$sub_celltype_v2 <- sc_pv$sub_celltype

sc_pv$sub_celltype[sc_pv$RNA_snn_res.0.2 == "6"] <- "Pericytes / Contractile PV"

sc_pv$sub_celltype_v2[
  sc_pv$sub_celltype == "Contractile PV"
] <- "HMCN2+ PV"
sc_pv$sub_celltype_v2[
  sc_pv$sub_celltype == "Endothelial-like PV"
] <- "PECAM1+ PV"
sc_pv$sub_celltype_v2[
  sc_pv$sub_celltype == "Fibroblast-like PV"
] <- "EGFLAM+ PV"
sc_pv$sub_celltype_v2[
  sc_pv$sub_celltype == "Immune associated PV"
] <- "CCL2+ PV"
sc_pv$sub_celltype_v2[
  sc_pv$sub_celltype == "Nerve-vascular modulated PV"
] <- "DPP6+ PV"
sc_pv$sub_celltype_v2[
  sc_pv$sub_celltype == "Pericytes / Contractile PV"
] <- "ELN+ PV"
sc_pv$sub_celltype_v2[
  sc_pv$sub_celltype == "Stem-like PV (?)"
] <- "TOP2A+ PV"

DimPlot(sc_pv, group.by = "sub_celltype_v2")

SaveSeuratRds(
  sc_pv,
  file = file.path(output_dir, "data/annotations/pv/sc_pv_annoted.Rds")
)

Idents(sc_pv) <- "sub_celltype_v2"

future::plan("multisession", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
DefaultAssay(object = sc_pv) <- "RNA"
sc_pv_markers <- FindAllMarkers(object = sc_pv, assay = "RNA", min.pct = 0.3)
future::plan("sequential")

write.table(
  sc_pv_markers,
  file = file.path(
    output_dir,
    "data/annotations/pv/pv_sub_celltypes_refined.csv"
  ),
  sep = ",",
  col.names = TRUE,
  row.names = FALSE
)





# Annotation epi ---------------------------------------------------------------
FeaturePlot(
  sc_epi,
  features = c("PAEP", "VTCN1", "LGR5", "SUFU"),
  min.cutoff = "q5",
  max.cutoff = "q95"
) +
  DimPlot(sc_epi, label = TRUE, group.by = "sub_celltype") +
  DimPlot(
    sc_epi,
    label = TRUE,
    group.by = "RNA_snn_res.0.5"
  ) +
  DimPlot(
    sc_epi,
    label = TRUE,
    group.by = "RNA_snn_res.0.4"
  ) +
  patchwork::plot_layout(design = "ABCD
                                 GF#E")
# Reannotation to cluster at 0.5 resolution ------------------------------------

sc_epi$sub_celltype_v2 <- sc_epi$sub_celltype

sc_epi$sub_celltype_v2[
  which(sc_epi$RNA_snn_res.0.5 == "1" |
    sc_epi$RNA_snn_res.0.5 == "6" |
    sc_epi$RNA_snn_res.0.5 == "7")
] <- "Luminal-VTCN1"
sc_epi$sub_celltype_v2[
  which(sc_epi$RNA_snn_res.0.5 == "2" |
    sc_epi$RNA_snn_res.0.5 == "3" |
    sc_epi$RNA_snn_res.0.5 == "5")
] <- "Glandular"
sc_epi$sub_celltype_v2[
  which(sc_epi$RNA_snn_res.0.5 == "0" |
    sc_epi$RNA_snn_res.0.5 == "4" |
    sc_epi$RNA_snn_res.0.5 == "8" |
    sc_epi$RNA_snn_res.0.5 == "13")
] <- "Glandular secretory-PAEP"
sc_epi$sub_celltype_v2[
  which(sc_epi$RNA_snn_res.0.5 == "10")
] <- "Ciliated"
sc_epi$sub_celltype_v2[
  which(sc_epi$RNA_snn_res.0.5 == "11")
] <- "CTNNA2+"
sc_epi$sub_celltype_v2[
  which(sc_epi$RNA_snn_res.0.5 == "12")
] <- "Cycling cells"
sc_epi$sub_celltype_v2[
  which(sc_epi$RNA_snn_res.0.5 == "9")
] <- "Unidentified Doublets"

DimPlot(sc_epi, group.by = "sub_celltype_v2")

# pre filtered object save
SaveSeuratRds(
  sc_epi,
  file = file.path(
    output_dir,
    "data/annotations/epithelial/sc_epithelial_annoted_pre_filter.Rds"
  )
)

# Filtered object save
sc_epi_filt <- subset(
  sc_epi,
  subset = sub_celltype_v2 == "Unidentified Doublets",
  invert = TRUE
)

SaveSeuratRds(
  sc_epi_filt,
  file = file.path(
    output_dir,
    "data/annotations/epithelial/sc_epithelial_annoted.Rds"
  )
)


Idents(sc_epi) <- "sub_celltype_v2"

future::plan("multisession", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
DefaultAssay(object = sc_epi) <- "RNA"
sc_epi_markers <- FindAllMarkers(object = sc_epi, assay = "RNA", min.pct = 0.3)
future::plan("sequential")

write.table(
  sc_epi_markers,
  file = file.path(
    output_dir,
    "data/annotations/epithelial/epithelial_sub_celltypes_refined.csv"
  ),
  sep = ",",
  col.names = TRUE,
  row.names = FALSE
)


tbl <- table(sc_epi$sub_celltype_v2, sc_epi$group)
write.table(
  tbl,
  sep = "\t",
  quote = FALSE,
  col.names = NA,
  file = file.path(
    output_dir,
    "data/annotations/contingencies_sc_epi_refined.tsv"
  )
)


for (sc in ls(pattern = paste0("^sc_"))) {
  print(table(get(sc)$sub_celltype, get(sc)$group))
}

# --------------------------------------
for (sc in ls(pattern = "^sc_")) {
  cat("\n", sc, ":\n", sep = "")
  tbl <- table(get(sc)$sub_celltype, get(sc)$group)
  write.table(
    tbl,
    sep = "\t",
    quote = FALSE,
    col.names = NA,
    file = file.path(
      output_dir,
      "data/annotations",
      paste0(
        "contingencies_",
        sc,
        ".tsv"
      )
    )
  )
}

DimPlot(sc_epi, group.by = "sub_celltype_v2") +
  FeaturePlot(
    sc_epi,
    features = c("CDH2"),
    min.cutoff = "q5",
    max.cutoff = "q95"
  ) +
  FeaturePlot(
    sc_epi,
    features = c("LGR5"),
    min.cutoff = "q5",
    max.cutoff = "q95"
  ) +
  VlnPlot(sc_epi,
    features = c("CDH2", "LGR5"),
    group.by = "sub_celltype_v2", pt.size = 0
  ) +
  patchwork::plot_layout(design = "ABC\nDDD")




# REANNOTATION SMC -------------------------------------------------------------

FeaturePlot(
  sc_smc,
  features = c("ACTA2", "MYH11", "ACTG2"),
  min.cutoff = "q5",
  max.cutoff = "q95"
)

sc_smc$sub_celltype_v2 <- sc_smc$sub_celltype

DimPlot(
  sc_smc,
  group.by = "RNA_snn_res.0.2",
  label = TRUE
)
# DIFERENTIAL EXPRESSION BETWEEN CLUSTERS IN RES 0.2
contractile_markers_0 <- FindMarkers(
  sc_smc,
  ident.1 = c("0"),
  ident.2 = c("4", "6", "2"),
  group.by = "RNA_snn_res.0.2",
  min.pct = 0.5
)
contractile_markers_2 <- FindMarkers(
  sc_smc,
  ident.1 = c("2"),
  ident.2 = c("4", "6", "0"),
  group.by = "RNA_snn_res.0.2",
  min.pct = 0.5
)
contractile_markers_4 <- FindMarkers(
  sc_smc,
  ident.1 = c("4"),
  ident.2 = c("0", "6", "2"),
  group.by = "RNA_snn_res.0.2",
  min.pct = 0.5
)
contractile_markers_6 <- FindMarkers(
  sc_smc,
  ident.1 = c("6"),
  ident.2 = c("4", "0", "2"),
  group.by = "RNA_snn_res.0.2",
  min.pct = 0.5
)

contractile_markers_0_4 <- FindMarkers(
  sc_smc,
  ident.1 = c("0", "4"),
  ident.2 = c("6", "2"),
  group.by = "RNA_snn_res.0.2",
  min.pct = 0.5
)
contractile_markers_2_6 <- FindMarkers(
  sc_smc,
  ident.1 = c("2", "6"),
  ident.2 = c("0", "4"),
  group.by = "RNA_snn_res.0.2",
  min.pct = 0.5
)

contractile_markers_2_v_4 <- FindMarkers(
  sc_smc,
  ident.1 = c("2"),
  ident.2 = c("4"),
  group.by = "RNA_snn_res.0.2",
  min.pct = 0.3
)
contractile_markers_0_v_4 <- FindMarkers(
  sc_smc,
  ident.1 = c("0"),
  ident.2 = c("4"),
  group.by = "RNA_snn_res.0.2",
  min.pct = 0.3
)

contractile_markers_0_4_v_all <- FindMarkers(
  sc_smc,
  ident.1 = c("0", "4"),
  ident.2 = c("1", "2", "3", "5", "6", "7"),
  group.by = "RNA_snn_res.0.2",
  min.pct = 0.3
)
contractile_markers_0_4_v_2_6 <- FindMarkers(
  sc_smc,
  ident.1 = c("0", "4"),
  ident.2 = c("2", "6"),
  group.by = "RNA_snn_res.0.2",
  min.pct = 0.3
)

contractile_markers_0_v_4 <- FindMarkers(
  sc_smc,
  ident.1 = c("0"),
  ident.2 = c("4"),
  group.by = "RNA_snn_res.0.2",
  min.pct = 0.3
)

contractile_markers_0_4 <- merge(
  contractile_markers_0_4,
  gene_annotations,
  by = "row.names",
  all = FALSE
)
contractile_markers_2_6 <- merge(
  contractile_markers_2_6,
  gene_annotations,
  by = "row.names",
  all = FALSE
)
contractile_markers_2 <- merge(
  contractile_markers_2,
  gene_annotations,
  by = "row.names",
  all = FALSE
)
contractile_markers_6 <- merge(
  contractile_markers_6,
  gene_annotations,
  by = "row.names",
  all = FALSE
)

contractile_markers_2_v_4 <- merge(
  contractile_markers_2_v_4,
  gene_annotations,
  by = "row.names",
  all = FALSE
)

writexl::write_xlsx(
  contractile_markers_0_4_v_all %>% rownames_to_column(var = "gene"),
  path = file.path(
    output_dir,
    "data/annotations/smc/contractile_markers_0_4_v_all.xlsx"
  ),
  col_names = TRUE
)
writexl::write_xlsx(
  contractile_markers_0_4_v_2_6 %>% rownames_to_column(var = "gene"),
  path = file.path(
    output_dir,
    "data/annotations/smc/contractile_markers_0_4_v_2_6.xlsx"
  ),
  col_names = TRUE
)
writexl::write_xlsx(
  contractile_markers_0_v_4 %>% rownames_to_column(var = "gene"),
  path = file.path(
    output_dir,
    "data/annotations/smc/contractile_markers_0_v_4.xlsx"
  ),
  col_names = TRUE
)

write.csv(
  contractile_markers_0_4,
  file = file.path(output_dir, "data/annotations/smc/contractile_markers_0_4.csv")
)
write.csv(
  contractile_markers_2_6,
  file = file.path(output_dir, "data/annotations/smc/contractile_markers_2_6.csv")
)
write.csv(
  contractile_markers_0,
  file = file.path(output_dir, "data/annotations/smc/contractile_markers_0.csv")
)
write.csv(
  contractile_markers_2,
  file = file.path(output_dir, "data/annotations/smc/contractile_markers_2.csv")
)
write.csv(
  contractile_markers_4,
  file = file.path(output_dir, "data/annotations/smc/contractile_markers_4.csv")
)
write.csv(
  contractile_markers_6,
  file = file.path(output_dir, "data/annotations/smc/contractile_markers_6.csv")
)
write.csv(
  contractile_markers_2_v_4,
  file = file.path(
    output_dir,
    "data/annotations/smc/contractile_markers_2_v_4.csv"
  )
)

sc_smc$tissue_abbr <- sc_smc$tissue_code
sc_smc$tissue_abbr[sc_smc$tissue_abbr == "LM"] <- "LE"
sc_smc$group_tissue <- paste(sc_smc$group, sc_smc$tissue_abbr, sep = "")

sc_smc$group_tissue_aggr <- sc_smc$group_tissue

sc_smc$group_tissue_aggr[
  which(sc_smc$group_tissue == "COEN" | sc_smc$group_tissue == "COMI")
] <- "CO"
sc_smc$group_tissue_aggr[
  which(sc_smc$group_tissue == "ADEN" | sc_smc$group_tissue == "ADMI")
] <- "ADUT"
sc_smc$group_tissue_aggr[
  which(sc_smc$group_tissue == "ADLE")
] <- "ADLE"


dittoBarPlot(
  sc_smc,
  group.by = "RNA_snn_res.0.2", var = "group_tissue_aggr", scale = "percent"
) +
  dittoBarPlot(
    sc_smc,
    group.by = "RNA_snn_res.0.2", var = "group_tissue_aggr", scale = "count"
  ) +
  DimPlot(
    sc_smc,
    group.by = "RNA_snn_res.0.2", label = TRUE, split.by = "group_tissue_aggr"
  ) +
  patchwork::plot_layout(widths = c(1, 1, 2), design = "ABC")


sc_smc$sub_celltype_v2[
  which(sc_smc$RNA_snn_res.0.2 == 0)
] <- "Nerve-associated SMC"
sc_smc$sub_celltype_v2[
  which(sc_smc$RNA_snn_res.0.2 == 2)
] <- "Contractile SMC"
sc_smc$sub_celltype_v2[
  which(sc_smc$RNA_snn_res.0.2 == 3 | sc_smc$RNA_snn_res.0.2 == 5)
] <- "Stimuli-responsive SMC"
sc_smc$sub_celltype_v2[
  which(sc_smc$RNA_snn_res.0.2 == 4)
] <- "Proliferative SMC"
sc_smc$sub_celltype_v2[
  which(sc_smc$RNA_snn_res.0.2 == 1)
] <- "Canonical SMC"
sc_smc$sub_celltype_v2[
  which(sc_smc$RNA_snn_res.0.2 == 6)
] <- "Stress-responsive"
sc_smc$sub_celltype_v2[
  which(sc_smc$RNA_snn_res.0.2 == 7)
] <- "Unidentified"



# Markers
sc_smc <- SetIdent(sc_smc, value = "sub_celltype_v2")
future::plan("multisession", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
DefaultAssay(object = sc_smc) <- "RNA"
sc_smc_markers <- FindAllMarkers(object = sc_smc, assay = "RNA", min.pct = 0.3)
future::plan("sequential")

top5 <- sc_smc_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  select(gene)

DotPlot(
  sc_smc,
  group.by = "sub_celltype_v2",
  features = top5$gene
) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

sel <- c(
  "STAMBPL1", "CARMN", "KANSL1", "PBX1", "KCNMA1", "CACNA1C", "CACNA1D",
  "HOMER1", "LMCD1", "LMDC1-AS", "HRH1", "OSMR",
  "TSHZ2", "SGIP1", "PID1",
  "RYR3", "FRY",
  "CLRN1",
  "ALCAM",
  "RCAN1",
  "SERPINE1",
  "COBLL1",
  "CRIM1"
)
can <- c("ACTG2", "MYL9", "TPM2", "ELOB", "COX8A", "IFITM2")
DotPlot(
  sc_smc,
  group.by = "sub_celltype_v2",
  features = sel
) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

DotPlot(
  sc_smc,
  group.by = "RNA_snn_res.0.2",
  features = c(sel)
) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  DotPlot(
    sc_smc,
    group.by = "RNA_snn_res.0.2",
    features = c(can)
  ) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  DimPlot(
    sc_smc,
    group.by = "RNA_snn_res.0.2",
    label = TRUE,
    label.size = 6
  ) +
  patchwork::plot_layout(
    design = "AC\nBC",
    heights = c(3, 1),
    widths = c(1, 2)
  )

FeaturePlot(
  features = c("RYR3", "TSHZ2", "FRY", "FILIP1"),
  min.cutoff = "q5",
  max.cutoff = "q95",
  object = sc_smc
) +
  DimPlot(sc_smc, group.by = "RNA_snn_res.0.2", label = TRUE) +
  DimPlot(
    object = sc_smc,
    group.by = "sub_celltype_v2",
    label = TRUE
  ) +
  VlnPlot(
    sc_smc,
    features = c("RYR3", "TSHZ2", "FRY"),
    pt.size = 0,
    group.by = "RNA_snn_res.0.2"
  ) +
  patchwork::plot_layout(design = "ABEF\nCDGG", widths = c(1, 1, 1, 1))



# Save seurat of old annotation before filtering:
# SaveSeuratRds(object = sc_smc, file="output/zoom-in-mito/data/annotations/smc/sc_smc_annoted_pre_filter.Rds")

sc_smc_filt <- subset(
  sc_smc,
  subset = sub_celltype_v2 == "Unidentified",
  invert = TRUE
)
SaveSeuratRds(
  object = sc_smc_filt,
  file = file.path(output_dir, "data/annotations/smc/sc_smc_annoted.Rds")
)


# REANNOTATION FIBROBLASTS -----------------------------------------------------

sc_stromal$sub_celltype_v2 <- sc_stromal$sub_celltype
unique(sc_stromal$sub_celltype)
sc_stromal$sub_celltype_v2[
  which(sc_stromal$sub_celltype == "Decidualization fibroblasts (?)")
] <- "Decidualization fibroblasts"
sc_stromal$sub_celltype_v2[
  which(sc_stromal$sub_celltype == "WNT Modulated Fibroblasts" |
    sc_stromal$sub_celltype == "Response Fibroblasts")
] <- "Stimuli Response Fibroblasts"

# Markers
sc_stromal <- SetIdent(
  sc_stromal,
  value = "sub_celltype_v2"
)
future::plan("multisession", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
DefaultAssay(object = sc_stromal) <- "RNA"
sc_stromal_markers <- FindAllMarkers(
  object = sc_stromal,
  assay = "RNA",
  min.pct = 0.3
)
future::plan("sequential")

top5 <- sc_stromal_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  select(gene)


dittoSeq::dittoBarPlot(
  sc_stromal,
  group.by = "sub_celltype_v2",
  var = "tissue_code"
) +
  DimPlot(
    sc_stromal,
    group.by = "sub_celltype_v2"
  ) +
  DimPlot(
    sc_stromal,
    group.by = "tissue_code"
  ) +
  DotPlot(
    sc_stromal,
    group.by = "sub_celltype_v2",
    features = unique(top5$gene)
  ) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  patchwork::plot_layout(design = "ABD\nACD")

SaveSeuratRds(
  object = sc_stromal,
  file = file.path(
    output_dir, "data/annotations/stromal/sc_stromal_annoted.Rds"
  )
)


# REANNOTATION IMMUNE


# Markers res 0.4
sc_immune <- SetIdent(sc_immune, value = "RNA_snn_res.0.4")
future::plan("multisession", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
DefaultAssay(object = sc_immune) <- "RNA"
sc_immune_clust04_markers <- FindAllMarkers(
  object = sc_immune,
  assay = "RNA",
  min.pct = 0.3
)
future::plan("sequential")

top5 <- sc_immune_clust04_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  select(gene)

DimPlot(sc_immune, group.by = "RNA_snn_res.0.4", label = TRUE) +
  FeaturePlot(sc_immune, features = "mitoRatio") +
  DotPlot(
    sc_immune,
    group.by = "RNA_snn_res.0.4",
    features = unique(top5$gene)
  ) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  patchwork::plot_layout(design = "AC\nBC", heights = c(2, 1))

sc_immune$sub_celltype_v2 <- sc_immune$sub_celltype

sc_immune$sub_celltype_v2[
  which(sc_immune$sub_celltype == "Granulocytes")
] <- "Mast cells"
sc_immune$sub_celltype_v2[
  which(sc_immune$sub_celltype == "T/NK Cells")
] <- "T-NK cells"
sc_immune$sub_celltype_v2[
  which(sc_immune$sub_celltype == "Unidentified Doublets (?)")
] <- "Unidentified Doublets"
DimPlot(
  sc_immune,
  group.by = "sub_celltype_v2",
  label = TRUE
)

table(sc_immune$sub_celltype_v2)
# Take out cells from unidentified types


# The original file has been renamed as sc_immune_annoted_pre_filter.Rds
SaveSeuratRds(
  object = sc_immune,
  file = file.path(
    output_dir,
    "data/annotations/immune/sc_immune_annoted_pre_filter.Rds"
  )
)

sc_immune <- subset(
  sc_immune,
  subset = sub_celltype_v2 == "Unidentified Doublets" |
    sub_celltype_v2 == "Unidentified Non Conclusive",
  invert = TRUE
)

DimPlot(
  sc_immune,
  group.by = "sub_celltype_v2"
)

SaveSeuratRds(
  object = sc_immune,
  file = file.path(output_dir, "data/annotations/immune/sc_immune_annoted.Rds")
)

# REANOTATION ENDOTHELIUM
sc_endo
sc_endo$sub_celltype_v2 <- sc_endo$sub_celltype

SaveSeuratRds(
  object = sc_endo,
  file = file.path(
    output_dir,
    "data/annotations/endothelial/sc_endothelial_annoted.Rds"
  )
)


# REANNOTATION PLOTS -----------------------------------------------------------

# Define plots in function -----------------------------------------------------
ann_plots <- function(sc_obj, celltype) {
  db <- dittoBarPlot(
    sc_obj,
    var = "patient",
    group.by = "sub_celltype_v2",
    split.by = "group"
  ) +
    dittoBarPlot(
      sc_obj,
      var = "patient",
      group.by = "sub_celltype_v2",
      scale = "count",
      split.by = "group"
    )

  ggsave(
    plot = db,
    filename = file.path(
      output_dir,
      "plots",
      celltype,
      "re_annotations/dittobarplot.png"
    ),
    device = "png",
    width = 10,
    height = 6,
    dpi = 400,
    create.dir = TRUE
  )

  dimpl <- DimPlot(
    sc_obj,
    group.by = "sub_celltype_v2",
    split.by = "group",
    label = TRUE
  )
  ggsave(
    plot = dimpl,
    filename = file.path(
      output_dir,
      "plots",
      celltype,
      "re_annotations/dimplot-group.png"
    ),
    device = "png",
    width = 10,
    height = 6,
    dpi = 400,
    create.dir = TRUE
  )

  dimpl <- DimPlot(
    sc_obj,
    group.by = "sub_celltype_v2",
    split.by = "protocol",
    label = TRUE
  )
  ggsave(
    plot = dimpl,
    filename = file.path(
      output_dir,
      "plots",
      celltype,
      "re_annotations/dimplot-pcol.png"
    ),
    device = "png",
    width = 14,
    height = 6,
    dpi = 400,
    create.dir = TRUE
  )

  dimpl <- DimPlot(
    sc_obj,
    group.by = "sub_celltype_v2",
    label = TRUE
  )
  ggsave(
    plot = dimpl,
    filename = file.path(
      output_dir,
      "plots",
      celltype,
      "re_annotations/dimplot.png"
    ),
    device = "png",
    width = 8,
    height = 6,
    dpi = 400,
    create.dir = TRUE
  )

  # QC plots:
  dimpl <- FeaturePlot(
    sc_obj,
    features = c("nCount_RNA", "mitoRatio", "scDblFinder_RNA.score"),
    min.cutoff = "q5",
    max.cutoff = "q95",
    ncol = 3
  )
  ggsave(
    plot = dimpl,
    filename = file.path(
      output_dir,
      "plots",
      celltype,
      "re_annotations/qcplot.png"
    ),
    device = "png",
    width = 9,
    height = 3,
    dpi = 400,
    bg = "white",
    create.dir = TRUE
  )
  # Define genes and extract expression data
}

# DimPlots for refined subtypes:
# endothelial
ann_plots(sc_endo, celltype = "endothelial")
# smc
ann_plots(sc_smc_filt, celltype = "smc")
# immune
ann_plots(sc_immune, celltype = "immune")
# pv
ann_plots(sc_pv, celltype = "pv")
# epithelial
ann_plots(sc_epi_filt, celltype = "epithelial")
# stromal
ann_plots(sc_stromal, celltype = "stromal")
