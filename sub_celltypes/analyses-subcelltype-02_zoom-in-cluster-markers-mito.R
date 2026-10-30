suppressPackageStartupMessages({
  library("Seurat")
  library("harmony")
  library("ggplot2")
  library("ggpubr")
  library("dittoSeq")
  library("dplyr")
  library("clustree")
  library("future")
  library("tidyverse")
  library("RhpcBLASctl")
})

Sys.setenv(
  OMP_NUM_THREADS = 24,
  MKL_NUM_THREADS = 24,
  OPENBLAS_NUM_THREADS = 24
)

RhpcBLASctl::blas_set_num_threads(24)

print(R.version)
print(.libPaths())

input_dir <- "paper_output/zoom-in-mito/"
output_dir <- "paper_output/zoom-in-mito/"


cell_list <- list(
  list(celltype = "pv", cell_nick = "pv"),
  list(celltype = "immune", cell_nick = "immune"),
  list(celltype = "endothelial", cell_nick = "endo"),
  list(celltype = "smc", cell_nick = "smc"),
  list(celltype = "epithelial", cell_nick = "epi"),
  list(celltype = "stromal", cell_nick = "stromal")
)

for (cell in cell_list) {
  celltype <- cell$celltype
  cell_nick <- cell$cell_nick
  ### RDSSeurat files to load
  sc_obj <- LoadSeuratRds(
    file.path(input_dir, "data", paste0("sc_", cell_nick, ".Rds"))
  )
  print("Seurat object on load: ")
  print(sc_obj)
  dir.create(file.path(output_dir, "plots", celltype), recursive = TRUE)
  dir.create(
    paste0(output_dir, "data/annotations/", celltype),
    recursive = TRUE
  )
  s_genes <- Seurat::cc.genes$s.genes
  g2m_genes <- Seurat::cc.genes$g2m.genes

  sc_obj <- Seurat::CellCycleScoring(
    sc_obj,
    s.features = s_genes,
    g2m.features = g2m_genes
  )

  sc_obj$cc_Phase <- sc_obj$Phase
  sc_obj$Phase <- sc_obj$S.Score - sc_obj$G2M.Score

  # we can observe how many cells have  net 0 mitoRatio
  # We calculate the cells which mitoRatio is abnormally large with isOutlier
  print(table(sc_obj$protocol))

  sc_obj <- NormalizeData(sc_obj)

  sc_obj <- FindVariableFeatures(object = sc_obj, nfeatures = 2000)

  print("Seurat object on found variable features: ")
  print(sc_obj)

  DefaultAssay(sc_obj) <- "RNA"

  sc_obj <- DietSeurat(sc_obj)

  sc_obj <- sc_obj %>%
    ScaleData(
      verbose = TRUE,
      vars.to.regress = c("nFeature_RNA", "mitoRatio", "Phase")
    ) %>%
    # 30 dimensions PCA
    # PCA to apply Harmony integration
    RunPCA(
      pc.genes = sc_obj@var.genes,
      npcs = 30,
      verbose = TRUE,
      seed.use = 42
    )

  print("Table of cells by protocol, group, patient, tissue: ")
  print(table(sc_obj$protocol, sc_obj$group, sc_obj$patient, sc_obj$tissue))


  sc_obj <- RunHarmony(
    object = sc_obj,
    group.by.vars = c("bioinformatics_id"),
    theta = 2,
    plot_convergence = TRUE,
    reduction.save = "harmony_sample_theta_2",
    assay.use = "RNA"
  )
  gc()

  sc_obj <- sc_obj %>%
    RunUMAP(
      reduction.name = "UMAP_sample_theta_2",
      reduction = "harmony_sample_theta_2",
      dims = 1:10,
      seed.use = 42
    ) %>%
    FindNeighbors(
      reduction = "harmony_sample_theta_2",
      dims = 1:10
    ) %>%
    FindClusters(
      resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
      random.seed = 42
    ) %>%
    identity()


  DefaultAssay(object = sc_obj) <- "RNA"

  clustree_out <- clustree(
    sc_obj,
    prefix = "RNA_snn_res.",
    exprs = c("scale.data"),
    assay = "RNA",
    show_axis = TRUE
  )

  clustree_out


  # CLUSTREE RESOLUTION--------------------------------------------------------
  resolution_ct <- "RNA_snn_res.0.2"
  #----------------------------------------------------------------------------
  # we settle at 0.2 resolution
  DimPlot(
    sc_obj,
    reduction = "UMAP_sample_theta_2",
    group.by = resolution_ct,
    pt.size = .8
  )
  if (celltype == "pv") {
    ggsave(file.path(output_dir, "plots", celltype, "test_real_env.png"))
  }

  sc_obj$scDblFinder_RNA.score

  DimPlot(
    sc_obj,
    reduction = "UMAP_sample_theta_2",
    group.by = resolution_ct,
    label = TRUE
  ) +
    DimPlot(
      sc_obj,
      reduction = "UMAP_sample_theta_2",
      group.by = "celltype",
      label = TRUE
    )
  dittoBarPlot(
    sc_obj,
    var = "bioinformatics_id",
    group.by = resolution_ct,
    split.by = "group"
  )

  DimPlot(sc_obj, group.by = "protocol") +
    dittoBarPlot(sc_obj, var = "protocol", group.by = resolution_ct) +
    FeaturePlot(
      sc_obj,
      features = c("scDblFinder_RNA.score", "mitoRatio")
    ) + DimPlot(sc_obj, group.by = resolution_ct, label = TRUE)

  clustering_stats <- ggarrange(
    dittoBarPlot(
      sc_obj,
      var = resolution_ct,
      group.by = "patient",
      scale = "count",
      legend.show = FALSE
    ),
    dittoBarPlot(
      sc_obj,
      var = resolution_ct,
      group.by = "patient",
      scale = "percent",
      legend.show = FALSE
    ),
    dittoBarPlot(
      sc_obj,
      var = resolution_ct,
      group.by = "group",
      scale = "count",
      legend.show = FALSE
    ),
    dittoBarPlot(
      sc_obj,
      var = resolution_ct,
      group.by = "group",
      scale = "percent"
    )
  )

  ggsave(
    filename = file.path(
      output_dir,
      "plots",
      celltype,
      "clustering_stats.png"
    ),
    device = "png",
    bg = "white",
    dpi = 400
  )

  table(sc_obj$patient)


  sc_obj <- SetIdent(sc_obj, value = resolution_ct)
  DimPlot(sc_obj)
  # Now we will find the markers for the clusters at resolution 0.3
  future::plan("multisession", workers = 8)
  options(future.globals.maxSize = 4000 * 1024^2)
  DefaultAssay(object = sc_obj) <- "RNA"
  sc_harmony_markers <- FindAllMarkers(
    object = sc_obj,
    assay = "RNA",
    min.pct = 0.3
  )
  future::plan("sequential")


  # Obtain Z-score from RNA data
  # zscore for each gene corresponds to mean divided by sd along all cells
  zscore_sc <- function(seurat_obj) {
    counts <- DelayedArray::DelayedArray(
      GetAssayData(seurat_obj, assay = "RNA", layer = "data")
    )
    counts <- (counts - DelayedMatrixStats::rowMeans2(counts)) / DelayedMatrixStats::rowSds(counts) # nolint
    counts[base::is.nan(counts)] <- 0
    counts
  }

  top5 <- sc_harmony_markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)


  ggsave(
    plot = clustree_out,
    file = file.path(
      output_dir,
      "plots",
      celltype,
      "clustree_out_harmony_sample.png"
    ),
    device = "png",
    width = 10,
    height = 15
  )

  htm <- DoHeatmap(
    object = sc_obj,
    features = top5$gene,
    group.by = resolution_ct,
    assay = "RNA",
    slot = "scale.data"
  )

  ggsave(
    filename = file.path(
      output_dir,
      "plots",
      celltype,
      "heatmap_cluster_markers.png"
    ),
    plot = htm,
    device = "png",
    width = 6,
    height = 8,
    bg = "white"
  )

  dp <- DotPlot(
    sc_obj,
    features = unique(top5$gene),
    group.by = resolution_ct,
    cols = "RdYlBu"
  ) +
    coord_flip()
  ggsave(
    filename = file.path(
      output_dir,
      "plots",
      celltype,
      "dotplot_cluster_markers.png"
    ),
    plot = dp,
    device = "png",
    width = 6,
    height = 8,
    bg = "white"
  )

  dimplot <- DimPlot(
    sc_obj,
    group.by = resolution_ct,
    reduction = "UMAP_sample_theta_2",
    label = TRUE
  ) + ggtitle(paste0(celltype, " ", resolution_ct))
  ggsave(
    filename = file.path(
      output_dir,
      "plots",
      celltype,
      "umap_cluster_markers.png"
    ),
    plot = dimplot,
    device = "png",
    width = 6,
    height = 6,
    bg = "white"
  )

  dp + dimplot
  ggsave(
    filename = file.path(
      output_dir,
      "plots",
      celltype,
      "dotplot_umap_cluster_markers.png"
    ),
    plot = dp + dimplot,
    device = "png",
    width = 12,
    height = 7,
    bg = "white"
  )

  SaveSeuratRds(
    object = sc_obj,
    file = file.path(
      output_dir,
      "/data/annotations/",
      celltype,
      paste0(
        "/sc_",
        cell_nick,
        "_clustered.Rds"
      )
    )
  )

  write.table(
    sc_harmony_markers,
    file = file.path(
      output_dir,
      "data/annotations",
      celltype,
      paste0(
        celltype,
        "_",
        resolution_ct,
        ".csv"
      )
    ),
    sep = ",",
    col.names = TRUE,
    row.names = FALSE
  )

  ### QC plots ----------------------------------------------------------------
  qc_plot <- FeaturePlot(
    sc_obj,
    features = c("mitoRatio", "scDblFinder_RNA.score"),
    reduction = "UMAP_sample_theta_2",
    label = TRUE
  )

  ggsave(
    filename = file.path(output_dir, "plots", celltype, "umap_qc.png"),
    plot = qc_plot,
    device = "png",
    width = 6,
    height = 6,
    bg = "white"
  )
}
