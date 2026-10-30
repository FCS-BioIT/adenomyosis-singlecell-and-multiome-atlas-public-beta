library(Seurat)
library(ggplot2)

markers_by_celltype <- list(
  Epithelium = list(
    Cycling = c("TOP2A", "ASPM"),
    Ciliated = c("DNAH12", "CFAP299"),
    CTNNA2 = c("CTNNA2", "RARRES1", "KIF26B"),
    Glandular = c("SUFU", "ADAMTS8"),
    Glandular_secretory = c("PAEP", "DPP4"),
    Luminal = c("LAMC2", "ENPP3")
  ),
  Stromal = list(
    Myofibroblast = c("MYH11", "KCNMA1"),
    Active_fibroblast = c("SCN7A", "EPHA3"),
    Stimuli_response_fibroblast = c("ADGRL3", "NKAIN3"),
    Stressed_fibroblast = c("GADD45B", "JUNB", "SOCS3"),
    Decidualized_fibroblast = c("IGFBP1", "LEFTY2", "KDR")
  ),
  SMC = list(
    Canonical = c("COX8A", "UQCR10"),
    Contractile = c("STAMBPL1", "RAP2C"),
    Nerve_associated = c("RYR3", "LINGO2"),
    Proliferative = c("PID1", "TSHZ2"),
    Stimuli_responsive = c("HRH1", "HOMER1"),
    Stress_responsive = c("SERPINE1", "CLRN1")
  ),
  Endothelial = list(
    Lymphatic = c("CCL21", "PKHD1L1"),
    Arterial = c("FBLN5", "BMX"),
    Venous = c("ACKR1", "TLL1")
  ),
  Immune = list(
    B_cells = c("MS4A1", "CD79A"),
    Macrophages = c("CD163", "MS4A6A"),
    Mast_cells = c("TPSAB1", "TPSB2"),
    T_NK_cells = c("CD247", "CD96")
  )
)


# Input/output directories
data_dir <- "paper_output/zoom-in-mito/data/annotations_conditions"
plot_dir <- "paper_output/zoom-in-mito/plots/dotplots_markers_subtypes"

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Match cell types with corresponding Seurat objects
seurat_files <- c(
  Epithelium  = "sc_epithelial_annoted_conditions.Rds",
  Stromal     = "sc_stromal_annoted_conditions.Rds",
  SMC         = "sc_smc_annoted_conditions.Rds",
  Endothelial = "sc_endothelial_annoted_conditions.Rds",
  Immune      = "sc_immune_annoted_conditions.Rds"
)

plots <- list()


for (celltype in names(seurat_files)) {
  message("Plotting: ", celltype)

  # Load Seurat object
  obj <- LoadSeuratRds(
    file.path(
      data_dir,
      seurat_files[[celltype]]
    )
  )

  # Markers corresponding to this major cell type
  features <- markers_by_celltype[[celltype]]

  # Keep only genes present in the Seurat object
  features <- lapply(
    features,
    \(x) x[x %in% rownames(obj)]
  )

  # Remove subgroups for which none of the genes are present
  features <- features[
    lengths(features) > 0
  ]

  # Total number of genes actually represented
  n_features <- sum(lengths(features))


  # Dynamic plot width
  # Minimum width = 6 inches
  # Each additional gene contributes ~0.6 inches
  plot_width <- max(
    6,
    3 + (1 * n_features)
  )

  message(
    "  Features: ", n_features,
    " | Width: ", round(plot_width, 1), " inches"
  )


  # Generate DotPlot
  p <- DotPlot(
    object = obj,
    features = features,
    group.by = "sub_celltype_v2"
  ) +
    RotatedAxis() +
    ggtitle(celltype) +
    theme(
      plot.title = element_text(
        hjust = 0.5
      )
    )


  # Store plot in memory
  plots[[celltype]] <- p


  # Save PNG
  ggsave(
    filename = file.path(
      plot_dir,
      paste0(
        "DotPlot_",
        celltype,
        ".png"
      )
    ),
    plot = p,
    width = plot_width,
    height = 6,
    dpi = 300,
    bg = "white"
  )


  # Save SVG
  ggsave(
    filename = file.path(
      plot_dir,
      paste0(
        "DotPlot_",
        celltype,
        ".svg"
      )
    ),
    plot = p,
    width = plot_width,
    height = 6,
    bg = "white",
    fix_text_size = FALSE
  )

  rm(obj)
  gc()
}
