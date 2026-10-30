# Publication-style marker dotplots for all annotated zoom-in populations.
#
# Marker sets were explicitly transcribed from
# `metadata/Markers_subpopulations.docx` (October 2026). Corrections confirmed
# by the project owner: CTNNA2, RARRES1, SCN7A, SOCS3, and IGFBP1. The document
# calls the endothelial Capillary group "VSMC"; this script deliberately keeps
# the current Seurat annotation label, "Capillary".

library(Seurat)
library(ggplot2)
library(svglite)

data_dir <- "paper_output/zoom-in-mito/data/annotations_conditions"
output_dir <- "paper_output/zoom-in-mito/plots/dotplots_paper"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Names must match `sub_celltype_v2` in the corresponding Seurat object.
marker_sets <- list(
  Epithelium = list(
    "Ciliated" = c("DNAH12", "CFAP299"),
    "CTNNA2+" = c("CTNNA2", "RARRES1", "KIF26B"),
    "Cycling cells" = c("TOP2A", "ASPM"),
    "Glandular" = c("SUFU", "ADAMTS8"),
    "Glandular secretory-PAEP" = c("PAEP", "DPP4"),
    "Luminal-VTCN1" = c("LAMC2", "ENPP3")
  ),
  Stromal = list(
    "Myofibroblasts" = c("MYH11", "KCNMA1"),
    "Active Fibroblasts" = c("SCN7A", "EPHA3"),
    "Stimuli Response Fibroblasts" = c("ADGRL3", "NKAIN3"),
    "Stressed Fibroblasts" = c("GADD45B", "JUNB", "SOCS3"),
    "Decidualization fibroblasts" = c("IGFBP1", "LEFTY2", "KDR")
  ),
  SMC = list(
    "Canonical SMC" = c("COX8A", "UQCR10"),
    "Contractile SMC" = c("STAMBPL1", "RAP2C"),
    "Nerve-associated SMC" = c("RYR3", "LINGO2"),
    "Proliferative SMC" = c("PID1", "TSHZ2"),
    "Stimuli-responsive SMC" = c("HRH1", "HOMER1"),
    "Stress-responsive" = c("SERPINE1", "CLRN1")
  ),
  Endothelial = list(
    "Lymphatic endothelial" = c("CCL21", "PKHD1L1"),
    "Arterial" = c("FBLN5", "BMX"),
    "Capillary" = c("MYH11", "RGS5"),
    "Venous" = c("ACKR1", "TLL1")
  ),
  Immune = list(
    "B Cells" = c("MS4A1", "CD79A"),
    "Macrophages" = c("CD163", "MS4A6A"),
    "Mast cells" = c("TPSAB1", "TPSB2"),
    "T-NK cells" = c("CD247", "CD96")
  ),
  Perivascular = list(
    "DPP6+ PV" = c("DPP6", "LINGO2"),
    "PECAM1+ PV" = c("PECAM1", "EMCN"),
    "TOP2A+ PV" = c("TOP2A", "MKI67"),
    "CCL2+ PV" = c("CCL2", "STEAP4"),
    "ELN+ PV" = c("RERGL", "ELN"),
    "EGFLAM+ PV" = c("EGFLAM", "GUCY1A2"),
    "HMCN2+ PV" = c("HMCN2", "MYH11")
  )
)

# Marker order on the x axis. This is separate from the row order above so the
# epithelial panel matches the requested publication-style layout.
marker_plot_order <- lapply(marker_sets, names)
marker_plot_order$Epithelium <- c(
  "Cycling cells", "Ciliated", "CTNNA2+", "Glandular",
  "Glandular secretory-PAEP", "Luminal-VTCN1"
)

seurat_files <- c(
  Epithelium = "sc_epithelial_annoted_conditions.Rds",
  Stromal = "sc_stromal_annoted_conditions.Rds",
  SMC = "sc_smc_annoted_conditions.Rds",
  Endothelial = "sc_endothelial_annoted_conditions.Rds",
  Immune = "sc_immune_annoted_conditions.Rds",
  Perivascular = "sc_pv_annoted_conditions.Rds"
)

# Optional command-line cell type(s), useful when refreshing a single panel.
requested_celltypes <- commandArgs(trailingOnly = TRUE)
if (length(requested_celltypes) > 0L) {
  unknown_celltypes <- setdiff(requested_celltypes, names(seurat_files))
  if (length(unknown_celltypes) > 0L) {
    stop("Unknown cell type(s): ", paste(unknown_celltypes, collapse = ", "))
  }
  seurat_files <- seurat_files[requested_celltypes]
}

# Display labels make epithelial rows read naturally while retaining the Seurat
# labels above as the marker-set keys.
display_labels <- c(
  "Ciliated" = "Ciliated epithelium",
  "CTNNA2+" = "CTNNA2+ epithelium",
  "Cycling cells" = "Cycling cells epithelium",
  "Glandular" = "Glandular epithelium",
  "Glandular secretory-PAEP" = "Glandular secretory epithelium",
  "Luminal-VTCN1" = "Luminal epithelium"
)

marker_table <- do.call(
  rbind,
  lapply(names(marker_sets), function(celltype) {
    sets <- marker_sets[[celltype]]
    data.frame(
      major_celltype = celltype,
      sub_celltype_v2 = rep(names(sets), lengths(sets)),
      marker = unlist(sets, use.names = FALSE),
      row.names = NULL,
      check.names = FALSE
    )
  })
)
write.csv(
  marker_table,
  file.path(output_dir, "marker_sets.csv"),
  row.names = FALSE,
  quote = TRUE
)

# Grey-to-purple/blue mirrors the reference publication dotplot.  Passing
# these directly to `DotPlot()` keeps Seurat's standard narrow colourbar.
expression_colours <- c("#D9D9D9", "#3E50A4")

for (celltype in names(seurat_files)) {
  message("Plotting: ", celltype)
  object <- LoadSeuratRds(file.path(data_dir, seurat_files[[celltype]]))
  DefaultAssay(object) <- "RNA"

  requested_sets <- marker_sets[[celltype]]
  available_sets <- lapply(
    requested_sets,
    function(genes) genes[genes %in% rownames(object)]
  )
  available_sets <- available_sets[lengths(available_sets) > 0]

  if (length(available_sets) == 0L) {
    warning("No requested markers are present for ", celltype, call. = FALSE)
    next
  }

  missing_markers <- setdiff(unlist(requested_sets, use.names = FALSE), rownames(object))
  if (length(missing_markers) > 0L) {
    warning(
      celltype, ": omitted markers absent from the object: ",
      paste(missing_markers, collapse = ", "),
      call. = FALSE
    )
  }

  object$dotplot_subtype <- as.character(object$sub_celltype_v2)
  subtype_order <- names(requested_sets)
  object$dotplot_subtype <- factor(
    object$dotplot_subtype,
    levels = subtype_order,
    labels = unname(ifelse(
      subtype_order %in% names(display_labels),
      display_labels[subtype_order],
      subtype_order
    ))
  )

  features <- unname(unlist(available_sets[marker_plot_order[[celltype]]]))
  n_features <- length(features)
  # Compact dimensions keep the dot matrix dense, as in the reference panel.
  # Text is deliberately larger so the axes remain readable at figure scale.
  plot_width <- max(8, 2.4 + 0.58 * n_features)
  # Keep the marker labels at 45 degrees, while tightening the vertical
  # subtype matrix so rows sit closer together.
  plot_height <- max(3.2, 0.4 + 0.42 * length(available_sets))

  plot <- DotPlot(
    object = object,
    features = features,
    group.by = "dotplot_subtype",
    cols = expression_colours,
    dot.min = 0
  ) +
    scale_size(range = c(0.08, 6.5), name = "Percent Expressed") +
    RotatedAxis() +
    theme_classic(base_size = 13) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(size = 13, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "right",
      legend.key.height = grid::unit(0.30, "cm"),
      legend.spacing.y = grid::unit(0.10, "cm"),
      plot.margin = margin(t = 8, r = 8, b = 4, l = 4)
    )

  stem <- file.path(output_dir, paste0("DotPlot_", celltype))
  ggsave(
    paste0(stem, ".png"), plot,
    width = plot_width, height = plot_height, dpi = 400, bg = "white"
  )
  ggsave(
    paste0(stem, ".svg"), plot,
    device = svglite::svglite,
    width = plot_width, height = plot_height, bg = "white",
    fix_text_size = FALSE
  )
}
