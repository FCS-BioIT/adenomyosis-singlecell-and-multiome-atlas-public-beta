# Plot selected subpopulation markers by condition and tissue.
#
# Feature plots show individual cells on the integrated UMAP. Expression boxplots
# match the Shiny tool: each observation is one cell's log1p-transformed RNA count.

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(svglite)

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

# Exact labels in `sub_celltype_v2`, corresponding to the marker sets above.
subtype_labels <- list(
  Epithelium = c(
    Cycling = "Cycling cells",
    Ciliated = "Ciliated",
    CTNNA2 = "CTNNA2+",
    Glandular = "Glandular",
    Glandular_secretory = "Glandular secretory-PAEP",
    Luminal = "Luminal-VTCN1"
  ),
  Stromal = c(
    Myofibroblast = "Myofibroblasts",
    Active_fibroblast = "Active Fibroblasts",
    Stimuli_response_fibroblast = "Stimuli Response Fibroblasts",
    Stressed_fibroblast = "Stressed Fibroblasts",
    Decidualized_fibroblast = "Decidualization fibroblasts"
  ),
  SMC = c(
    Canonical = "Canonical SMC",
    Contractile = "Contractile SMC",
    Nerve_associated = "Nerve-associated SMC",
    Proliferative = "Proliferative SMC",
    Stimuli_responsive = "Stimuli-responsive SMC",
    Stress_responsive = "Stress-responsive"
  ),
  Endothelial = c(
    Lymphatic = "Lymphatic endothelial",
    Arterial = "Arterial",
    Venous = "Venous"
  ),
  Immune = c(
    B_cells = "B Cells",
    Macrophages = "Macrophages",
    Mast_cells = "Mast cells",
    T_NK_cells = "T-NK cells"
  )
)

data_dir <- "paper_output/zoom-in-mito/data/annotations_conditions"
output_dir <- "paper_output/zoom-in-mito/plots/subpopulation_markers_by_condition_tissue"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seurat_files <- c(
  Epithelium = "sc_epithelial_annoted_conditions.Rds",
  Stromal = "sc_stromal_annoted_conditions.Rds",
  SMC = "sc_smc_annoted_conditions.Rds",
  Endothelial = "sc_endothelial_annoted_conditions.Rds",
  Immune = "sc_immune_annoted_conditions.Rds"
)

tissues <- c(
  Endometrium = "condition_ema_emco_adeno_en",
  Myometrium = "condition_ema_emco_adeno_mi"
)
condition_levels <- c("ADENO", "EMA", "EMCO")
condition_colours <- c(ADENO = "#dc267f", EMA = "#ffb000", EMCO = "#648fff")

save_plot <- function(plot, stem, width, height) {
  ggsave(
    filename = paste0(stem, ".png"), plot = plot,
    width = width, height = height, dpi = 300, bg = "white"
  )
  ggsave(
    filename = paste0(stem, ".svg"), plot = plot,
    device = svglite::svglite, width = width, height = height,
    bg = "white", fix_text_size = FALSE
  )
}

subset_tissue <- function(object, condition_column) {
  metadata <- object[[]]
  keep <- !is.na(metadata[[condition_column]]) &
    metadata[[condition_column]] %in% condition_levels
  object <- subset(object, cells = rownames(metadata)[keep])
  object$plot_condition <- factor(
    as.character(object[[]][[condition_column]]),
    levels = condition_levels
  )
  object
}

cell_log1p_expression <- function(object, genes, target_subtype) {
  cell_expression <- FetchData(
    object,
    vars = c("sub_celltype_v2", "plot_condition", genes),
    layer = "counts",
    clean = "none"
  )

  cell_expression |>
    filter(sub_celltype_v2 == target_subtype) |>
    pivot_longer(
      cols = all_of(genes), names_to = "gene", values_to = "expression"
    ) |>
    filter(
      !is.na(sub_celltype_v2), !is.na(plot_condition)
    ) |>
    mutate(
      exprs_log1p = log1p(expression),
      contrast = factor(plot_condition, levels = condition_levels)
    )
}

plot_cell_boxplots <- function(
    cell_expression, target_subtype, genes, condition_column) {
  expression_data <- cell_expression |>
    filter(gene %in% genes)

  if (nrow(expression_data) == 0L) {
    warning("No cell-level data for ", target_subtype, call. = FALSE)
    return(NULL)
  }

  # This is deliberately the same boxplot recipe as the Shiny tool:
  # per-cell log1p RNA counts and ggplot2's standard 1.5 x IQR whiskers.
  ggplot(
    expression_data,
    aes(x = sub_celltype_v2, y = exprs_log1p, fill = contrast)
  ) +
    geom_boxplot(outlier.shape = NA, staplewidth = 0.6, linewidth = 0.4) +
    facet_wrap(vars(gene), scales = "free_y", nrow = 1) +
    scale_fill_manual(
      name = condition_column,
      values = condition_colours,
      breaks = condition_levels,
      labels = condition_levels
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    ggtitle(label = target_subtype)
}

marker_validation <- list()
cell_counts <- list()

for (major_celltype in names(seurat_files)) {
  message("Loading ", major_celltype)
  object <- LoadSeuratRds(file.path(data_dir, seurat_files[[major_celltype]]))
  DefaultAssay(object) <- "RNA"

  marker_sets <- markers_by_celltype[[major_celltype]]
  marker_validation[[major_celltype]] <- tibble(
    major_celltype = major_celltype,
    marker_set = rep(names(marker_sets), lengths(marker_sets)),
    gene = unlist(marker_sets, use.names = FALSE)
  ) |>
    mutate(present_in_object = gene %in% rownames(object))

  for (tissue in names(tissues)) {
    message("  ", tissue)
    tissue_object <- subset_tissue(object, tissues[[tissue]])
    featureplot_dir <- file.path(
      output_dir, "featureplots", tolower(major_celltype), tolower(tissue)
    )
    boxplot_dir <- file.path(
      output_dir, "boxplots", tolower(major_celltype), tolower(tissue)
    )
    dir.create(featureplot_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(boxplot_dir, recursive = TRUE, showWarnings = FALSE)

    cell_counts[[paste(major_celltype, tissue, sep = "_")]] <- tissue_object[[]] |>
      count(plot_condition, name = "n_cells") |>
      mutate(major_celltype = major_celltype, tissue = tissue, .before = 1)

    available_genes <- marker_validation[[major_celltype]] |>
      filter(present_in_object) |>
      pull(gene) |>
      unique()
    for (marker_set in names(marker_sets)) {
      genes <- marker_sets[[marker_set]]
      available_set_genes <- genes[genes %in% rownames(tissue_object)]
      if (length(available_set_genes) == 0L) {
        warning(
          "No requested genes are present for ", major_celltype, "/", marker_set,
          call. = FALSE
        )
        next
      }

      feature_plot <- FeaturePlot(
        tissue_object,
        features = available_set_genes,
        reduction = "UMAP_sample_theta_2",
        split.by = "plot_condition",
        keep.scale = "feature",
        # Do not raise true zeroes to a quantile cutoff. This is important for
        # sparse condition/tissue panels, where q05 can make every cell equal.
        min.cutoff = 0,
        max.cutoff = "q99",
        order = TRUE,
        pt.size = 0.35,
        ncol = 3
      ) +
        plot_annotation(title = paste(major_celltype, marker_set, "—", tissue)) &
        NoAxes() &
        theme(plot.title = element_text(hjust = 0.5))

      featureplot_stem <- file.path(featureplot_dir, marker_set)
      save_plot(
        feature_plot, featureplot_stem,
        width = 12,
        height = max(4, 3.8 * length(available_set_genes))
      )

      target_subtype <- subtype_labels[[major_celltype]][[marker_set]]
      expression_boxplot <- plot_cell_boxplots(
        cell_log1p_expression(
          tissue_object, available_set_genes, target_subtype
        ),
        target_subtype = target_subtype,
        genes = available_set_genes,
        condition_column = tissues[[tissue]]
      )
      if (!is.null(expression_boxplot)) {
        save_plot(
          expression_boxplot, file.path(boxplot_dir, marker_set),
          width = max(7, 3.2 * length(available_set_genes)), height = 5
        )
      }
    }
  }
}

bind_rows(marker_validation) |>
  write.table(
    file.path(output_dir, "marker_validation.tsv"),
    sep = "\t", row.names = FALSE,
    quote = FALSE
  )
bind_rows(cell_counts) |>
  write.table(
    file.path(output_dir, "cell_counts_by_tissue_condition.tsv"),
    sep = "\t",
    row.names = FALSE, quote = FALSE
  )
