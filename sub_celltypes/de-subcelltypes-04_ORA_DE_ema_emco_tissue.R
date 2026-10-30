# Over-representation analysis for one existing tissue-specific EMA versus EMCO
# DE result. Set EMA_EMCO_TISSUE to EN or MI. This script never recomputes DE.
# Run via the project conda environment: conda run -n adenomyosis Rscript <wrapper>.

suppressPackageStartupMessages({
  library(WebGestaltR)
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(cowplot)
  library(openxlsx)
  library(ggplot2)
})

tissue_code_use <- Sys.getenv("EMA_EMCO_TISSUE")
if (!tissue_code_use %in% c("EN", "MI")) {
  stop("EMA_EMCO_TISSUE must be either 'EN' or 'MI'.")
}
input_dir <- "paper_output/zoom-in-mito"
output_root <- "paper_output/de-zoom-minus-rb-mt-nopat-no23-nosn"
contrast_name <- "EMA_EMCO"
celltypes <- c("stromal", "endothelial", "epithelial", "pv", "immune", "smc")
de_outdir <- file.path(output_root, "data", "DE", tissue_code_use, contrast_name)
de_file <- file.path(de_outdir, paste0("DE_", contrast_name, "_", tissue_code_use, "_celltypes.Rds"))
if (!file.exists(de_file)) stop("Run tissue-specific DE first: ", de_file)
de_results <- readRDS(de_file)

analyses <- list(
  GO = "geneontology_Biological_Process_noRedundant",
  KEGG = "pathway_KEGG"
)
conditions <- c("EMA", "EMCO")

for (analysis_type in names(analyses)) {
  enrich_database <- analyses[[analysis_type]]
  ora_outdir <- file.path(
    output_root, "data", paste0("ORA_lfc0_", analysis_type), tissue_code_use, contrast_name
  )
  plot_outdir <- file.path(
    output_root, "plots", paste0("ORA_lfc0_", analysis_type), tissue_code_use, contrast_name
  )
  dir.create(ora_outdir, recursive = TRUE, showWarnings = FALSE)

  all_celltypes_workbook <- createWorkbook()
  all_celltypes_has_sheets <- FALSE

  for (celltype in celltypes) {
    if (!celltype %in% names(de_results)) {
      message("No ", celltype, " DE object in ", de_file)
      next
    }
    celltype_de <- de_results[[celltype]]
    if (!nrow(celltype_de)) {
      message("No DE rows for ", celltype, "; skipping ORA.")
      next
    }

    seurat <- LoadSeuratRds(file.path(
      input_dir, "data", "annotations", celltype, paste0("sc_", celltype, "_annoted.Rds")
    ))
    reference_genes <- rownames(seurat[["RNA"]])
    reference_genes <- setdiff(reference_genes, grep("^MT-|^RPS|^RPL", reference_genes, value = TRUE))
    significant_de <- celltype_de %>%
      filter(p_val_adj < 0.05, !is.na(avg_log2FC), gene %in% reference_genes)
    if (!nrow(significant_de)) {
      message("No significant DE genes for ", celltype, "; skipping ORA.")
      next
    }

    celltype_workbook <- createWorkbook()
    celltype_has_sheets <- FALSE
    for (direction in c("UP", "DOWN")) {
      condition_label <- if (direction == "UP") conditions[1] else conditions[2]
      ora_by_subtype <- list()

      for (subtype_contrast in unique(significant_de$diffexp)) {
        genes_use <- significant_de %>% filter(
          diffexp == subtype_contrast,
          if (direction == "UP") avg_log2FC > 0 else avg_log2FC < 0
        )
        if (!nrow(genes_use)) next

        ora_result <- tryCatch(
          WebGestaltR::WebGestaltR(
            enrichMethod = "ORA",
            organism = "hsapiens",
            minNum = 10,
            maxNum = 300,
            isOutput = FALSE,
            fdrThr = 0.05,
            enrichDatabase = enrich_database,
            interestGene = genes_use$gene,
            interestGeneType = "genesymbol",
            referenceGene = reference_genes,
            referenceGeneType = "genesymbol"
          ),
          error = function(error) {
            message("ORA failed for ", celltype, " / ", subtype_contrast, ": ", error$message)
            NULL
          }
        )
        if (is.data.frame(ora_result) && nrow(ora_result)) {
          ora_result$condition <- condition_label
          ora_result$cluster <- paste("cluster", subtype_contrast)
          ora_by_subtype[[subtype_contrast]] <- ora_result
        }
      }

      if (!length(ora_by_subtype)) {
        message("No ", direction, " ORA results for ", celltype, " / ", analysis_type)
        next
      }
      saveRDS(
        ora_by_subtype,
        file.path(
          ora_outdir,
          paste0(celltype, "_ORA_lfc0_cluster_list_", condition_label, "_", enrich_database, "_", tissue_code_use, ".RDS")
        )
      )
      ora_table <- bind_rows(ora_by_subtype, .id = "column_label") %>% filter(FDR < 0.05)
      csv_file <- file.path(
        ora_outdir,
        paste0(celltype, "_ORA_lfc0_dataframe_", condition_label, "_", enrich_database, "_", tissue_code_use, ".csv")
      )
      write.csv(ora_table, csv_file, row.names = FALSE)

      worksheet <- paste0("ORA_", condition_label)
      addWorksheet(celltype_workbook, worksheet)
      writeData(celltype_workbook, worksheet, ora_table)
      celltype_has_sheets <- TRUE
      addWorksheet(all_celltypes_workbook, paste0(celltype, "_", condition_label))
      writeData(all_celltypes_workbook, paste0(celltype, "_", condition_label), ora_table)
      all_celltypes_has_sheets <- TRUE

      if (!nrow(ora_table)) next
      for (cluster_label in unique(ora_table$cluster)) {
        plot_data <- ora_table %>%
          filter(cluster == cluster_label) %>%
          slice_max(enrichmentRatio, n = 10)
        subtype <- str_replace(cluster_label, "^cluster\\\\s+", "")
        ora_plot <- ggplot(plot_data, aes(
          x = reorder(description, enrichmentRatio), y = enrichmentRatio, size = size, fill = FDR
        )) +
          geom_point(shape = 21, colour = "black") +
          scale_size(range = c(5, 12)) +
          coord_flip() +
          scale_fill_distiller(palette = "Spectral", direction = 1) +
          theme_cowplot() +
          theme(axis.text.y = element_text(size = 14, face = "bold"), axis.title.y = element_blank()) +
          labs(
            title = paste(subtype, contrast_name, condition_label, analysis_type, tissue_code_use, sep = " - "),
            x = "Functional annotation", y = "Enrichment ratio"
          ) +
          scale_x_discrete(labels = function(x) str_wrap(x, width = 30))

        celltype_plot_dir <- file.path(plot_outdir, celltype)
        dir.create(celltype_plot_dir, recursive = TRUE, showWarnings = FALSE)
        filename <- paste0(
          "ORA_lfc0_dotplot_cluster_", str_replace_all(subtype, " ", "_"), "_", condition_label,
          "_", enrich_database, "_", tissue_code_use
        )
        for (extension in c("png", "svg")) {
          ggsave(
            file.path(celltype_plot_dir, paste0(filename, ".", extension)), ora_plot,
            width = 9, height = 2 + nrow(plot_data) * 0.8, dpi = 400, bg = "white", scale = 1.2
          )
        }
      }
    }

    if (celltype_has_sheets) {
      saveWorkbook(
        celltype_workbook,
        file.path(ora_outdir, paste0(celltype, "_ORA_lfc0_results_", enrich_database, "_", tissue_code_use, ".xlsx")),
        overwrite = TRUE
      )
    }
  }

  if (all_celltypes_has_sheets) {
    saveWorkbook(
      all_celltypes_workbook,
      file.path(
        output_root, "data",
        paste0("ORA_lfc0_all_celltypes_", contrast_name, "_", enrich_database, "_", tissue_code_use, ".xlsx")
      ),
      overwrite = TRUE
    )
  }
}

message("ORA completed for ", contrast_name, " in ", tissue_code_use)
