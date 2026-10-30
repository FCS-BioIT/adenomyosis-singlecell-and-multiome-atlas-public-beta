# LIBRARIES -------------------------------------------------------------------
library(WebGestaltR)
library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
library(openxlsx)
library(ggplot2)

# DIRECTORIES -----------------------------------------------------------------
input_dir <- "paper_output/zoom-in-mito/"
output_dir <- "paper_output/de-zoom-minus-rb-mt-nopat-no23-nosn/"

# INPUTS ----------------------------------------------------------------------
celltypes <- c("stromal", "endothelial", "epithelial", "pv", "immune", "smc")
de_results_path <- file.path(output_dir, "data", "DE", "DE_ADENO_EMA_EMCO_celltypes.Rds")
de_results_all <- readRDS(de_results_path)

ora_jobs <- list(
  list(contrast = "EMA_EMCO", tissue = NULL, de_key = "ema_emco"),
  list(contrast = "ADENO_EMCO", tissue = "MI", de_key = "adeno_emco_mi"),
  list(contrast = "ADENO_EMCO", tissue = "EN", de_key = "adeno_emco_en"),
  list(contrast = "ADENO_EMA", tissue = "MI", de_key = "adeno_ema_mi"),
  list(contrast = "ADENO_EMA", tissue = "EN", de_key = "adeno_ema_en")
)

analyses <- list(
  GO = "geneontology_Biological_Process_noRedundant",
  KEGG = "pathway_KEGG"
)

# ORA -------------------------------------------------------------------------
for (job in ora_jobs) {
  contrast <- job$contrast
  tissue <- job$tissue
  de_key <- job$de_key
  tissue_label <- if (is.null(tissue)) "all_tissues" else tissue
  conditions <- unlist(strsplit(contrast, "_"))

  if (!de_key %in% names(de_results_all)) {
    message("DE result ", de_key, " is not in ", de_results_path)
    next
  }

  de_res_all <- de_results_all[[de_key]]

  for (celltype in celltypes) {
    message("Running ORA lfc0: ", contrast, " ", tissue_label, " ", celltype)

    if (!celltype %in% names(de_res_all)) {
      message("Celltype ", celltype, " is not in the DE results")
      next
    }

    sc_obj <- LoadSeuratRds(
      file.path(
        input_dir,
        "data",
        "annotations",
        celltype,
        paste0("sc_", celltype, "_annoted.Rds")
      )
    )

    reference_genes <- rownames(sc_obj[["RNA"]])
    genes_mito <- grep("^MT-", reference_genes, value = TRUE)
    genes_ribo <- grep("^RPS|^RPL", reference_genes, value = TRUE)
    reference_genes <- setdiff(reference_genes, unique(c(genes_mito, genes_ribo)))

    de_sig <- de_res_all[[celltype]] %>%
      filter(
        p_val_adj < 0.05,
        !is.na(avg_log2FC),
        gene %in% reference_genes
      )

    clusters_markers <- unique(de_sig$diffexp)
    if (length(clusters_markers) == 0) {
      message("No significant DE genes for ", celltype, " in ", contrast, " ", tissue_label)
      next
    }

    for (analysis_type in names(analyses)) {
      enrich_db <- analyses[[analysis_type]]
      output_dir_local <- file.path(
        output_dir,
        "data",
        paste0("ORA_lfc0_", analysis_type),
        tissue_label,
        contrast
      )
      wb <- createWorkbook()
      has_worksheets <- FALSE

      for (direction in c("UP", "DOWN")) {
        condition_label <- ifelse(direction == "UP", conditions[1], conditions[2])
        ora_list <- list()

        for (cluster_name in clusters_markers) {
          if (direction == "UP") {
            genes_use <- de_sig %>%
              filter(diffexp == cluster_name, avg_log2FC > 0)
          } else {
            genes_use <- de_sig %>%
              filter(diffexp == cluster_name, avg_log2FC < 0)
          }

          if (nrow(genes_use) == 0) {
            message("No ", direction, " genes for cluster ", cluster_name)
            next
          }

          ora_results <- NULL
          tryCatch(
            {
              ora_results <- WebGestaltR::WebGestaltR(
                enrichMethod = "ORA",
                organism = "hsapiens",
                minNum = 10,
                maxNum = 300,
                isOutput = FALSE,
                fdrThr = 0.05,
                enrichDatabase = enrich_db,
                interestGene = as.vector(genes_use$gene),
                interestGeneType = "genesymbol",
                referenceGene = reference_genes,
                referenceGeneType = "genesymbol"
              )
            },
            error = function(e) {
              message("Error in cluster ", cluster_name, ": ", e$message)
            }
          )

          if (is.data.frame(ora_results) && nrow(ora_results) > 0) {
            ora_results$condition <- condition_label
            ora_results$cluster <- paste("cluster", cluster_name)
            ora_list[[cluster_name]] <- ora_results
          }
        }

        if (length(ora_list) == 0) {
          message("ORA list was empty for ", celltype, " ", condition_label, " ", analysis_type)
          next
        }

        dir.create(output_dir_local, recursive = TRUE, showWarnings = FALSE)

        rds_path <- file.path(output_dir_local, gsub(
          x = paste0(
            celltype,
            "_ORA_lfc0_cluster_list_",
            condition_label,
            "_",
            analysis_type,
            "_",
            tissue_label,
            ".RDS"
          ),
          pattern = " ",
          replacement = "_"
        ))
        saveRDS(ora_list, file = rds_path)

        ora_df <- bind_rows(ora_list, .id = "column_label") %>%
          filter(FDR < 0.05)

        csv_path <- file.path(output_dir_local, gsub(
          x = paste0(
            celltype,
            "_ORA_lfc0_dataframe_",
            condition_label,
            "_",
            analysis_type,
            "_",
            tissue_label,
            ".csv"
          ),
          pattern = " ",
          replacement = "_"
        ))
        write.csv(ora_df, file = csv_path, row.names = FALSE)

        addWorksheet(wb, sheetName = paste0("ORA_", condition_label))
        has_worksheets <- TRUE
        writeData(wb, sheet = paste0("ORA_", condition_label), ora_df)

        if (nrow(ora_df) == 0) {
          next
        }

        for (cluster_label in unique(ora_df$cluster)) {
          plot_df <- ora_df %>%
            filter(cluster == cluster_label)

          sig_data <- plot_df %>%
            slice_max(enrichmentRatio, n = 10)

          sub_celltype <- stringr::str_replace(cluster_label, "^cluster\\s+", "")

          plot <- sig_data %>%
            ggplot(aes(
              x = reorder(description, enrichmentRatio),
              y = enrichmentRatio,
              size = size,
              fill = FDR
            )) +
            geom_point(shape = 21, color = "black") +
            scale_size(range = c(5, 12)) +
            coord_flip() +
            scale_fill_distiller(palette = "Spectral", direction = 1) +
            theme_cowplot() +
            theme(
              axis.text.y = element_text(size = 14, face = "bold"),
              axis.text.x = element_text(size = 10),
              axis.title.y = element_blank()
            ) +
            xlab("Functional annotation") +
            ylab("Enrichment ratio") +
            scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) +
            ggtitle(
              paste(
                sub_celltype,
                "-",
                contrast,
                "-\n",
                condition_label,
                "-",
                analysis_type,
                "-",
                tissue_label
              )
            )

          plot_dir <- file.path(
            output_dir,
            "plots",
            paste0("ORA_lfc0_", analysis_type),
            tissue_label,
            contrast,
            celltype
          )
          dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

          filename_base <- gsub(
            x = paste0(
              "ORA_lfc0_dotplot_cluster_",
              sub_celltype,
              "_",
              condition_label,
              "_",
              analysis_type,
              "_",
              tissue_label
            ),
            pattern = " ",
            replacement = "_"
          )

          ggsave(
            file.path(plot_dir, paste0(filename_base, ".png")),
            plot = plot,
            width = 9,
            height = 2 + length(sig_data$enrichmentRatio) * 0.8,
            dpi = 400,
            bg = "white",
            scale = 1.2
          )
          ggsave(
            file.path(plot_dir, paste0(filename_base, ".svg")),
            plot = plot,
            width = 9,
            height = 2 + length(sig_data$enrichmentRatio) * 0.8,
            dpi = 400,
            bg = "white",
            scale = 1.2
          )
        }
      }

      excel_path <- file.path(
        output_dir_local,
        paste0(celltype, "_ORA_lfc0_results_", analysis_type, "_", tissue_label, ".xlsx")
      )
      if (has_worksheets) {
        saveWorkbook(wb, file = excel_path, overwrite = TRUE)
      }
    }
  }

  for (analysis_type in names(analyses)) {
    wb <- createWorkbook()
    has_worksheets <- FALSE

    for (celltype in celltypes) {
      for (condition_label in conditions) {
        file_path <- file.path(
          output_dir,
          "data",
          paste0("ORA_lfc0_", analysis_type),
          tissue_label,
          contrast,
          gsub(
            x = paste0(
              celltype,
              "_ORA_lfc0_dataframe_",
              condition_label,
              "_",
              analysis_type,
              "_",
              tissue_label,
              ".csv"
            ),
            pattern = " ",
            replacement = "_"
          )
        )

        if (!file.exists(file_path)) {
          next
        }

        df <- read.csv(file_path)
        if (nrow(df) == 0) {
          next
        }

        sheet_name <- paste0(celltype, "_", condition_label)
        addWorksheet(wb, sheetName = sheet_name)
        has_worksheets <- TRUE
        writeData(wb, sheet = sheet_name, df)
      }
    }

    out_path <- file.path(
      output_dir,
      "data",
      paste0(
        "ORA_lfc0_all_celltypes_",
        contrast,
        "_",
        analysis_type,
        "_",
        tissue_label,
        ".xlsx"
      )
    )
    if (has_worksheets) {
      saveWorkbook(wb, file = out_path, overwrite = TRUE)
      message("Saved Excel: ", out_path)
    } else {
      message("No ORA lfc0 CSV files found for ", contrast, " ", analysis_type, " ", tissue_label)
    }
  }
}
