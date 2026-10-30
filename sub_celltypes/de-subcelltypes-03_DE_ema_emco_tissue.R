# Tissue-specific EMA versus EMCO differential expression.
#
# This is intentionally parameterised only by EMA_EMCO_TISSUE (EN or MI), so
# the two wrapper scripts can be run independently without recomputing the
# all-tissues or ADENO contrasts.
# Run via the project conda environment: conda run -n adenomyosis Rscript <wrapper>.

suppressPackageStartupMessages({
  library(Seurat)
  library(tibble)
  library(openxlsx)
  library(dplyr)
  library(purrr)
  library(future)
})

tissue_code_use <- Sys.getenv("EMA_EMCO_TISSUE")
if (!tissue_code_use %in% c("EN", "MI")) {
  stop("EMA_EMCO_TISSUE must be either 'EN' or 'MI'.")
}

input_dir <- "paper_output/zoom-in-mito"
output_root <- "paper_output/de-zoom-minus-rb-mt-nopat-no23-nosn"

contrast_name <- "EMA_EMCO"
de_outdir <- file.path(output_root, "data", "DE", tissue_code_use, contrast_name)
dir.create(de_outdir, recursive = TRUE, showWarnings = FALSE)

workers <- as.integer(Sys.getenv("DE_WORKERS", unset = "16"))
if (is.na(workers) || workers < 1L) stop("DE_WORKERS must be a positive integer.")
report_only <- tolower(Sys.getenv("DE_REPORT_ONLY", unset = "false")) %in% c("1", "true", "yes")
if (!report_only) plan(multisession, workers = workers)
options(future.globals.maxSize = 8000 * 1024^2)

updated_md <- openxlsx::read.xlsx("metadata/metadata_2026.xlsx")
stopifnot(!anyDuplicated(updated_md$bioinformatics_id))

celltypes <- c("pv", "immune", "endothelial", "smc", "epithelial", "stromal")

initialize_sample_report <- function(report_file) {
  dir.create(dirname(report_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(
    paste(
      c(
        "record_type", "main_celltype", "subcelltype", "contrast",
        "condition", "samples", "n_cells", "covariates", "status"
      ),
      collapse = "\t"
    ),
    report_file
  )
}

append_sample_report <- function(report_file, main_celltype, ident_1, ident_2,
                                 seurat, valid_latent_vars, status) {
  condition_labels <- sub("^.*_", "", c(ident_1, ident_2))
  subtype <- sub(paste0("_", condition_labels[1], "$"), "", ident_1)
  sample_rows <- lapply(c(ident_1, ident_2), function(ident) {
    cells <- Cells(seurat)[seurat$cluster_group == ident]
    samples <- sort(unique(as.character(seurat$bioinformatics_id[cells])))
    data.frame(
      record_type = "condition",
      main_celltype = main_celltype,
      subcelltype = subtype,
      contrast = paste(condition_labels, collapse = "_vs_"),
      condition = sub("^.*_", "", ident),
      samples = paste(samples, collapse = ";"),
      n_cells = length(cells),
      covariates = "",
      status = status,
      stringsAsFactors = FALSE
    )
  })
  covariate_row <- data.frame(
    record_type = "covariates",
    main_celltype = main_celltype,
    subcelltype = subtype,
    contrast = paste(condition_labels, collapse = "_vs_"),
    condition = "",
    samples = "",
    n_cells = sum(vapply(sample_rows, `[[`, integer(1), "n_cells")),
    covariates = paste(valid_latent_vars, collapse = ";"),
    status = status,
    stringsAsFactors = FALSE
  )
  write.table(
    do.call(rbind, c(sample_rows, list(covariate_row))),
    report_file,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    append = TRUE
  )
}

shared_subtypes <- function(seurat) {
  ema_idents <- paste(
    sort(unique(as.character(seurat$sub_celltype_v2[seurat$condition_ema_emco == "EMA"]))),
    "EMA",
    sep = "_"
  )
  emco_idents <- paste(
    sort(unique(as.character(seurat$sub_celltype_v2[seurat$condition_ema_emco == "EMCO"]))),
    "EMCO",
    sep = "_"
  )
  intersect(sub("_EMA$", "", ema_idents), sub("_EMCO$", "", emco_idents))
}

run_mast <- function(seurat, ident_1, ident_2, latent_vars, report_file, main_celltype) {
  cells_use <- Cells(seurat)[seurat$cluster_group %in% c(ident_1, ident_2)]
  covariates <- seurat@meta.data[cells_use, latent_vars, drop = FALSE]
  valid_latent_vars <- latent_vars[vapply(
    covariates,
    function(x) length(unique(stats::na.omit(x))) >= 2L,
    logical(1)
  )]
  n_cells_1 <- sum(seurat$cluster_group == ident_1)
  n_cells_2 <- sum(seurat$cluster_group == ident_2)
  report_status <- if (n_cells_1 >= 10L && n_cells_2 >= 10L) {
    "tested"
  } else {
    "skipped_fewer_than_10_cells"
  }
  append_sample_report(
    report_file = report_file,
    main_celltype = main_celltype,
    ident_1 = ident_1,
    ident_2 = ident_2,
    seurat = seurat,
    valid_latent_vars = valid_latent_vars,
    status = report_status
  )
  if (report_status != "tested") {
    message("Skipping ", ident_1, " vs ", ident_2, ": fewer than 10 cells in one group.")
    return(NULL)
  }
  message(
    ident_1, " vs ", ident_2, "; latent variables: ",
    if (length(valid_latent_vars)) paste(valid_latent_vars, collapse = ", ") else "none"
  )

  markers <- FindMarkers(
    seurat,
    ident.1 = ident_1,
    ident.2 = ident_2,
    assay = "RNA",
    min.pct = 0.25,
    test.use = "MAST",
    latent.vars = valid_latent_vars,
    random.seed = 1
  )
  markers$diffexp <- paste0(ident_1, "-vs-", ident_2)
  rownames_to_column(markers, var = "gene")
}

prepare_celltype <- function(celltype) {
  seurat <- LoadSeuratRds(file.path(
    input_dir, "data", "annotations", celltype, paste0("sc_", celltype, "_annoted.Rds")
  ))

  metadata_index <- match(seurat$bioinformatics_id, updated_md$bioinformatics_id)
  if (anyNA(metadata_index)) {
    stop(
      "Missing metadata for: ",
      paste(unique(seurat$bioinformatics_id[is.na(metadata_index)]), collapse = ", ")
    )
  }
  metadata_columns <- setdiff(colnames(updated_md), "bioinformatics_id")
  metadata_add <- updated_md[metadata_index, metadata_columns, drop = FALSE]
  rownames(metadata_add) <- Cells(seurat)
  seurat@meta.data[, metadata_columns] <- metadata_add

  # Final-DE exclusions: no snGEX samples and no patient 23 cells.
  # Mitochondrial/ribosomal features are removed before running MAST.
  seurat <- subset(
    seurat,
    subset = lib_type %in% c("scGEX", "snMU") & patient != 23
  )
  genes <- rownames(seurat[["RNA"]])
  genes_keep <- setdiff(genes, grep("^MT-|^RPS|^RPL", genes, value = TRUE))
  seurat <- subset(seurat, features = genes_keep)

  seurat$condition_ema_emco <- dplyr::case_when(
    seurat$group == "AD" ~ "EMA",
    seurat$group == "CO" ~ "EMCO",
    TRUE ~ NA_character_
  )
  seurat <- subset(
    seurat,
    subset = tissue_code == tissue_code_use & condition_ema_emco %in% c("EMA", "EMCO")
  )
  seurat$condition_ema_emco <- factor(seurat$condition_ema_emco, levels = c("EMCO", "EMA"))
  seurat$cluster_group <- paste(seurat$sub_celltype_v2, seurat$condition_ema_emco, sep = "_")
  Idents(seurat) <- "cluster_group"
  seurat
}

sample_report_file <- function(celltype) {
  file.path(
    de_outdir,
    "sample_reports",
    paste0("samples_covariates_", contrast_name, "_", tissue_code_use, "_", celltype, ".txt")
  )
}

write_celltype_sample_report <- function(celltype) {
  message("Writing sample report for ", contrast_name, " in ", tissue_code_use, ": ", celltype)
  seurat <- prepare_celltype(celltype)
  report_file <- sample_report_file(celltype)
  initialize_sample_report(report_file)
  subtypes <- shared_subtypes(seurat)
  if (!length(subtypes)) {
    message("No shared subcelltypes in ", celltype, " for ", tissue_code_use, " EMA vs EMCO.")
    return(invisible(NULL))
  }
  walk(subtypes, function(subtype) {
    ident_1 <- paste0(subtype, "_EMA")
    ident_2 <- paste0(subtype, "_EMCO")
    cells_use <- Cells(seurat)[seurat$cluster_group %in% c(ident_1, ident_2)]
    latent_vars <- c("nCount_RNA", "current_horm_treatm", "lib_type")
    covariates <- seurat@meta.data[cells_use, latent_vars, drop = FALSE]
    valid_latent_vars <- latent_vars[vapply(
      covariates,
      function(x) length(unique(stats::na.omit(x))) >= 2L,
      logical(1)
    )]
    status <- if (sum(seurat$cluster_group == ident_1) >= 10L &&
      sum(seurat$cluster_group == ident_2) >= 10L) {
      "tested"
    } else {
      "skipped_fewer_than_10_cells"
    }
    append_sample_report(
      report_file, celltype, ident_1, ident_2, seurat, valid_latent_vars, status
    )
  })
  invisible(report_file)
}

run_celltype_de <- function(celltype) {
  message("Running ", contrast_name, " in ", tissue_code_use, ": ", celltype)
  seurat <- prepare_celltype(celltype)

  message("Cells retained by condition and subtype:")
  print(table(seurat$sub_celltype_v2, seurat$condition_ema_emco))

  subtypes <- shared_subtypes(seurat)
  if (!length(subtypes)) {
    stop("No shared subcelltypes in ", celltype, " for ", tissue_code_use, " EMA vs EMCO.")
  }
  report_file <- sample_report_file(celltype)
  initialize_sample_report(report_file)

  # Tissue is deliberately not a latent variable: it is fixed by this script.
  markers <- map_dfr(subtypes, function(subtype) {
    run_mast(
      seurat,
      ident_1 = paste0(subtype, "_EMA"),
      ident_2 = paste0(subtype, "_EMCO"),
      latent_vars = c("nCount_RNA", "current_horm_treatm", "lib_type"),
      report_file = report_file,
      main_celltype = celltype
    )
  })
  if (!nrow(markers)) {
    return(tibble(
      gene = character(), p_val = numeric(), avg_log2FC = numeric(),
      pct.1 = numeric(), pct.2 = numeric(), p_val_adj = numeric(),
      diffexp = character()
    ))
  }
  filter(markers, !gene %in% c("HBA2", "HBB"))
}

if (report_only) {
  walk(celltypes, write_celltype_sample_report)
  message("Completed sample reports for ", contrast_name, " in ", tissue_code_use, ". Outputs: ", de_outdir)
  quit(save = "no")
}

set.seed(42)
de_results <- setNames(lapply(celltypes, run_celltype_de), celltypes)

saveRDS(
  de_results,
  file.path(de_outdir, paste0("DE_", contrast_name, "_", tissue_code_use, "_celltypes.Rds"))
)

workbook <- createWorkbook()
for (celltype in names(de_results)) {
  addWorksheet(workbook, celltype)
  writeData(workbook, celltype, de_results[[celltype]])
}
saveWorkbook(
  workbook,
  file.path(de_outdir, paste0(contrast_name, "_", tissue_code_use, "_tables.xlsx")),
  overwrite = TRUE
)

if (!report_only) plan(sequential)
message("Completed ", contrast_name, " in ", tissue_code_use, ". Outputs: ", de_outdir)
