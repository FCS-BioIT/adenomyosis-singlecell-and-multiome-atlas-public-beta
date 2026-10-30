# Differential Expression Pipeline ---------------------------------------------
# Final-DE variant: excludes snGEX and patient 23 from every contrast.

library(Seurat)
library(harmony)
library(tibble)
library(stringi)
library(stringr)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(future)

plan(multisession, workers = 16)
options(future.globals.maxSize = 8000 * 1024^2)

# Define consistent directories
input_dir <- "paper_output/zoom-in-mito/"
output_dir <- "paper_output/de-zoom-minus-rb-mt-nopat-no23-nosn/"

# METADATA LOAD:
updated_md <- openxlsx::read.xlsx(xlsxFile = "metadata/metadata_2026.xlsx")

stopifnot(!anyDuplicated(updated_md$bioinformatics_id))

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

append_sample_report <- function(report_file, main_celltype, ident_list,
                                 seurat, valid_latent_vars, status) {
  if (is.null(report_file)) {
    return(invisible(NULL))
  }

  condition_labels <- sub("^.*_", "", ident_list)
  subtype <- sub(paste0("_", condition_labels[1], "$"), "", ident_list[1])
  sample_rows <- lapply(ident_list, function(ident) {
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
  invisible(NULL)
}

# DEFINE DE FUNCTION -----------------------------------------------------------
diffexp_conditions <- function(
    ident_list,
    latent_vars = c("nCount_RNA", "current_horm_treatm", "tissue", "lib_type"),
    seurat,
    report_file = NULL,
    main_celltype = NA_character_) {
  # print(ident_list)

  if (length(unique(as.character(seurat@meta.data[, "protocol"]))) == 1) {
    latent_vars <- latent_vars
    print("USING REDUCED MODEL FOR 1 PROTOCOL ONLY")
  } else {
    latent_vars <- latent_vars
  }
  print("CHECK USED COVARIATES: ")

  cells_use <- rownames(seurat@meta.data)[
    seurat@meta.data$cluster_group %in% ident_list
  ]

  print(paste0("Levels in the cells covariates in celltype ", ident_list[1], ":"))
  seurat@meta.data[cells_use, latent_vars, drop = FALSE] %>%
    summarise(across(everything(), ~ length(unique(na.omit(.x))))) %>%
    print()

  md_use <- seurat@meta.data[cells_use, latent_vars, drop = FALSE]

  valid_latent_vars <- latent_vars[
    sapply(md_use, function(v) {
      length(unique(na.omit(v))) >= 2
    })
  ]

  print("Valid latent variables used:")
  print(valid_latent_vars)


  n_cells_1 <- sum(seurat@meta.data[, "cluster_group"] == ident_list[1])
  n_cells_2 <- sum(seurat@meta.data[, "cluster_group"] == ident_list[2])
  report_status <- if (n_cells_1 >= 10L && n_cells_2 >= 10L) {
    "tested"
  } else {
    "skipped_fewer_than_10_cells"
  }
  append_sample_report(
    report_file = report_file,
    main_celltype = main_celltype,
    ident_list = ident_list,
    seurat = seurat,
    valid_latent_vars = valid_latent_vars,
    status = report_status
  )

  if (n_cells_1 >= 10L && n_cells_2 >= 10L) {
    markers_df <- FindMarkers(
      seurat,
      ident.1 = ident_list[1],
      ident.2 = ident_list[2],
      assay = "RNA",
      min.pct = 0.25,
      test.use = "MAST",
      latent.vars = valid_latent_vars,
      random.seed = 1
    )
    markers_df$diffexp <- rep(
      paste(
        ident_list[1],
        "-vs-",
        ident_list[2],
        collapse = "",
        sep = ""
      ),
      nrow(markers_df)
    )
    markers_df <- markers_df %>%
      rownames_to_column(var = "gene")
    return(markers_df)
  } else {
    print("One of the two cell groups has fewer than 10 cells")
    print(paste0("N cells group ", ident_list[1]))
    print(n_cells_1)
    print(paste0("N cells group ", ident_list[2]))
    print(n_cells_2)
    return(NULL)
  }
}

set.seed(42)
# VARIABLES:
celltypes <- list(
  pv = c("pv", "pv"),
  immune = c("immune", "immune"),
  endothelial = c("endothelial", "endo"),
  smc = c("smc", "smc"),
  epithelial = c("epithelial", "epi"),
  stromal = c("stromal", "stromal")
)

comparisons <- list(
  ADLE_ADUT = c("ADLE", "ADUT"),
  ADLE_CO = c("ADLE", "CO"),
  ADUT_CO = c("ADUT", "CO")
)

# Technologies
tech_sc <- "scGEX"
tech_mu <- "snMU"
tech_sn <- "snGEX"
tech_nu <- c("snGEX", "snMU")

# ------------------------------------------------------------------
# Que hay de distinto entre los úteros para desarrollar la enfermedad (therefore, no lesiones)
# EMA vs EMCO
de_ema_emco <- function(
    celltype,
    cell_nick,
    tech = NULL) {
  print(celltype)
  print(tech)

  sc_obj <- LoadSeuratRds(
    file.path(
      input_dir,
      "data/annotations",
      celltype,
      paste0(
        "sc_",
        celltype,
        "_annoted.Rds"
      )
    )
  )

  # ADD METADATA
  idx <- match(sc_obj$bioinformatics_id, updated_md$bioinformatics_id)
  if (anyNA(idx)) {
    stop(
      "Missing metadata for: ",
      paste(unique(sc_obj$bioinformatics_id[is.na(idx)]), collapse = ", ")
    )
  }
  md_cols <- setdiff(colnames(updated_md), "bioinformatics_id")
  new_md <- updated_md[idx, md_cols, drop = FALSE]
  rownames(new_md) <- rownames(sc_obj@meta.data)

  sc_obj@meta.data[, md_cols] <- new_md

  # Final-DE exclusions apply before contrast-specific filtering.
  sc_obj <- subset(sc_obj, subset = lib_type != "snGEX" & patient != 23)

  # Important: patient should be categorical, age numeric
  sc_obj$patient <- factor(sc_obj$patient)
  sc_obj$group <- factor(sc_obj$group)
  sc_obj$protocol <- factor(sc_obj$protocol)
  sc_obj$lib_type <- factor(sc_obj$lib_type)
  sc_obj$age <- as.numeric(sc_obj$age)

  if (!is.null(tech)) {
    sc_obj <- subset(sc_obj, subset = lib_type %in% tech)
  }


  sc_obj$tissue_abbr <- sc_obj$tissue_code
  sc_obj$tissue_abbr[sc_obj$tissue_abbr == "LM"] <- "LE"
  sc_obj$group_tissue <- paste(sc_obj$group, sc_obj$tissue_abbr, sep = "")

  sc_obj$cluster_group <- paste(sc_obj$sub_celltype_v2, sc_obj$group, sep = "_")

  print(table(sc_obj$cluster_group))

  # remove mt rb genes
  genes <- rownames(sc_obj[["RNA"]])
  genes_mito <- grep("^MT-", genes, value = TRUE)
  genes_ribo <- grep("^RPS|^RPL", genes, value = TRUE)
  genes_remove <- unique(c(genes_mito, genes_ribo))
  genes_keep <- setdiff(genes, genes_remove)
  sc_obj <- subset(sc_obj, features = genes_keep)

  # Do not use lesions
  sc_obj <- subset(sc_obj, tissue != "lesion")

  sc_obj <- subset(sc_obj, lib_type != "snGEX")

  Idents(sc_obj) <- "cluster_group"

  cond1 <- "AD"
  cond2 <- "CO"

  pre_ident.1 <- unique(sc_obj$cluster_group[which(sc_obj$group == cond1)])
  post_ident.2 <- unique(sc_obj$cluster_group[which(sc_obj$group == cond2)])

  pre_ident.1 <- pre_ident.1[order(pre_ident.1)]
  post_ident.2 <- post_ident.2[order(post_ident.2)]

  idents_merger <- t(rbind(pre_ident.1, post_ident.2))
  idents_merger <- split(idents_merger, seq(nrow(idents_merger)))

  # defining latent vars:
  latent_vars_global <- c(
    "nCount_RNA",
    "current_horm_treatm",
    "tissue",
    "lib_type"
  )

  sc_conditions_cell_type_markers <- purrr::map_dfr(
    idents_merger,
    function(x) diffexp_conditions(ident_list = x, seurat = sc_obj, latent_vars = latent_vars_global)
  )

  sc_conditions_cell_type_markers <- sc_conditions_cell_type_markers[
    which(!sc_conditions_cell_type_markers$gene %in% c("HBA2", "HBB")),
  ]

  return(sc_conditions_cell_type_markers)
}


dir.create(paste0(output_dir, "data/DE"), recursive = TRUE)
# difexp_EMA_EMCO <- lapply(celltypes, function(x) de_ema_emco(x[1], x[2])) # DO NOT SEPARATE BY TECH


#
# saveRDS(
#   object = difexp_EMA_EMCO,
#   file = file.path(output_dir, "data/DE/DE_AD_CO_celltypes.Rds")
# )


# DE Multicondition ------------------------------------------------------------
de_multicond <- function(
    celltype,
    cell_nick,
    comparison,
    tissue = NULL,
    tech = NULL) {
  # Working on what celltype
  print(paste0("Celltype: ", celltype))

  # Load seurat object
  sc_obj <- LoadSeuratRds(
    file.path(
      input_dir,
      "data/annotations", celltype, paste0("sc_", celltype, "_annoted.Rds")
    )
  )

  # filter technology
  if (!is.null(tech)) {
    sc_obj <- subset(sc_obj, subset = lib_type %in% tech)
  }

  # remove mt rb genes
  genes <- rownames(sc_obj[["RNA"]])
  genes_mito <- grep("^MT-", genes, value = TRUE)
  genes_ribo <- grep("^RPS|^RPL", genes, value = TRUE)
  genes_remove <- unique(c(genes_mito, genes_ribo))
  genes_keep <- setdiff(genes, genes_remove)
  sc_obj <- subset(sc_obj, features = genes_keep)

  # ADD METADATA
  idx <- match(sc_obj$bioinformatics_id, updated_md$bioinformatics_id)
  if (anyNA(idx)) {
    stop(
      "Missing metadata for: ",
      paste(unique(sc_obj$bioinformatics_id[is.na(idx)]), collapse = ", ")
    )
  }
  md_cols <- setdiff(colnames(updated_md), "bioinformatics_id")
  new_md <- updated_md[idx, md_cols, drop = FALSE]
  rownames(new_md) <- rownames(sc_obj@meta.data)

  sc_obj@meta.data[, md_cols] <- new_md

  # Final-DE exclusions apply before contrast-specific filtering.
  sc_obj <- subset(sc_obj, subset = lib_type != "snGEX" & patient != 23)


  # contrast metadata
  sc_obj$tissue_abbr <- sc_obj$tissue_code
  sc_obj$tissue_abbr[sc_obj$tissue_abbr == "LM"] <- "LE"
  sc_obj$group_tissue <- paste(sc_obj$group, sc_obj$tissue_abbr, sep = "")

  sc_obj$group_tissue_aggr <- sc_obj$group_tissue
  table(sc_obj$group_tissue_aggr)
  sc_obj$group_tissue_aggr[
    sc_obj$group_tissue_aggr %in% c("COEN", "COMI")
  ] <- "EMCO"

  sc_obj$group_tissue_aggr[
    sc_obj$group_tissue_aggr %in% c("ADEN", "ADMI")
  ] <- "EMA"

  sc_obj$group_tissue_aggr[sc_obj$group_tissue_aggr == "ADLE"] <- "ADENO"


  # Filter by comparison if not balanced lib type
  if (setequal(as.character(comparison), c("ADENO", "EMCO"))) {
    sc_obj <- subset(sc_obj, subset = lib_type != "snGEX")
  }
  # remove patient that was reclassified
  if ("ADENO" %in% as.character(comparison)) {
    sc_obj <- subset(sc_obj, subset = patient != 23)
  }

  # Filter by tissue
  comparator <- as.character(comparison[2])
  if (is.null(tissue)) {
    sc_obj <- subset(
      sc_obj,
      subset =
        (
          lib_type %in% c("scGEX", "snMU")
        )
    )

    print(paste0("All tissues selected"))
    print("The lib types selected are: ")
    print(table(sc_obj$lib_type))
    print("The tissue code selected are: ")
    print(table(sc_obj$tissue_code))
    print("The groups and tissue are: ")
    print(table(sc_obj$group_tissue_aggr))
    print(table(sc_obj$patient))
    print("Tissue code by lib type considered: ")
    print(table(sc_obj$lib_type, sc_obj$tissue_code))
  } else if (tissue == "MIO") {
    sc_obj <- subset(
      sc_obj,
      subset =
        (
          lib_type == "scGEX" &
            tissue_code == "LM" &
            group_tissue_aggr == "ADENO"
        ) |
          (
            lib_type %in% c("snGEX", "snMU") &
              tissue_code == "LE" &
              group_tissue_aggr == "ADENO"
          ) |
          (
            tissue_code == "MI" &
              group_tissue_aggr == comparator
          )
    )
    print(paste0("The tissue is: ", tissue))
    print("The lib types selected are: ")
    print(table(sc_obj$lib_type))
    print("The tissue code selected are: ")
    print(table(sc_obj$tissue_code))
    print("The groups and tissue are: ")
    print(table(sc_obj$group_tissue_aggr))
    print(table(sc_obj$patient))
    print("Tissue code by lib type considered: ")
    print(table(sc_obj$lib_type, sc_obj$tissue_code))
  } else if (tissue == "ENDO") {
    sc_obj <- subset(
      sc_obj,
      subset =
        (
          tissue_code == "LE" &
            group_tissue_aggr == "ADENO"
        ) |
          (
            tissue_code == "EN" &
              group_tissue_aggr == comparator
          )
    )
    print(paste0("The tissue is: ", tissue))
    print("The lib types selected are: ")
    print(table(sc_obj$lib_type))
    print("The tissue code selected are: ")
    print(table(sc_obj$tissue_code))
    print("The groups and tissue are: ")
    print(table(sc_obj$group_tissue_aggr))
    print(table(sc_obj$patient))
    print("Tissue code by lib type considered: ")
    print(table(sc_obj$lib_type, sc_obj$tissue_code))
  }

  print(table(sc_obj$tissue_code, sc_obj$group_tissue_aggr))

  sc_obj <- subset(
    sc_obj,
    subset = group_tissue_aggr == comparison[1] | group_tissue_aggr == comparison[2]
  )

  print(table(sc_obj$group_tissue_aggr))

  sc_obj$cluster_group <- paste(
    sc_obj$sub_celltype_v2, sc_obj$group_tissue_aggr,
    sep = "_"
  )

  tissue_label <- if (is.null(tissue)) "all_tissues" else tissue
  report_file <- file.path(
    output_dir,
    "data",
    "DE",
    "sample_reports",
    paste0(
      "samples_covariates_",
      tissue_label,
      "_",
      paste(as.character(comparison), collapse = "_"),
      "_",
      celltype,
      ".txt"
    )
  )
  initialize_sample_report(report_file)

  Idents(sc_obj) <- "cluster_group"

  cond1 <- comparison[1]
  cond2 <- comparison[2]

  pre_ident.1 <- unique(sc_obj$cluster_group[which(sc_obj$group_tissue_aggr == cond1)])
  post_ident.2 <- unique(sc_obj$cluster_group[which(sc_obj$group_tissue_aggr == cond2)])

  pre_name <- sub(paste0("_", cond1, "$"), "", pre_ident.1)
  post_name <- sub(paste0("_", cond2, "$"), "", post_ident.2)

  names(pre_ident.1) <- pre_name
  names(post_ident.2) <- post_name

  shared <- sort(intersect(names(pre_ident.1), names(post_ident.2)))

  missing_in_cond1 <- setdiff(names(post_ident.2), names(pre_ident.1))
  missing_in_cond2 <- setdiff(names(pre_ident.1), names(post_ident.2))

  if (length(missing_in_cond1) > 0) {
    message("Missing in ", cond1, ": ", paste(missing_in_cond1, collapse = ", "))
  }

  if (length(missing_in_cond2) > 0) {
    message("Missing in ", cond2, ": ", paste(missing_in_cond2, collapse = ", "))
  }

  if (length(shared) == 0) {
    stop(glue::glue("No shared subcelltypes for {cond1} vs {cond2}."))
  }

  idents_merger <- Map(c, pre_ident.1[shared], post_ident.2[shared])
  names(idents_merger) <- shared

  print("Ident list being passed to diffexp_conditions():")

  print(idents_merger)

  latent_vars_lesion <- if (is.null(tissue)) {
    c(
      "nCount_RNA",
      "tissue",
      "current_horm_treatm",
      "lib_type"
    )
  } else {
    c(
      "nCount_RNA",
      "current_horm_treatm",
      "lib_type"
    )
  }

  print("Latent variables used: ")

  sapply(latent_vars_lesion, FUN = function(x) {
    print(x)
  })

  sc_conditions_cell_type_markers <- purrr::map_dfr(
    idents_merger,
    function(x) {
      diffexp_conditions(
        ident_list = x,
        seurat = sc_obj,
        latent_vars = latent_vars_lesion,
        report_file = report_file,
        main_celltype = celltype
      )
    }
  )

  sc_conditions_cell_type_markers <- sc_conditions_cell_type_markers[which(!sc_conditions_cell_type_markers$gene %in% c("HBA2", "HBB")), ]

  return(sc_conditions_cell_type_markers)
}


comparisons <- list(
  adeno_ema_mi = list(
    conds = factor(c("ADENO", "EMA"), levels = c("EMA", "ADENO")),
    tissue = "MIO"
  ),
  adeno_ema_en = list(
    conds = factor(c("ADENO", "EMA"), levels = c("EMA", "ADENO")),
    tissue = "ENDO"
  ),
  adeno_emco_mi = list(
    conds = factor(c("ADENO", "EMCO"), levels = c("EMCO", "ADENO")),
    tissue = "MIO"
  ),
  adeno_emco_en = list(
    conds = factor(c("ADENO", "EMCO"), levels = c("EMCO", "ADENO")),
    tissue = "ENDO"
  )
)

# Lesiones vs tejido sano por regiones uterinas
# ADENO VS EMA MYOMETRIUM ------------------------------------------------------
print("ADENO vs EMA MYOMETRIUM")
diffexp_ADENO_EMA_MI <- lapply(
  celltypes, function(x) {
    de_multicond(
      x[1],
      x[2],
      comparisons$adeno_ema_mi$conds,
      tissue = comparisons$adeno_ema_mi$tissue
    )
  }
)
print("ADENO vs EMCO MYOMETRIUM")
# ADENO VS EMCO MYOMETRIUM -----------------------------------------------------
diffexp_ADENO_EMCO_MI <- lapply(
  celltypes, function(x) {
    de_multicond(
      x[1],
      x[2],
      comparisons$adeno_emco_mi$conds,
      tissue = comparisons$adeno_emco_mi$tissue
    )
  }
)
print("ADENO vs EMA ENDOMETRIUM")
# ADENO VS EMA ENDOMETRIUM -----------------------------------------------------
diffexp_ADENO_EMA_EN <- lapply(
  celltypes, function(x) {
    de_multicond(
      x[1],
      x[2],
      comparisons$adeno_ema_en$conds,
      tissue = comparisons$adeno_ema_en$tissue
    )
  }
)
print("ADENO vs EMA ENDOMETRIUM")
# ADENO VS EMCO ENDOMETRIUM ----------------------------------------------------
diffexp_ADENO_EMCO_EN <- lapply(
  celltypes, function(x) {
    de_multicond(
      x[1],
      x[2],
      comparisons$adeno_emco_en$conds,
      tissue = comparisons$adeno_emco_en$tissue
    )
  }
)
print("EMA vs EMA ALL TISSUES")
# EMA VS EMCO ALL TISSUES
diffexp_EMA_EMCO <- lapply(
  celltypes, function(x) {
    de_multicond(
      x[1],
      x[2],
      c("EMA", "EMCO")
    )
  }
)


# Save results -----------------------------------------------------------------
save_diffexp_xlsx <- function(diffexp_contrast, contrast_name, tissue = NULL) {
  # Create a new workbook
  wb <- createWorkbook()

  for (sheet_name in names(diffexp_contrast)) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet = sheet_name, diffexp_contrast[[sheet_name]])
  }
  if (!is.null(tissue)) {
    if (length(tissue) == 1) {
      saveWorkbook(
        wb,
        file = file.path(
          output_dir,
          "data/DE/",
          paste0(
            contrast_name,
            "_",
            tissue,
            "_tables.xlsx"
          )
        ),
        overwrite = TRUE
      )
    } else {
      saveWorkbook(
        wb,
        file = file.path(
          output_dir,
          "data/DE/",
          paste0(
            contrast_name,
            "_",
            tissue[1],
            "_",
            tissue[2],
            "_tables.xlsx"
          )
        ),
        overwrite = TRUE
      )
    }
  } else {
    saveWorkbook(
      wb,
      file = file.path(
        output_dir,
        "data/DE/",
        paste0(
          contrast_name,
          "_tables.xlsx"
        )
      ),
      overwrite = TRUE
    )
  }
}

# AD CO
# diffexp_AD_CO <- readRDS("output/de-zoom-separate-tech/data/DE/DE_AD_CO_celltypes.Rds")
#
# save_diffexp_xlsx(diffexp_AD_CO$sc, contrast_name = "AD_CO", tech = NULL)

# LESION
save_diffexp_xlsx(diffexp_ADENO_EMA_MI, contrast_name = "ADENO_EMA", tissue = "MI")
save_diffexp_xlsx(diffexp_ADENO_EMCO_MI, contrast_name = "ADENO_EMCO", tissue = "MI")
save_diffexp_xlsx(diffexp_ADENO_EMA_EN, contrast_name = "ADENO_EMA", tissue = "EN")
save_diffexp_xlsx(diffexp_ADENO_EMCO_EN, contrast_name = "ADENO_EMCO", tissue = "EN")
save_diffexp_xlsx(diffexp_EMA_EMCO, contrast_name = "EMA_EMCO", tissue = NULL)



saveRDS(
  list(
    adeno_ema_mi = diffexp_ADENO_EMA_MI,
    adeno_emco_mi = diffexp_ADENO_EMCO_MI,
    adeno_ema_en = diffexp_ADENO_EMA_EN,
    adeno_emco_en = diffexp_ADENO_EMCO_EN,
    ema_emco = diffexp_EMA_EMCO
  ),
  file = file.path(output_dir, "data/DE/DE_ADENO_EMA_EMCO_celltypes.Rds")
)


# end-of-file
