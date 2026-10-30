# Add EMA/EMCO/ADENO condition metavariables to annotated zoom-in objects.
#
# Tissue-specific variables are defined by anatomical tissue code for DE-eligible
# cells only: EN contains endometrium (EN) and endometrial lesions (LE); MI
# contains myometrium (MI), myometrial lesions (LM), and LE samples from the
# multiome cohort. snGEX cells and patient 23 cells are excluded from DE and
# therefore retain NA for these variables.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(openxlsx)
})

input_files <- c(
  smc = "paper_output/zoom-in-mito/data/annotations/smc/sc_smc_annoted.Rds",
  endothelial = "paper_output/zoom-in-mito/data/annotations/endothelial/sc_endothelial_annoted.Rds",
  immune = "paper_output/zoom-in-mito/data/annotations/immune/sc_immune_annoted.Rds",
  stromal = "paper_output/zoom-in-mito/data/annotations/stromal/sc_stromal_annoted.Rds",
  pv = "paper_output/zoom-in-mito/data/annotations/pv/sc_pv_annoted.Rds",
  epithelial = "paper_output/zoom-in-mito/data/annotations/epithelial/sc_epithelial_annoted.Rds"
)
output_dir <- "paper_output/zoom-in-mito/data/annotations_conditions"

updated_metadata <- openxlsx::read.xlsx("metadata/metadata_2026.xlsx")
if (anyDuplicated(updated_metadata$bioinformatics_id)) {
  stop("metadata_2026.xlsx has duplicated bioinformatics_id values.")
}

condition_levels <- c("EMCO", "EMA", "ADENO")
condition_metavariables_to_report <- c(
  "condition_ema_emco_adeno_en",
  "condition_ema_emco_adeno_mi"
)
multiome_sample_ids <- c(
  "MU08ADEN", "MU08ADMI", "MU08ADLE", "MU11COEN", "MU11COMI",
  "MU12COEN", "MU12COMI", "MU13ADEN", "MU13ADMI", "MU13ADLE",
  "MU14ADEN", "MU14ADMI", "MU14ADLE", "MU18COEN", "MU18COMI"
)

add_condition_metavariables <- function(input_file, celltype) {
  seu <- LoadSeuratRds(input_file)

  metadata_index <- match(
    seu$bioinformatics_id,
    updated_metadata$bioinformatics_id
  )
  if (anyNA(metadata_index)) {
    stop(
      "Missing updated metadata for ", celltype, ": ",
      paste(unique(seu$bioinformatics_id[is.na(metadata_index)]), collapse = ", ")
    )
  }
  metadata_columns <- setdiff(colnames(updated_metadata), "bioinformatics_id")
  metadata_add <- updated_metadata[metadata_index, metadata_columns, drop = FALSE]
  rownames(metadata_add) <- Cells(seu)
  seu@meta.data[, metadata_columns] <- metadata_add

  eligible_for_de <-
    !is.na(seu$lib_type) & seu$lib_type != "snGEX" &
      !is.na(seu$patient) & as.character(seu$patient) != "23"

  condition <- case_when(
    eligible_for_de & seu$group == "AD" & seu$tissue_code %in% c("EN", "MI") ~ "EMA",
    eligible_for_de & seu$group == "CO" & seu$tissue_code %in% c("EN", "MI") ~ "EMCO",
    eligible_for_de & seu$group == "AD" & seu$tissue_code %in% c("LE", "LM") ~ "ADENO",
    TRUE ~ NA_character_
  )
  seu$condition_ema_emco_adeno <- factor(condition, levels = condition_levels)

  condition_en <- if_else(
    seu$tissue_code %in% c("EN", "LE"),
    condition,
    NA_character_
  )
  seu$condition_ema_emco_adeno_en <- factor(condition_en, levels = condition_levels)

  multiome_le_for_mi <-
    seu$bioinformatics_id %in% multiome_sample_ids & seu$tissue_code == "LE"
  condition_mi <- if_else(
    seu$tissue_code %in% c("MI", "LM") | multiome_le_for_mi,
    condition,
    NA_character_
  )
  seu$condition_ema_emco_adeno_mi <- factor(condition_mi, levels = condition_levels)

  seu
}

print_tissue_condition_membership <- function(seu, celltype) {
  for (metavariable in condition_metavariables_to_report) {
    values <- as.character(seu@meta.data[[metavariable]])
    samples_by_condition <- lapply(condition_levels, function(condition_level) {
      sort(unique(as.character(seu$bioinformatics_id[
        !is.na(values) & values == condition_level
      ])))
    })
    names(samples_by_condition) <- condition_levels
    max_samples <- max(lengths(samples_by_condition), 0L)
    condition_samples <- as.data.frame(lapply(samples_by_condition, function(samples) {
      length(samples) <- max_samples
      samples
    }), check.names = FALSE)

    message("Sample membership for ", celltype, " / ", metavariable, ":")
    print(condition_samples, row.names = FALSE, right = FALSE)
  }
}

# Complete this reporting pass before any condition-annotated object is saved.
# The objects are deliberately reloaded during the save pass to avoid retaining
# all six large Seurat objects in memory at once.
for (celltype in names(input_files)) {
  seu <- add_condition_metavariables(input_files[[celltype]], celltype)
  print_tissue_condition_membership(seu, celltype)
  rm(seu)
  gc(verbose = FALSE)
}

message("All EN and MI sample-membership reports have been printed. Beginning save pass.")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (celltype in names(input_files)) {
  seu <- add_condition_metavariables(input_files[[celltype]], celltype)
  output_file <- file.path(output_dir, paste0("sc_", celltype, "_annoted_conditions.Rds"))
  SaveSeuratRds(seu, output_file)
  message("Saved ", celltype, ": ", output_file)
  rm(seu)
  gc(verbose = FALSE)
}
