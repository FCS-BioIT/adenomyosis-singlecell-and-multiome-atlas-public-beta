

library(Seurat)
library(ggplot2)
library(svglite)

data_dir <- "paper_output/zoom-in-mito/data/annotations_conditions"
output_dir <- "paper_output/zoom-in-mito/plots/featureplots"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seu <- LoadSeuratRds(
  file=file.path(
    data_dir,
    "sc_epithelial_annoted_conditions.Rds")
  )

feats <- c("MMP26", "ADAMTS8", "COL1A2", "TFPI2", "SERPINA5")

plt_en <- FeaturePlot(
  seu, 
  features=feats, 
  min.cutoff="q5", 
  max.cutoff="q95",
  split.by="condition_ema_emco_adeno_en",
  ncol=2
  )

ggsave(
  plot=plt_en,
  filename=file.path(output_dir, "epithelium_feats_en.png"),
  width=12, 
  height=length(feats)*5
  )

ggsave(
  plot=plt_en,
  filename=file.path(output_dir, "epithelium_feats_en.svg"),
  width=12, 
  height=length(feats)*5,
  fix_text_size=FALSE
  )


feats <- c("VIM", "SOX4")

plt_en <- FeaturePlot(
  seu, 
  features=feats, 
  min.cutoff="q5", 
  max.cutoff="q95",
  split.by="condition_ema_emco_adeno_en",
  ncol=2
  )

ggsave(
  plot=plt_en,
  filename=file.path(output_dir, "epithelium_feats_en_2.png"),
  width=12, 
  height=length(feats)*5
  )

ggsave(
  plot=plt_en,
  filename=file.path(output_dir, "epithelium_feats_en_2.svg"),
  width=12, 
  height=length(feats)*5,
  fix_text_size=FALSE
  )


feats <- c(
  "SOX9",
  "LEF1",
  "TCF7L2",
  "TCF7",
  "SOX13",
  "CTCF",
  "GRHL2",
  "TCF7L1"
  )

plt_en <- FeaturePlot(
  seu, 
  features=feats, 
  min.cutoff="q5", 
  max.cutoff="q95",
  split.by="condition_ema_emco_adeno_en",
  ncol=2
  )

ggsave(
  plot=plt_en,
  filename=file.path(output_dir, "epithelium_feats_en_3.png"),
  width=12, 
  height=length(feats)*5
  )

ggsave(
  plot=plt_en,
  filename=file.path(output_dir, "epithelium_feats_en_3.svg"),
  width=12, 
  height=length(feats)*5,
  fix_text_size=FALSE
  )




feats <- c(
  "DCN", "LUM", "ACTA2", "TAGLN", "PECAM1", "VWF"
  )

plt_en <- FeaturePlot(
  seu, 
  features=feats, 
  min.cutoff="q5", 
  max.cutoff="q95",
  split.by="condition_ema_emco_adeno_en",
  ncol=2
  )

ggsave(
  plot=plt_en,
  filename=file.path(output_dir, "epithelium_feats_en_contaminants.png"),
  width=12, 
  height=length(feats)*5
  )

ggsave(
  plot=plt_en,
  filename=file.path(output_dir, "epithelium_feats_en_contaminants.svg"),
  width=12, 
  height=length(feats)*5,
  fix_text_size=FALSE
  )


umap <- DimPlot(seu, split.by="lib_type")

ggsave(
  plot=umap,
  filename=file.path(output_dir, "epithelium_umap_tec.png"),
  width=15, 
  height=5
  )

ggsave(
  plot=umap,
  filename=file.path(output_dir, "epithelium_umap_tec.svg"),
  width=15, 
  height=5,
  fix_text_size=FALSE
  )

umap <- DimPlot(seu, split.by="protocol")

ggsave(
  plot=umap,
  filename=file.path(output_dir, "epithelium_umap_protocol.png"),
  width=15, 
  height=5
  )

ggsave(
  plot=umap,
  filename=file.path(output_dir, "epithelium_umap_protocol.svg"),
  width=15, 
  height=5,
  fix_text_size=FALSE
  )
