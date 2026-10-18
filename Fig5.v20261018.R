# Figure 5: tumor microenvironment, C1Q+ TAMs, and spatial validation
#
# Expected repository layout:
#   data/fig5/
#   figures/       optional output directory
#   source_data/   optional source-data directory
#
# The script defines fig5a, fig5d, fig5e_left, fig5e_right, fig5g,
# fig5g_right, and fig5h. It does not force a working directory, install
# packages, or save figures automatically.

required_packages <- c(
  "here", "Seurat", "dplyr", "tidyr", "ggplot2", "patchwork",
  "readr", "STDistance"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install the missing packages before running Figure 5: ",
       paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(STDistance)
})

set.seed(20261018)
fig5_data_dir <- here::here("data", "fig5")
fig5_file <- function(name) file.path(fig5_data_dir, name)

tissue_from_id <- function(ids) {
  code <- gsub("[0-9]", "", ids)
  ifelse(code == "SN", "GC", ifelse(code == "ST", "GC",
    ifelse(code == "PT", "AGC", ifelse(code == "PTH", "ALi",
      ifelse(code == "PN", "AGC", "APB")))))
}

# -------------------------------------------------------------------------
# Figure 5a: paired TME composition in AGC and GC
# -------------------------------------------------------------------------
tme_object <- readRDS(fig5_file("tme_cells.rds"))
tme_meta <- as.data.frame(tme_object@meta.data, stringsAsFactors = FALSE)
tme_meta$paper_annotation <- as.character(tme_meta$paper_annotation)
tme_meta$tissue_type <- tissue_from_id(tme_meta$orig.ident)
tme_meta$tissue_type[tme_meta$orig.ident == "HCC1T"] <- "HCC"

tme_counts <- as.data.frame.matrix(table(tme_meta$paper_annotation, tme_meta$orig.ident))
tme_percent <- sweep(tme_counts, 2, colSums(tme_counts), "/") * 100
fig5a_sample <- as.data.frame(t(tme_percent), stringsAsFactors = FALSE)
fig5a_sample$orig_ident <- rownames(fig5a_sample)
fig5a_sample$tissue_type <- tme_meta$tissue_type[match(fig5a_sample$orig_ident, tme_meta$orig.ident)]
fig5a_data <- fig5a_sample %>%
  pivot_longer(c(-orig_ident, -tissue_type), names_to = "cell_type", values_to = "percent") %>%
  mutate(
    cell_type = ifelse(cell_type == "Monocyte cell", "Myeloid cell", cell_type),
    patient = gsub("N|T|S", "", orig_ident),
    condition = ifelse(grepl("N", orig_ident), "normal", "tumor"),
    group = ifelse(tissue_type == "AGC", ifelse(condition == "normal", "AGCn", "AGC"),
                   ifelse(condition == "normal", "GCn", "GC")),
    comparison = tissue_type
  ) %>%
  filter(tissue_type %in% c("AGC", "GC"),
         cell_type %in% c("Myeloid cell", "Endothelial cell", "Fibroblast cell"),
         !(tissue_type == "AGC" & orig_ident == "AGC093T2"))

make_tme_panel <- function(data, comparison, cell_type, ymax = 45) {
  panel_data <- data %>% filter(comparison == !!comparison, cell_type == !!cell_type)
  colors <- if (comparison == "AGC") c(AGCn = "#FFD27F", AGC = "#FFA500") else
    c(GCn = "#D7E7FF", GC = "#0079F4")
  ggplot(panel_data, aes(group, percent, fill = group)) +
    geom_boxplot(width = 0.28, outlier.shape = NA, show.legend = FALSE) +
    geom_line(aes(group = patient), color = "#00BFFF", linetype = "dashed", linewidth = 0.25, show.legend = FALSE) +
    geom_point(size = 0.9, show.legend = FALSE) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    labs(title = cell_type, x = NULL, y = "% among all TME cells") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(size = 9), axis.text.x = element_text(size = 8),
          axis.title.y = element_text(size = 8))
}
fig5a_cell_types <- c("Myeloid cell", "Endothelial cell", "Fibroblast cell")
fig5a <- (
  make_tme_panel(fig5a_data, "AGC", fig5a_cell_types[1]) |
    make_tme_panel(fig5a_data, "AGC", fig5a_cell_types[2]) |
    make_tme_panel(fig5a_data, "AGC", fig5a_cell_types[3])
) / (
  make_tme_panel(fig5a_data, "GC", fig5a_cell_types[1]) |
    make_tme_panel(fig5a_data, "GC", fig5a_cell_types[2]) |
    make_tme_panel(fig5a_data, "GC", fig5a_cell_types[3])
)

# -------------------------------------------------------------------------
# Shared helpers for C1Q+ TAM score and spatial distance panels
# -------------------------------------------------------------------------
extract_c1q_score <- function(ucell_file, spatial_file, tissue, sample_id) {
  ucell <- readRDS(ucell_file)
  spatial <- readRDS(spatial_file)
  spatial@meta.data$predicted_celltype[spatial@meta.data$predicted_celltype == "Monocyte cell"] <- "Myeloid cell"
  ids <- rownames(spatial@meta.data)[spatial@meta.data$predicted_celltype == "Myeloid cell"]
  ids <- intersect(ids, rownames(ucell@meta.data))
  data.frame(
    cell_id = ids, orig_ident = sample_id, tissue_type = tissue,
    C1Q_TAM_score = as.numeric(ucell@meta.data[ids, "C1Q.TAM_UCell"]),
    stringsAsFactors = FALSE
  )
}

score_boxplot <- function(data, y_label, colors, y_limits, y_breaks) {
  ggplot(data, aes(tissue_type, .data[[y_label]])) +
    geom_jitter(aes(color = tissue_type), width = 0.2, size = 1.5, alpha = 0.6) +
    geom_boxplot(aes(fill = tissue_type), width = 0.45, outlier.shape = NA, alpha = 0.7) +
    scale_color_manual(values = colors) + scale_fill_manual(values = colors) +
    scale_y_continuous(limits = y_limits, breaks = y_breaks, expand = c(0, 0)) +
    ylab(y_label) + xlab(NULL) + theme_classic() + theme(legend.position = "none")
}

# -------------------------------------------------------------------------
# Figure 5d: C1Q+ TAM score in AGC102 primary and liver-metastatic samples
# -------------------------------------------------------------------------
fig5d_data <- bind_rows(
  extract_c1q_score(fig5_file("AGC102T_ucell.rds"), fig5_file("AGC102T_spatial.rds"), "AGC", "AGC102T"),
  extract_c1q_score(fig5_file("AGC102M_ucell.rds"), fig5_file("AGC102M_spatial.rds"), "ALi", "AGC102M")
) %>% filter(!is.na(C1Q_TAM_score))
fig5d <- score_boxplot(fig5d_data, "C1Q_TAM_score",
                       c(AGC = "#FFD700", ALi = "#eb1d37"), c(-0.01, 0.4), seq(0, 0.4, 0.1)) +
  ylab("C1Q+ TAM score")

get_spatial_coordinates <- function(object) {
  coordinates <- GetTissueCoordinates(object)
  if (!all(c("x", "y") %in% names(coordinates))) names(coordinates)[1:2] <- c("x", "y")
  coordinates$barcode <- rownames(coordinates)
  coordinates
}

nearest_distance_by_status <- function(ucell_file, spatial_file, tissue, sample_id, allowed_types, y_max) {
  ucell <- readRDS(ucell_file)
  spatial <- readRDS(spatial_file)
  ids <- intersect(rownames(spatial@meta.data), rownames(ucell@meta.data))
  metadata <- spatial@meta.data[ids, , drop = FALSE]
  scores <- as.numeric(ucell@meta.data[ids, "C1Q.TAM_UCell"])
  names(scores) <- ids
  afp <- as.numeric(GetAssayData(spatial, assay = "Spatial", layer = "counts")["AFP", ids]) > 0
  status <- rep("others", length(ids)); names(status) <- ids
  status[metadata$predicted_celltype == "Monocyte cell" & scores > 0.001] <- "high C1Q+ TAM"
  status[metadata$predicted_celltype == "Monocyte cell" & scores <= 0.001] <- "low C1Q+ TAM"
  if (length(allowed_types) == 0) {
    status[afp] <- "AFP+ cell"
  } else {
    status[afp & metadata$predicted_celltype %in% allowed_types] <- "AFP+ cell"
  }
  keep <- status %in% c("high C1Q+ TAM", "low C1Q+ TAM", "AFP+ cell")
  coordinates <- get_spatial_coordinates(spatial)[ids[keep], , drop = FALSE]
  coordinates$C1QTAM_status <- status[rownames(coordinates)]
  distances <- bind_rows(lapply(c("high C1Q+ TAM", "low C1Q+ TAM"), function(reference) {
    subset_data <- coordinates[coordinates$C1QTAM_status %in% c(reference, "AFP+ cell"), , drop = FALSE]
    result <- calculate_nearest_distances(
      spatial_data = subset_data, reference_type = reference,
      target_types = "AFP+ cell", x_col = "x", y_col = "y",
      id_col = "barcode", type_col = "C1QTAM_status"
    )
    result$type <- reference
    result
  }))
  distances$distance <- as.numeric(distances[["AFP+ cell"]])
  distances$orig_ident <- sample_id
  distances$tissue_type <- tissue
  list(data = distances, y_max = y_max)
}

# -------------------------------------------------------------------------
# Figure 5e: spatial distance to AFP+ cells for AGC102
# -------------------------------------------------------------------------
fig5e_left_result <- nearest_distance_by_status(
  fig5_file("AGC102T_ucell.rds"), fig5_file("AGC102T_spatial.rds"), "AGC", "AGC102T",
  c("Epithelial cell", "Endothelial cell", "Fibroblast cell"), 125
)
fig5e_right_result <- nearest_distance_by_status(
  fig5_file("AGC102M_ucell.rds"), fig5_file("AGC102M_spatial.rds"), "ALi", "AGC102M",
  c("Epithelial cell", "Hepatocytes", "Fibroblast cell"), 20
)
make_distance_plot <- function(result) {
  data <- result$data
  data$type <- factor(data$type, levels = c("high C1Q+ TAM", "low C1Q+ TAM"))
  box_cols <- c("high C1Q+ TAM" = "#3e5584", "low C1Q+ TAM" = "#469c88")
  ggplot(data, aes(type, distance)) +
    geom_jitter(aes(color = type), width = 0.2, size = 1.2, alpha = 0.6) +
    geom_boxplot(aes(fill = type), outlier.shape = NA, alpha = 0.7) +
    scale_color_manual(values = box_cols) + scale_fill_manual(values = box_cols) +
    scale_y_continuous(limits = c(-0.01, result$y_max), expand = c(0, 0)) +
    ylab("Nearest distance from AFP+ cell to TAM") + xlab(NULL) + theme_classic() +
    theme(legend.position = "none")
}
fig5e_left <- make_distance_plot(fig5e_left_result)
fig5e_right <- make_distance_plot(fig5e_right_result)

# -------------------------------------------------------------------------
# Figure 5g-I: C1Q+ TAM score in P10 primary and liver-metastatic samples
# -------------------------------------------------------------------------
extract_myeloid_score <- function(ucell_file, spatial_file, tissue, sample_id) {
  ucell <- readRDS(ucell_file); spatial <- readRDS(spatial_file)
  metadata <- spatial@meta.data
  metadata$first_type[metadata$first_type == "Monocyte cell"] <- "Myeloid cell"
  metadata$second_type[metadata$second_type == "Monocyte cell"] <- "Myeloid cell"
  merged_type <- ifelse(!is.na(metadata$second_type), metadata$second_type,
                        ifelse(!is.na(metadata$first_type), metadata$first_type, "Unknown"))
  ids <- intersect(rownames(metadata)[merged_type == "Myeloid cell"], rownames(ucell@meta.data))
  data.frame(barcode = ids, C1Q_TAM_UCell = as.numeric(ucell@meta.data[ids, "C1Q.TAM_UCell"]),
             tissue_type = tissue, sample = sample_id, stringsAsFactors = FALSE)
}
fig5g_data <- bind_rows(
  extract_myeloid_score(fig5_file("P10T_ucell.rds"), fig5_file("P10T_spatial.rds"), "AGC", "P10T"),
  extract_myeloid_score(fig5_file("P10TH_ucell.rds"), fig5_file("P10TH_spatial.rds"), "ALi", "P10TH")
) %>% mutate(tissue_type = factor(tissue_type, levels = c("AGC", "ALi")))
fig5g <- score_boxplot(fig5g_data, "C1Q_TAM_UCell",
                       c(AGC = "#FFD700", ALi = "#eb1d37"), c(-0.01, 0.6), seq(0, 0.6, 0.2)) +
  ylab("C1Q+ TAM score")

fig5g_right_result <- nearest_distance_by_status(
  fig5_file("P10T_ucell.rds"), fig5_file("P10T_spatial.rds"), "AGC", "P10T",
  character(), 30
)
fig5g_right <- make_distance_plot(fig5g_right_result)

# -------------------------------------------------------------------------
# Figure 5h: organoid validation of C1Q+ TAM score
# -------------------------------------------------------------------------
organoid_data <- read.delim(fig5_file("organoid_ucell.tsv"), check.names = FALSE,
                            stringsAsFactors = FALSE)
fig5h_data <- organoid_data[, c("cell.id", "orig.ident", "COMPLEMENT.hi_UCell")]
names(fig5h_data) <- c("cell_id", "orig_ident", "C1Q_TAM_score")
fig5h_data$sample_source <- ifelse(fig5h_data$orig_ident == "Org.Neg", "GC organoid", "AGC organoid")
fig5h_data$organoid_status <- ifelse(fig5h_data$orig_ident == "Org.Neg", "AFP-", "AFP+")
fig5h_data$group <- factor(ifelse(fig5h_data$orig_ident == "Org.Neg", "GC.org", "AGC.org"),
                           levels = c("GC.org", "AGC.org"))
fig5h_data$C1Q_TAM_score <- as.numeric(fig5h_data$C1Q_TAM_score)
fig5h <- ggplot(fig5h_data, aes(group, C1Q_TAM_score)) +
  geom_jitter(aes(color = group), width = 0.22, size = 1.2, alpha = 0.60, stroke = 0) +
  geom_boxplot(aes(fill = group), width = 0.42, outlier.shape = NA, alpha = 0.72, color = "black", linewidth = 0.55) +
  scale_color_manual(values = c("GC.org" = "#57c3f4", "AGC.org" = "#FFA500")) +
  scale_fill_manual(values = c("GC.org" = "#57c3f4", "AGC.org" = "#FFA500")) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2), expand = c(0, 0)) +
  labs(x = NULL, y = "C1Q+ TAM score") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none", axis.text = element_text(color = "black"),
        axis.title.y = element_text(color = "black"), axis.line = element_line(color = "black", linewidth = 0.55),
        axis.ticks = element_line(color = "black", linewidth = 0.55), plot.margin = margin(5, 5, 5, 5))

fig5a
fig5d
fig5e_left
fig5e_right
fig5g
fig5g_right
fig5h
