# Figure 3: epigenetic programs and pathway evolution along pseudotime
#
# Expected repository layout:
#   data/fig3/
#   figures/       optional output directory
#   source_data/   optional source-data directory
#
# The script defines fig3a, fig3g, and the source-data objects used by the
# plots. Figure 3b-c are represented by the three reported pseudotime cutoffs.
# No working directory is forced and no figures are saved automatically.

required_packages <- c(
  "here", "Seurat", "monocle", "dplyr", "ggplot2", "ggrepel",
  "patchwork", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install the missing packages before running Figure 3: ",
       paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(monocle)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

set.seed(20261018)
fig3_data_dir <- here::here("data", "fig3")
fig3_file <- function(name) file.path(fig3_data_dir, name)

malignant_cells <- readRDS(fig3_file("malignant_cells.rds"))
trajectory_cds <- readRDS(fig3_file("trajectory_cds.rds"))

trajectory_meta <- as.data.frame(pData(trajectory_cds), stringsAsFactors = FALSE)
trajectory_meta$cell.id <- if ("cell.id" %in% names(trajectory_meta)) {
  as.character(trajectory_meta$cell.id)
} else {
  rownames(trajectory_meta)
}
trajectory_meta$tissue_type <- gsub("[0-9]", "", trajectory_meta$orig.ident)
trajectory_meta$tissue_type <- ifelse(
  trajectory_meta$tissue_type %in% c("SN", "ST"), "GC",
  ifelse(trajectory_meta$tissue_type %in% c("PT", "PN"), "AGC",
         ifelse(trajectory_meta$tissue_type == "PTH", "ALi", "APB"))
)
trajectory_meta$diff_line <- ifelse(
  trajectory_meta$State %in% "1", "root-line",
  ifelse(trajectory_meta$State %in% as.character(2:6), "localized-line",
         ifelse(trajectory_meta$State %in% as.character(7:12), "invasive-line", "GC-line"))
)

# -------------------------------------------------------------------------
# Figure 3a: chromatin, histone, and DNA-methylation programs
# -------------------------------------------------------------------------
fig3a_genes <- list(
  CEBPA = "CEBPA",
  `Chromatin remodeling complex genes` = c("SMARCA4", "SMARCB1", "SMARCA2", "ARID1A", "ARID1B"),
  `Histone modification genes` = c("EP300", "CREBBP", "EZH2"),
  `DNA methylation genes` = c("DNMT3A", "DNMT3B", "DNMT1", "TET1", "TET2", "TET3")
)
all_fig3a_genes <- unique(unlist(fig3a_genes))
expression_matrix <- FetchData(malignant_cells, vars = intersect(all_fig3a_genes, rownames(malignant_cells)))
expression_matrix$cell.id <- rownames(expression_matrix)

score_column <- function(data, genes) {
  genes <- intersect(genes, names(data))
  if (length(genes) == 0) return(rep(0, nrow(data)))
  rowSums(data[, genes, drop = FALSE], na.rm = TRUE)
}
fig3a_scores <- data.frame(
  cell.id = expression_matrix$cell.id,
  CEBPA = if ("CEBPA" %in% names(expression_matrix)) expression_matrix$CEBPA else 0,
  `Chromatin remodeling complex genes` = score_column(expression_matrix, fig3a_genes[[2]]),
  `Histone modification genes` = score_column(expression_matrix, fig3a_genes[[3]]),
  `DNA methylation genes` = score_column(expression_matrix, fig3a_genes[[4]]),
  check.names = FALSE
)
fig3a_plot_data <- inner_join(
  trajectory_meta[, c("cell.id", "Pseudotime", "orig.ident", "tissue_type")],
  fig3a_scores, by = "cell.id"
) %>%
  filter(tissue_type != "GC", is.finite(Pseudotime))

fig3a_variables <- names(fig3a_genes)
fig3a_smooth <- bind_rows(lapply(fig3a_variables, function(variable) {
  fit <- loess(fig3a_plot_data[[variable]] ~ fig3a_plot_data$Pseudotime, span = 0.3)
  data.frame(
    Pseudotime = fig3a_plot_data$Pseudotime,
    Expression = as.numeric(predict(fit)),
    DisplayLabel = variable,
    stringsAsFactors = FALSE
  )
})) %>% arrange(DisplayLabel, Pseudotime)

fig3a_colors <- c(
  CEBPA = "#E64B35FF",
  `Chromatin remodeling complex genes` = "#00A087FF",
  `Histone modification genes` = "#4DBBD5FF",
  `DNA methylation genes` = "#F39B7FFF"
)
fig3a_linetypes <- setNames(rep("dotdash", length(fig3a_colors)), names(fig3a_colors))
fig3a_linetypes[["CEBPA"]] <- "solid"
fig3a <- ggplot(fig3a_smooth, aes(Pseudotime, Expression,
                                  color = DisplayLabel, linetype = DisplayLabel)) +
  geom_vline(xintercept = 2.34, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = 20.4, color = "#FFA500", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = 32.9, color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = fig3a_colors) +
  scale_linetype_manual(values = fig3a_linetypes) +
  labs(x = "Pseudotime", y = "Expression Level") +
  scale_y_continuous(limits = c(0, 4), breaks = 0:4, expand = c(0, 0)) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
        legend.position = "none")

# -------------------------------------------------------------------------
# Figure 3b-c: reported pseudotime transition points
# -------------------------------------------------------------------------
fig3b_c_cutoffs <- data.frame(
  transition = c("first transition", "AFP switch", "late transition"),
  pseudotime = c(2.34, 20.4, 32.9)
)

# -------------------------------------------------------------------------
# Figure 3g: AGC101 MYC pathway evolution panels
# -------------------------------------------------------------------------
read_gmt_simple <- function(path) {
  lines <- readLines(path, warn = FALSE)
  rows <- strsplit(lines, "\t", fixed = TRUE)
  bind_rows(lapply(rows, function(x) {
    if (length(x) < 3) return(NULL)
    data.frame(term = x[1], gene = x[-c(1, 2)], stringsAsFactors = FALSE)
  }))
}

infercnv_matrix <- readRDS(fig3_file("infercnv_gene_matrix.rds"))
hallmark <- read_gmt_simple(fig3_file("hallmark_pathways.gmt"))
hallmark_genes <- unique(hallmark$gene[hallmark$term == "HALLMARK_MYC_TARGETS_V1"])
myc_genes <- intersect(hallmark_genes, rownames(infercnv_matrix))
if (length(myc_genes) == 0) stop("No MYC pathway genes overlap the inferCNV matrix.")

ucell_scores <- read.table(fig3_file("ucell_scores.tsv"), sep = "\t", header = TRUE,
                           stringsAsFactors = FALSE, check.names = FALSE)
ucell_scores$cell.barcode <- rownames(ucell_scores)
malignant_meta <- malignant_cells@meta.data
malignant_meta$cell.barcode <- rownames(malignant_meta)
ucell_meta <- merge(malignant_meta, ucell_scores, by = "cell.barcode", all = FALSE)
ucell_meta <- ucell_meta[ucell_meta$orig.ident %in% c("P16N", "P16T", "P16TH"), , drop = FALSE]

evolution_order <- c(
  "root-line", "GC-line", "localized-line", "AFP- invasive-line",
  "AFP+ invasive-line", "AFP+ ALi invasive-line", "AFP- ALi invasive-line"
)
myc_colors <- c(
  "root-line" = "#EEEED1", "GC-line" = "#0079f4", "localized-line" = "#6acea6",
  "AFP- invasive-line" = "#FF7F00", "AFP+ invasive-line" = "#FFD700",
  "AFP+ ALi invasive-line" = "#3e8b26AD", "AFP- ALi invasive-line" = "#E64B35B2"
)
origin_colors <- c(Primary = "#4DBBD5B2", Metastasis = "#E64B35B2")

ucell_column <- "HALLMARK_MYC_TARGETS_V1_UCell"
if (!ucell_column %in% names(ucell_meta)) {
  stop("Missing UCell column: ", ucell_column)
}
score_vector <- colMeans(infercnv_matrix[myc_genes, , drop = FALSE])
infercnv_scores <- data.frame(
  cell.id = gsub("[.]", "-", names(score_vector)),
  infercnv.mean = as.numeric(score_vector),
  stringsAsFactors = FALSE
)
myc_ucell <- ucell_meta[, c("cell.barcode", ucell_column), drop = FALSE]
names(myc_ucell) <- c("cell.id", "UCell_score")
myc_meta <- trajectory_meta
myc_meta$AFP.type <- "others"
myc_meta$AFP.type[myc_meta$diff_line == "GC-line"] <- "GC-line"
myc_meta$AFP.type[myc_meta$diff_line == "root-line"] <- "root-line"
myc_meta$AFP.type[myc_meta$diff_line == "localized-line"] <- "localized-line"
myc_meta$AFP.type[myc_meta$diff_line == "invasive-line"] <- "AFP- invasive-line"
myc_meta$AFP.type[myc_meta$diff_line == "invasive-line" & myc_meta$Pseudotime > 20.4] <- "AFP+ invasive-line"
myc_meta$AFP.type[myc_meta$tissue_type == "ALi" & myc_meta$diff_line == "invasive-line"] <- "AFP- ALi invasive-line"
myc_meta$AFP.type[myc_meta$tissue_type == "ALi" & myc_meta$diff_line == "invasive-line" & myc_meta$Pseudotime > 20.4] <- "AFP+ ALi invasive-line"

myc_data <- myc_meta %>%
  merge(myc_ucell, by = "cell.id", all = FALSE) %>%
  merge(infercnv_scores, by = "cell.id", all = FALSE) %>%
  filter(AFP.type %in% evolution_order, is.finite(Pseudotime))
myc_data <- myc_data[!(myc_data$AFP.type == "localized-line" & myc_data$orig.ident == "P16TH"), , drop = FALSE]
myc_data <- myc_data[!(myc_data$AFP.type == "root-line" & myc_data$orig.ident == "P16TH"), , drop = FALSE]

safe_wilcox <- function(x, y) {
  if (length(x) < 3 || length(y) < 3) return(1)
  result <- tryCatch(wilcox.test(x, y)$p.value, error = function(e) 1)
  if (is.nan(result)) 1 else result
}
get_significance <- function(data, variable) {
  baseline <- data[[variable]][data$global_order == 1]
  max_value <- max(data[[variable]], na.rm = TRUE)
  offset <- (max_value - min(data[[variable]], na.rm = TRUE)) * 0.05
  data %>% group_by(global_order) %>%
    summarise(value_median = median(.data[[variable]], na.rm = TRUE),
              p_value = if (global_order[1] == 1) 1 else safe_wilcox(.data[[variable]], baseline),
              .groups = "drop") %>%
    filter(global_order != 1) %>%
    mutate(y_position = max_value + offset,
           label = case_when(p_value < 0.001 ~ "**", p_value < 0.05 ~ "*", TRUE ~ ""),
           sig_color = case_when(p_value >= 0.05 ~ NA_character_,
                                 value_median > median(baseline, na.rm = TRUE) ~ "red",
                                 value_median < median(baseline, na.rm = TRUE) ~ "blue",
                                 TRUE ~ NA_character_)) %>%
    filter(label != "")
}

build_myc_panels <- function(data) {
  data <- data %>% mutate(
    AFP.type = factor(AFP.type, levels = evolution_order),
    Origin = ifelse(AFP.type %in% evolution_order[1:5], "Primary", "Metastasis")
  ) %>% filter(!is.na(Origin))
  bins <- data %>% group_by(Origin) %>% group_split() %>% lapply(function(x) {
    n_bins <- ifelse(unique(x$Origin) == "Primary", 15, 5)
    x %>% mutate(bin = ntile(Pseudotime, min(n_bins, nrow(x))))
  }) %>% bind_rows()
  centroids <- bins %>% group_by(Origin, bin) %>%
    summarise(mean_X = mean(UCell_score, na.rm = TRUE),
              mean_Y = mean(infercnv.mean, na.rm = TRUE), cell_count = n(), .groups = "drop") %>%
    arrange(desc(Origin), bin) %>% mutate(global_order = row_number())
  final_bins <- bins %>% left_join(centroids %>% select(Origin, bin, global_order), by = c("Origin", "bin"))
  common_theme <- theme_classic() + theme(
    axis.text = element_text(size = 7, color = "black"), axis.title = element_text(size = 8),
    plot.title = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 6),
    legend.key.size = unit(0.2, "cm"), plot.margin = margin(t = 2, r = 5, b = 2, l = 5)
  )
  no_x_theme <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
  p1 <- ggplot() + geom_path(data = centroids, aes(mean_X, mean_Y), color = "black", linewidth = 0.6,
                             alpha = 0.7, arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
    geom_point(data = centroids, aes(mean_X, mean_Y, fill = Origin, size = cell_count), shape = 21,
               color = "white", stroke = 0.4, show.legend = FALSE) +
    geom_text_repel(data = centroids, aes(mean_X, mean_Y, label = global_order), size = 2.2,
                    max.overlaps = Inf) + scale_fill_manual(values = origin_colors) +
    scale_size_continuous(range = c(2, 4)) + scale_x_continuous(limits = c(0.1, 0.5), breaks = seq(0.1, 0.5, 0.1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.25, 0.45), breaks = seq(0.25, 0.45, 0.1), expand = c(0, 0)) +
    labs(x = "Transcriptional Activity", y = "Genomic Instability") + common_theme
  p2 <- ggplot(count(final_bins, global_order, AFP.type), aes(factor(global_order), n, fill = AFP.type)) +
    geom_col(width = 0.6) + scale_fill_manual(values = myc_colors) + scale_y_continuous(limits = c(0, 90), breaks = seq(0, 90, 30), expand = c(0, 0)) +
    labs(x = NULL, y = "Count") + common_theme + no_x_theme
  p3 <- ggplot(final_bins, aes(factor(global_order), infercnv.mean)) + geom_boxplot(aes(fill = Origin), width = 0.6, outlier.size = 0.1, linewidth = 0.2, fatten = 1) +
    scale_fill_manual(values = origin_colors) + scale_y_continuous(limits = c(0.2, 0.6), breaks = seq(0.2, 0.6, 0.2), expand = c(0, 0)) + labs(x = NULL, y = "inferCNV") + common_theme + no_x_theme + theme(legend.position = "none")
  sig3 <- get_significance(final_bins, "infercnv.mean")
  if (nrow(sig3) > 0) p3 <- p3 + geom_text(data = sig3, aes(factor(global_order), y_position, label = label, color = sig_color), inherit.aes = FALSE, size = 3.5, show.legend = FALSE) + scale_color_identity()
  p4 <- ggplot(final_bins, aes(factor(global_order), UCell_score)) + geom_boxplot(aes(fill = Origin), width = 0.6, outlier.size = 0.1, linewidth = 0.2, fatten = 1) +
    scale_fill_manual(values = origin_colors) + scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2), expand = c(0, 0)) + labs(x = NULL, y = "UCell") + common_theme + no_x_theme + theme(legend.position = "none")
  sig4 <- get_significance(final_bins, "UCell_score")
  if (nrow(sig4) > 0) p4 <- p4 + geom_text(data = sig4, aes(factor(global_order), 0.55, label = label, color = sig_color), inherit.aes = FALSE, size = 3.5, show.legend = FALSE) + scale_color_identity()
  p5 <- ggplot(final_bins, aes(factor(global_order), Pseudotime, fill = Origin)) + geom_boxplot(width = 0.6, outlier.size = 0.1, linewidth = 0.2, fatten = 1) +
    scale_fill_manual(values = origin_colors) + scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10), expand = c(0, 0)) + labs(x = "Bin Index", y = "Pseudotime") + common_theme + theme(legend.position = "none")
  list(plot = (p1 / p2 / p3 / p4 / p5) + plot_layout(heights = c(2.5, 1, 1, 1, 1), guides = "collect") & theme(legend.position = "bottom"), bins = final_bins, centroids = centroids)
}

fig3g_parts <- build_myc_panels(myc_data)
fig3g <- fig3g_parts$plot
fig3g_source_bins <- fig3g_parts$bins
fig3g_source_centroids <- fig3g_parts$centroids

fig3a
fig3b_c_cutoffs
fig3g
