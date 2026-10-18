# Figure 2: integrated trajectory, AFP-switch, and marker-expression analysis
#
# Expected repository layout:
#   data/fig2/
#   figures/       optional output directory
#   source_data/   optional source-data directory
#
# The script defines fig2b_data, fig2c, and fig2d. It does not force a
# working directory, install packages, or save figures automatically.

required_packages <- c(
  "here", "Seurat", "monocle", "dplyr", "tidyr", "ggplot2", "patchwork"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install the missing packages before running Figure 2: ",
       paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(monocle)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(20261018)
fig2_data_dir <- here::here("data", "fig2")
fig2_file <- function(name) file.path(fig2_data_dir, name)

malignant_cells <- readRDS(fig2_file("malignant_cells.rds"))
trajectory_cds <- readRDS(fig2_file("trajectory_cds.rds"))

cell_meta <- malignant_cells@meta.data
cell_meta$diff_line <- ifelse(
  cell_meta$State %in% "1", "root-line",
  ifelse(cell_meta$State %in% as.character(2:6), "localized-line",
         ifelse(cell_meta$State %in% as.character(7:12), "invasive-line", "GC-line"))
)

afp_expression <- FetchData(malignant_cells, vars = c("AFP", "FN1"))
afp_expression$cell.id <- rownames(afp_expression)
cell_meta$cell.id <- rownames(cell_meta)

# -------------------------------------------------------------------------
# Figure 2b: pseudotime groups and AFP-positive invasive cells
# -------------------------------------------------------------------------
trajectory_meta <- as.data.frame(pData(trajectory_cds), stringsAsFactors = FALSE)
trajectory_meta$cell.id <- if ("cell.id" %in% names(trajectory_meta)) {
  as.character(trajectory_meta$cell.id)
} else {
  rownames(trajectory_meta)
}
trajectory_meta$tissue_type <- gsub("[0-9]", "", trajectory_meta$orig.ident)
trajectory_meta$tissue_type <- dplyr::case_when(
  trajectory_meta$tissue_type %in% c("SN", "ST", "PN", "PT") ~ "AGC",
  trajectory_meta$tissue_type == "PTH" ~ "ALi",
  TRUE ~ "GC"
)
trajectory_meta$diff_line <- dplyr::case_when(
  as.character(trajectory_meta$State) == "1" ~ "root-line",
  as.character(trajectory_meta$State) %in% as.character(2:6) ~ "localized-line",
  as.character(trajectory_meta$State) %in% as.character(7:12) ~ "invasive-line",
  TRUE ~ "GC-line"
)
trajectory_meta <- left_join(
  trajectory_meta,
  afp_expression[, c("cell.id", "AFP")],
  by = "cell.id"
)
trajectory_meta$AFP.type <- "others"
trajectory_meta$AFP.type[trajectory_meta$diff_line == "GC-line"] <- "GC-line"
trajectory_meta$AFP.type[trajectory_meta$diff_line == "root-line"] <- "root-line"
trajectory_meta$AFP.type[trajectory_meta$diff_line == "invasive-line"] <- "AFP- invasive-line"
trajectory_meta$AFP.type[trajectory_meta$diff_line == "localized-line"] <- "localized-line"
trajectory_meta$AFP.type[trajectory_meta$AFP > 0 & trajectory_meta$diff_line == "invasive-line"] <- "AFP+ invasive-line"
trajectory_meta$AFP.type[trajectory_meta$tissue_type == "ALi" & trajectory_meta$diff_line == "invasive-line"] <- "ALi invasive-line"

fig2b_levels <- c(
  "root-line", "GC-line", "localized-line", "AFP- invasive-line",
  "AFP+ invasive-line", "ALi invasive-line"
)
fig2b_data <- trajectory_meta %>%
  filter(AFP.type %in% fig2b_levels, !is.na(Pseudotime)) %>%
  transmute(
    cell_id = cell.id,
    pseudotime = as.numeric(Pseudotime),
    group = factor(AFP.type, levels = fig2b_levels),
    tissue_type,
    state = as.character(State),
    diff_line,
    AFP_positive = AFP > 0
  )

# -------------------------------------------------------------------------
# Figure 2c: DEGs shared with the ALi invasive signature
# -------------------------------------------------------------------------
afp_switch_cutoff <- 20.33
afp_positive <- afp_expression$AFP > 0
names(afp_positive) <- afp_expression$cell.id
cell_meta$AFP_diff_line <- "others"
cell_meta$AFP_diff_line[cell_meta$diff_line == "GC-line"] <- "GC-line"
cell_meta$AFP_diff_line[cell_meta$diff_line == "root-line"] <- "root-line"
cell_meta$AFP_diff_line[cell_meta$diff_line == "localized-line"] <- "localized-line"
cell_meta$AFP_diff_line[cell_meta$diff_line == "invasive-line"] <- "before invasive-line"
cell_meta$AFP_diff_line[cell_meta$diff_line == "invasive-line" &
                         cell_meta$Pseudotime > afp_switch_cutoff] <- "after invasive-line"
cell_meta$AFP_diff_line[cell_meta$tissue_type == "ALi" &
                         cell_meta$diff_line == "invasive-line"] <- "ALi invasive-line"
malignant_cells@meta.data$AFP_diff_line <- cell_meta$AFP_diff_line
Idents(malignant_cells) <- "AFP_diff_line"

marker_for <- function(group) {
  result <- FindMarkers(
    malignant_cells, ident.1 = group, ident.2 = "root-line",
    only.pos = TRUE, logfc.threshold = 1
  )
  result$gene <- rownames(result)
  result[result$p_val_adj < 0.01, , drop = FALSE]
}

fig2c_markers <- list(
  `GC-line` = marker_for("GC-line"),
  `localized-line` = marker_for("localized-line"),
  `before invasive-line` = marker_for("before invasive-line"),
  `after invasive-line` = marker_for("after invasive-line"),
  `ALi invasive-line` = marker_for("ALi invasive-line")
)
ali_genes <- fig2c_markers[["ALi invasive-line"]]$gene
fig2c_groups <- c("GC-line", "localized-line", "after invasive-line", "ALi invasive-line")
fig2c_plot_data <- tibble(
  diff_line = factor(fig2c_groups, levels = fig2c_groups),
  Up_DEGs = vapply(fig2c_groups, function(x) nrow(fig2c_markers[[x]]), integer(1))
) %>%
  mutate(
    share_percentage = round(100 * vapply(fig2c_groups, function(x) {
      length(intersect(fig2c_markers[[x]]$gene, ali_genes)) / length(ali_genes)
    }, numeric(1)), 1),
    shared_pct = vapply(fig2c_groups, function(x) {
      length(intersect(fig2c_markers[[x]]$gene, ali_genes)) /
        nrow(fig2c_markers[[x]])
    }, numeric(1))
  ) %>%
  mutate(
    shared = round(Up_DEGs * shared_pct, 0),
    not_shared = Up_DEGs - shared
  ) %>%
  pivot_longer(c(shared, not_shared), names_to = "gene_type", values_to = "gene_count")

fig2c_scale_factor <- 6000 / 100
fig2c <- ggplot(fig2c_plot_data, aes(x = diff_line)) +
  geom_col(aes(y = gene_count, fill = gene_type),
           position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.2) +
  geom_line(aes(y = share_percentage * fig2c_scale_factor, group = 1),
            color = "#4D4D4D", linewidth = 0.6) +
  geom_point(aes(y = share_percentage * fig2c_scale_factor), color = "#B2182B",
             fill = "white", size = 2.5, stroke = 1, shape = 21) +
  geom_text(aes(y = gene_count, label = gene_count, group = gene_type),
            position = position_dodge(width = 0.7), vjust = -0.1, size = 3,
            color = "black") +
  geom_text(data = fig2c_plot_data %>% distinct(diff_line, share_percentage),
            aes(y = share_percentage * fig2c_scale_factor,
                label = paste0(share_percentage, "%")),
            nudge_y = 250, vjust = 0.2, size = 3.5, color = "black") +
  scale_y_continuous(
    name = "No. of Up-regulated Genes (vs. root-line)", limits = c(0, 6500),
    breaks = seq(0, 6000, 2000), expand = c(0, 0),
    sec.axis = sec_axis(~ . / fig2c_scale_factor,
                        name = "Overlap with ALi signature (%)")
  ) +
  scale_fill_manual(
    breaks = c("shared", "not_shared"),
    values = c(shared = "#B2182B", not_shared = "#E0E0E0"),
    labels = c("Shared Genes", "Not-shared Genes")
  ) +
  labs(x = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_text(size = 10, face = "plain", color = "black"),
    axis.title.y.right = element_text(size = 10, face = "plain", color = "black", angle = 90),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.y.right = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1),
    panel.border = element_rect(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.6, color = "black"),
    legend.position = c(0.12, 0.92),
    legend.text = element_text(size = 8, color = "black"),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key.size = unit(4, "mm"),
    plot.margin = margin(10, 10, 5, 10, "mm")
  )

# -------------------------------------------------------------------------
# Figure 2d: marker expression along pseudotime
# -------------------------------------------------------------------------
fig2d_genes <- c("AFP", "APOA1", "VTN", "GKN1", "MUC5AC", "TFF2")
fig2d_expression <- FetchData(malignant_cells, vars = fig2d_genes)
fig2d_expression$cell_id <- rownames(fig2d_expression)
fig2d_data <- inner_join(
  trajectory_meta[, c("cell.id", "Pseudotime")],
  fig2d_expression, by = c("cell.id" = "cell_id")
) %>%
  filter(!is.na(Pseudotime)) %>%
  pivot_longer(all_of(fig2d_genes), names_to = "gene", values_to = "expression") %>%
  group_by(gene) %>%
  arrange(Pseudotime, .by_group = TRUE) %>%
  mutate(loess_fitted = as.numeric(predict(loess(expression ~ Pseudotime, data = cur_data(), span = 0.4)))) %>%
  ungroup() %>%
  transmute(
    cell_id, gene, pseudotime = as.numeric(Pseudotime),
    expression = as.numeric(expression), loess_fitted
  )

fig2d_plots <- lapply(fig2d_genes, function(gene_name) {
  gene_data <- fig2d_data %>% filter(gene == gene_name)
  ggplot(gene_data, aes(pseudotime, expression)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_line(aes(y = loess_fitted),
              color = ifelse(gene_name %in% c("AFP", "APOA1", "VTN"), "#941a1a", "#0079f4"),
              linewidth = 1.1) +
    geom_vline(xintercept = afp_switch_cutoff, color = "#FFA500", linetype = "dashed") +
    labs(title = gene_name, x = "Pseudotime", y = "Expression") +
    theme_classic()
})
names(fig2d_plots) <- fig2d_genes
fig2d <- patchwork::wrap_plots(fig2d_plots, nrow = 2)

fig2b_data
fig2c
fig2d
