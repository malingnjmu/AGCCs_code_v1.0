# Figure 1: integrated analysis and plotting script
required_packages <- c(
  "here", "ggplot2", "ggrepel", "dplyr", "data.table", "patchwork",
  "maftools", "Seurat", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install the missing packages before running Figure 1: ",
       paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(data.table)
  library(patchwork)
  library(maftools)
  library(Seurat)
  library(scales)
})

set.seed(20261018)

fig1_data_dir <- here::here("data", "fig1")
fig1_figure_dir <- here::here("figures")
fig1_source_data_dir <- here::here("source_data")

fig1_file <- function(name) file.path(fig1_data_dir, name)

fig1_colors <- c(ALi = "#eb1d37", AGC = "#FFA500", GC = "#0079f4")
fig1_groups <- c("ALi", "AGC", "GC")
fig1_theme <- theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title.y = element_text(color = "black", size = 10, margin = margin(r = 5)),
    axis.title.x = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    legend.position = "none",
    panel.grid = element_blank()
  )

# -------------------------------------------------------------------------
# Figure 1c: mutation prevalence in in-house GC versus AGC
# -------------------------------------------------------------------------
agc_clinical <- read.table(fig1_file("agc_clinical.tsv"), header = TRUE,
                           stringsAsFactors = FALSE, sep = "\t")
agc_maf <- read.maf(
  fig1_file("agc_mutations.maf"),
  vc_nonSyn = c("Missense", "Nonsense", "Splice site", "Frameshift", "In-frame indel"),
  clinicalData = agc_clinical
)
gc_clinical <- read.table(fig1_file("gc_clinical.tsv"), header = TRUE,
                          stringsAsFactors = FALSE, sep = "\t")
gc_maf <- read.maf(
  fig1_file("gc_mutations.maf"),
  vc_nonSyn = c("Missense", "Nonsense", "Splice site", "Frameshift", "In-frame indel"),
  clinicalData = gc_clinical
)

cancer_genes_data <- read.table(fig1_file("top_cancer_genes.tsv"), sep = "\t",
                                header = TRUE, stringsAsFactors = FALSE)
cancer_genes <- cancer_genes_data$Gene

get_prevalence <- function(maf_object, genes) {
  summary_data <- maftools::getGeneSummary(maf_object)
  summary_data <- summary_data[Hugo_Symbol %in% genes,
                               .(Hugo_Symbol, MutatedSamples)]
  missing_genes <- setdiff(genes, summary_data$Hugo_Symbol)
  if (length(missing_genes) > 0) {
    summary_data <- rbind(
      summary_data,
      data.table(Hugo_Symbol = missing_genes, MutatedSamples = 0)
    )
  }
  total_samples <- as.numeric(
    maf_object@summary$summary[maf_object@summary$ID == "Samples"]
  )
  summary_data[, Prevalence := MutatedSamples / total_samples * 100]
  list(data = summary_data, total = total_samples)
}

agc_result <- get_prevalence(agc_maf, cancer_genes)
gc_result <- get_prevalence(gc_maf, cancer_genes)
fig1c_data <- data.frame(
  Gene = cancer_genes,
  AGC_Freq = agc_result$data$Prevalence[match(cancer_genes, agc_result$data$Hugo_Symbol)],
  GC_Freq = gc_result$data$Prevalence[match(cancer_genes, gc_result$data$Hugo_Symbol)],
  AGC_Mut = agc_result$data$MutatedSamples[match(cancer_genes, agc_result$data$Hugo_Symbol)],
  GC_Mut = gc_result$data$MutatedSamples[match(cancer_genes, gc_result$data$Hugo_Symbol)]
)
fig1c_data$p_value <- vapply(seq_len(nrow(fig1c_data)), function(i) {
  contingency <- matrix(c(
    fig1c_data$AGC_Mut[i], agc_result$total - fig1c_data$AGC_Mut[i],
    fig1c_data$GC_Mut[i], gc_result$total - fig1c_data$GC_Mut[i]
  ), nrow = 2)
  fisher.test(contingency)$p.value
}, numeric(1))
fig1c_data$fdr <- p.adjust(fig1c_data$p_value, method = "fdr")
fig1c_data$Significance <- ifelse(fig1c_data$fdr < 0.05, "FDR < 0.05", "FDR >= 0.05")

fig1c_colors <- c("FDR < 0.05" = "#D62728", "FDR >= 0.05" = "#1F77B4")
fig1c_theme <- theme_bw() + theme(
  panel.grid = element_blank(),
  panel.border = element_rect(linewidth = 0.6, color = "black"),
  axis.ticks = element_line(linewidth = 0.4, color = "black"),
  axis.text = element_text(size = 15, color = "black"),
  axis.title = element_text(size = 15, color = "black"),
  legend.text = element_text(size = 10),
  legend.title = element_blank(),
  aspect.ratio = 1,
  axis.title.x = element_blank()
)
fig1c_main <- ggplot(fig1c_data, aes(GC_Freq, AGC_Freq)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_point(aes(color = Significance), size = 1.3) +
  geom_text_repel(aes(label = ifelse(fdr < 0.05, Gene, "")), size = 3,
                  max.overlaps = Inf, min.segment.length = 0,
                  segment.color = "grey50", segment.size = 0.3, box.padding = 0.5) +
  annotate("rect", xmin = -0.1, xmax = 20, ymin = -0.1, ymax = 20,
           fill = NA, color = "black", linetype = "dotted", linewidth = 0.5) +
  scale_color_manual(values = fig1c_colors) +
  scale_x_continuous(limits = c(-0.5, 80), breaks = seq(0, 80, 20), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.5, 80), breaks = seq(0, 80, 20), expand = c(0, 0)) +
  labs(y = "Prevalence in AGC (%)") + fig1c_theme +
  theme(legend.position = c(0.15, 0.93), legend.background = element_blank())
fig1c_zoom <- fig1c_data %>% filter(GC_Freq <= 20, AGC_Freq <= 20)
fig1c_inset <- ggplot(fig1c_zoom, aes(GC_Freq, AGC_Freq)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_point(aes(color = Significance), size = 2) +
  geom_text_repel(aes(label = Gene), size = 3, max.overlaps = Inf,
                  min.segment.length = 0, segment.color = "black", segment.size = 0.35,
                  box.padding = 0.6, point.padding = 0.3, force = 10, direction = "both") +
  scale_color_manual(values = fig1c_colors) +
  scale_x_continuous(limits = c(-0.5, 20), breaks = seq(0, 20, 5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.5, 20), breaks = seq(0, 20, 5), expand = c(0, 0)) +
  fig1c_theme + theme(axis.title.y = element_blank(), legend.position = "none")
fig1c <- fig1c_main | fig1c_inset

# -------------------------------------------------------------------------
# Figure 1d/e: tumor mutation burden and SCNA burden
# -------------------------------------------------------------------------
tmb_data <- read.table(fig1_file("tmb_summary.tsv"), sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE)
scna_data <- read.table(fig1_file("scna_summary.tsv"), sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE, check.names = FALSE)
tmb_data$TMB <- tmb_data$mutation_counts / 38
scna_data$Type <- factor(scna_data$Type, levels = fig1_groups)
group_colors <- fig1_colors

fig1d_data <- subset(tmb_data, mutation_counts < 600)
fig1d_data$Type <- factor(fig1d_data$Type, levels = fig1_groups)
fig1d <- ggplot(fig1d_data, aes(Type, TMB, fill = Type)) +
  geom_jitter(aes(color = Type), width = 0.25, size = 1.2, alpha = 0.6, stroke = 0) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, color = "black", linewidth = 0.6) +
  scale_fill_manual(values = group_colors) + scale_color_manual(values = group_colors) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.05, 20), breaks = seq(0, 20, 5)) +
  labs(y = "Tumor mutation burden") + fig1_theme

fig1e <- ggplot(scna_data, aes(Type, Value, fill = Type)) +
  geom_jitter(aes(color = Type), width = 0.25, size = 1.2, alpha = 0.6, stroke = 0) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, color = "black", linewidth = 0.6) +
  scale_fill_manual(values = group_colors) + scale_color_manual(values = group_colors) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                expand = c(0, 0), limits = c(1, 10000)) +
  labs(y = "Burden of SCNAs") + fig1_theme

# Statistical summaries formerly split across two auxiliary scripts.
fig1_wilcox <- function(data, value, figure_name, filter_high_mutations = FALSE) {
  if (filter_high_mutations && "mutation_counts" %in% names(data)) {
    data <- data[data$mutation_counts < 600, , drop = FALSE]
  }
  a <- data[[value]][data$Type == "ALi"]
  b <- data[[value]][data$Type == "AGC"]
  c <- data[[value]][data$Type == "GC"]
  data.frame(
    Figure = figure_name,
    Comparison = c("ALi vs AGC", "AGC vs GC"),
    ALi_n = c(length(a), NA), AGC_n = c(length(b), length(b)), GC_n = c(NA, length(c)),
    ALi_median = c(median(a), NA), AGC_median = c(median(b), median(b)), GC_median = c(NA, median(c)),
    Wilcoxon_p = c(wilcox.test(a, b, exact = FALSE)$p.value,
                   wilcox.test(b, c, exact = FALSE)$p.value),
    check.names = FALSE
  )
}
fig1_statistics <- rbind(
  fig1_wilcox(tmb_data, "TMB", "Figure 1d"),
  fig1_wilcox(scna_data, "Value", "Figure 1e")
)

# -------------------------------------------------------------------------
# Figure 1h: inferCNV and CytoTRACE by cancer type
# -------------------------------------------------------------------------
malignant_object <- readRDS(fig1_file("malignant_cells.rds"))
infercnv_scores <- readRDS(fig1_file("infercnv_scores.rds"))
malignant_meta <- malignant_object@meta.data
malignant_meta$inferCNV_normlized <- NULL
malignant_meta <- merge(malignant_meta, infercnv_scores, by = "cell.id", all = FALSE)

fig1h_left_data <- data.frame(
  Cell_ID = malignant_meta$cell.id,
  `Cancer type` = factor(malignant_meta$tissue_type, levels = fig1_groups),
  inferCNV = malignant_meta$inferCNV_normlized,
  check.names = FALSE
)
fig1h_left <- ggplot(fig1h_left_data, aes(`Cancer type`, inferCNV, fill = `Cancer type`)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, color = "black", linewidth = 0.6) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(limits = c(0.2, 0.7), breaks = seq(0.2, 0.7, 0.1), expand = c(0, 0)) +
  labs(x = NULL, y = "inferCNV") + fig1_theme

fig1h_right_data <- data.frame(
  Cell_ID = rownames(malignant_object@meta.data),
  `Cancer type` = factor(malignant_object@meta.data$tissue_type, levels = fig1_groups),
  CytoTRACE = malignant_object@meta.data$CytoTRACE,
  check.names = FALSE
)
fig1h_right <- ggplot(fig1h_right_data, aes(`Cancer type`, CytoTRACE, fill = `Cancer type`)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, color = "black", linewidth = 0.6) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = NULL, y = "CytoTRACE") + fig1_theme

# -------------------------------------------------------------------------
# Figure 1i: SCNA burden versus CytoTRACE / inferCNV
# -------------------------------------------------------------------------
scna_cytotrace <- read.table(fig1_file("sample_cytotrace_scna.tsv"), sep = "\t",
                             header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
scna_cytotrace$`Cancer type` <- scna_cytotrace$Type
fig1i_left_data <- scna_cytotrace[, c("WES_ID", "orig.ident", "Cancer type", "CytoTRACE.score", "Value")]
names(fig1i_left_data) <- c("WES_ID", "scRNA_sample", "Cancer type", "Mean_CytoTRACE_score", "SCNA_burden")
fig1i_left_cor <- cor.test(fig1i_left_data$Mean_CytoTRACE_score,
                           fig1i_left_data$SCNA_burden, method = "spearman")
fig1i_left_label <- paste0("R = ", format(round(fig1i_left_cor$estimate, 2), nsmall = 2),
                           "\n", if (fig1i_left_cor$p.value < 0.001) "P < 0.001" else
                             paste0("P = ", format(fig1i_left_cor$p.value, digits = 3)))
fig1i_left <- ggplot(fig1i_left_data, aes(Mean_CytoTRACE_score, SCNA_burden)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +
  geom_point(aes(color = `Cancer type`), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = group_colors) +
  annotate("text", x = 0.02, y = Inf, label = fig1i_left_label, hjust = 0, vjust = 1.2, size = 3) +
  scale_x_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.3), expand = c(0, 0)) +
  scale_y_log10(limits = c(1, 10000), breaks = c(1, 10, 100, 1000, 10000),
                labels = scales::trans_format("log10", scales::math_format(10^.x)), expand = c(0, 0)) +
  labs(x = "Mean CytoTRACE score", y = "Burden of SCNAs") +
  theme_classic(base_size = 8) + theme(axis.text = element_text(color = "black", size = 7),
                                       axis.title = element_text(size = 8),
                                       legend.position = "none")

sample_mapping <- read.table(fig1_file("sample_mapping.tsv"), sep = "\t", header = TRUE,
                             stringsAsFactors = FALSE, check.names = FALSE)
fig1i_right_meta <- merge(malignant_meta, infercnv_scores, by = "cell.id", all = FALSE)
fig1i_right_data <- fig1i_right_meta %>%
  filter(tissue_type != "GC") %>%
  group_by(orig.ident) %>%
  summarise(
    CytoTRACE = mean(CytoTRACE, na.rm = TRUE),
    Mean_inferCNV_score = mean(inferCNV_normlized, na.rm = TRUE),
    `Cancer type` = first(tissue_type), Cell_count = n(), .groups = "drop"
  ) %>%
  left_join(sample_mapping, by = "orig.ident")
fig1i_right_cor <- cor.test(fig1i_right_data$CytoTRACE,
                            fig1i_right_data$Mean_inferCNV_score, method = "spearman")
fig1i_right_label <- paste0("R = ", format(round(fig1i_right_cor$estimate, 2), nsmall = 2),
                            "\n", if (fig1i_right_cor$p.value < 0.001) "P < 0.001" else
                              paste0("P = ", format(fig1i_right_cor$p.value, digits = 3)))
fig1i_right <- ggplot(fig1i_right_data, aes(CytoTRACE, Mean_inferCNV_score)) +
  geom_smooth(method = "lm", color = "black", se = FALSE, linewidth = 0.5) +
  geom_point(aes(color = `Cancer type`), size = 1.5, alpha = 0.8, stroke = 0) +
  scale_color_manual(values = group_colors) +
  annotate("text", x = 0.02, y = Inf, label = fig1i_right_label, hjust = 0, vjust = 1.2, size = 3) +
  scale_x_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.3), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.25, 0.6), breaks = seq(0.25, 0.55, 0.1), expand = c(0, 0)) +
  labs(x = NULL, y = "Mean inferCNV score") +
  theme_classic(base_size = 8) + theme(axis.text = element_text(color = "black", size = 7),
                                       axis.title = element_text(size = 8),
                                       legend.position = "none")

# Objects are intentionally left in the global environment for inspection.
fig1c
fig1d
fig1e
fig1h_left
fig1h_right
fig1i_left
fig1i_right
