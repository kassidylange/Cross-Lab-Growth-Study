##code for image analyses to be used in R studio

```{r}
# install packages
install.packages("readxl")
install.packages("ggplot2")
install.packages("car")
install.packages("dplyr")
install.packages("ggpubr")
install.packages("patchwork")
```

```{r}
# pull in libs
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(patchwork)
  library(car)

# update paths & output directory

data_path <- "CLGE data.xlsx"
out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)

# load and prep data

data <- read_excel(data_path, sheet = "data")
colnames(data) <- c("Lab", "Species", "Colony_ID", "Growth_cm", "Polyps_change")

data <- data %>%
  mutate(
    Lab = factor(Lab),
    Species = factor(Species)
  )


# test assumptions + run statistical analyses based on assumptions (parametric vs nonparametric)

message("\n--- assumption checks: Growth_cm (Species × Lab) ---")
model_g <- aov(Growth_cm ~ Species * Lab, data = data)
shapiro_g <- shapiro.test(residuals(model_g))
levene_g  <- car::leveneTest(Growth_cm ~ Species * Lab, data = data)
print(shapiro_g); print(levene_g)

if (shapiro_g$p.value > 0.05 & levene_g$`Pr(>F)`[1] > 0.05) {
  message("Growth_cm: assumptions met- complete a two-way ANOVA")
  print(summary(model_g))
  print(TukeyHSD(model_g))
} else {
  message("Growth_cm: assumptions violated- complete kruskal–wallis + pairwise wilcoxon") #assumptions were violated for growth data
  print(kruskal.test(Growth_cm ~ Species, data = data))
  print(kruskal.test(Growth_cm ~ Lab, data = data))
  print(pairwise.wilcox.test(data$Growth_cm, interaction(data$Species, data$Lab), p.adjust.method = "BH")) #significant differences found between SPECIES but not between LABS
}

message("\n--- assumption checks: Polyps_change (Species × Lab) ---")
model_p <- aov(Polyps_change ~ Species * Lab, data = data)
shapiro_p <- shapiro.test(residuals(model_p))
levene_p  <- car::leveneTest(Polyps_change ~ Species * Lab, data = data)
print(shapiro_p); print(levene_p)

if (shapiro_p$p.value > 0.05 & levene_p$`Pr(>F)`[1] > 0.05) {
  message("Polyps_change: assumptions met- complete two-way ANOVA") #assumptions were met for percent change in polyp count
  print(summary(model_p)) #significant difference found between SPECIES but not between LABS
  print(TukeyHSD(model_p))
} else {
  message("Polyps_change: assumptions violated- complete kruskal–wallis + pairwise wilcoxon")
  print(kruskal.test(Polyps_change ~ Species, data = data))
  print(kruskal.test(Polyps_change ~ Lab, data = data))
  print(pairwise.wilcox.test(data$Polyps_change, interaction(data$Species, data$Lab), p.adjust.method = "BH"))
}

# assign color palettes (shades per lab) by species for pretty plots

lab_colors <- list(
  "Muricea pendula" = c("GAL" = "grey30", "GNV" = "grey60", "HML" = "grey80"),
  "Swiftia exserta" = c("GAL" = "#cc7a00", "GNV" = "#e69f00", "HML" = "#ffb84d"),
  "Thesea nivea"          = c("GAL" = "#6a51a3", "GNV" = "#9467BD", "HML" = "#c5b0d5")
)

# establish shared y-limits to standardize plots for easy comparisons

y_lim_growth <- range(data$Growth_cm, na.rm = TRUE)
y_lim_polyps <- range(data$Polyps_change, na.rm = TRUE)

# make merged plots with significance indicators and biological thresholds for both datasets

make_species_plot_growth <- function(sp) {
  sp_data <- data %>% filter(Species == sp, !is.na(Growth_cm)) %>% droplevels()
  p <- ggplot(sp_data, aes(Lab, Growth_cm, fill = Lab)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.6) +
    geom_boxplot(color = "black", outlier.shape = 16, outlier.size = 1.5) +
    scale_fill_manual(values = lab_colors[[sp]]) +
    coord_cartesian(ylim = y_lim_growth) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.major.x = element_blank(),
          legend.title = element_blank(),
          legend.position = "right") +
    labs(title = bquote(italic(.(sp))), x = NULL, y = NULL)

  if (dplyr::n_distinct(sp_data$Lab) < 2) return(p + labs(subtitle = "Only one lab has data"))

  pw <- try({
    sp_data %>%
      rstatix::pairwise_wilcox_test(Growth_cm ~ Lab, p.adjust.method = "BH") %>%
      rstatix::add_xy_position(x = "Lab", data = sp_data) %>%
      dplyr::mutate(
        p.signif = dplyr::case_when(
          p.adj <= 0.001 ~ "***",
          p.adj <= 0.01  ~ "**",
          p.adj <= 0.05  ~ "*",
          TRUE ~ "ns"
        )
      ) %>% dplyr::filter(!is.na(p.adj))
  }, silent = TRUE)
  if (inherits(pw, "try-error") || nrow(pw) == 0) return(p)

  p + ggpubr::stat_pvalue_manual(pw, label = "p.signif", tip.length = 0.01, hide.ns = TRUE)
}

make_species_plot_polyps <- function(sp) {
  sp_data <- data %>% filter(Species == sp, !is.na(Polyps_change)) %>% droplevels()
  p <- ggplot(sp_data, aes(Lab, Polyps_change, fill = Lab)) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "red", linewidth = 0.6) +
    geom_boxplot(color = "black", outlier.shape = 16, outlier.size = 1.5) +
    scale_fill_manual(values = lab_colors[[sp]]) +
    coord_cartesian(ylim = y_lim_polyps) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.major.x = element_blank(),
          legend.title = element_blank(),
          legend.position = "right") +
    labs(title = bquote(italic(.(sp))), x = NULL, y = NULL)

  if (dplyr::n_distinct(sp_data$Lab) < 2) return(p + labs(subtitle = "Only one lab has data"))

  pw <- try({
    sp_data %>%
      rstatix::pairwise_wilcox_test(Polyps_change ~ Lab, p.adjust.method = "BH") %>%
      rstatix::add_xy_position(x = "Lab", data = sp_data) %>%
      dplyr::mutate(
        p.signif = dplyr::case_when(
          p.adj <= 0.001 ~ "***",
          p.adj <= 0.01  ~ "**",
          p.adj <= 0.05  ~ "*",
          TRUE ~ "ns"
        )
      ) %>% dplyr::filter(!is.na(p.adj))
  }, silent = TRUE)
  if (inherits(pw, "try-error") || nrow(pw) == 0) return(p)

  p + ggpubr::stat_pvalue_manual(pw, label = "p.signif", tip.length = 0.01, hide.ns = TRUE)
}


# output individual species plots

p_muricea_g <- make_species_plot_growth("Muricea pendula")
p_swiftia_g <- make_species_plot_growth("Swiftia exserta")
p_thesea_g  <- make_species_plot_growth("Thesea nivea")  

p_muricea_p <- make_species_plot_polyps("Muricea pendula")
p_swiftia_p <- make_species_plot_polyps("Swiftia exserta")
p_thesea_p  <- make_species_plot_polyps("Thesea nivea")


# combine with shared legend + axis labels using ggpubr

row_g <- ggpubr::ggarrange(p_muricea_g, p_swiftia_g, p_thesea_g,
                           ncol = 3, nrow = 1,
                           align = "v",
                           common.legend = TRUE, legend = "right")

final_row_g <- ggpubr::annotate_figure(
  row_g,
  left   = ggpubr::text_grob("Growth Rate (cm/branch/yr)", rot = 90, size = 12),
  bottom = ggpubr::text_grob("Lab", size = 12)
)

row_p <- ggpubr::ggarrange(p_muricea_p, p_swiftia_p, p_thesea_p,
                           ncol = 3, nrow = 1,
                           align = "v",
                           common.legend = TRUE, legend = "right")

final_row_p <- ggpubr::annotate_figure(
  row_p,
  left   = ggpubr::text_grob("%\u0394 Polyps over 6 mo.", rot = 90, size = 12),
  bottom = ggpubr::text_grob("Lab", size = 12)
)


# export high-res figures to output folder in R

ggsave(file.path(out_dir, "Growth_three_species_row.png"),
       final_row_g, width = 12, height = 4.5, dpi = 300)

ggsave(file.path(out_dir, "Polyps_three_species_row.png"),
       final_row_p, width = 12, height = 4.5, dpi = 300)


# output group summaries (means, SD, n)

summary_tables <- list(
  growth = data %>%
    group_by(Species, Lab) %>%
    summarise(
      mean = mean(Growth_cm, na.rm = TRUE),
      sd   = sd(Growth_cm, na.rm = TRUE),
      n    = dplyr::n(),
      .groups = "drop"
    ) %>%
    mutate(metric = "Growth_cm"),

  polyps = data %>%
    group_by(Species, Lab) %>%
    summarise(
      mean = mean(Polyps_change, na.rm = TRUE),
      sd   = sd(Polyps_change, na.rm = TRUE),
      n    = dplyr::n(),
      .groups = "drop"
    ) %>%
    mutate(metric = "Polyps_change")
)

summary_out <- bind_rows(summary_tables$growth, summary_tables$polyps) %>%
  relocate(metric, Species, Lab, mean, sd, n)

write.csv(summary_out, file.path(out_dir, "group_summaries_growth_polyps.csv"), row.names = FALSE)

message("\n files saved in: ", normalizePath(out_dir))
message(" - Growth_three_species_row.png")
message(" - Polyps_three_species_row.png")
message(" - group_summaries_growth_polyps.csv\n")


#  export summary tables for statistical tests 

# Growth_cm

# re-run tests and save their summaries
growth_anova <- tryCatch({
  aov(Growth_cm ~ Species * Lab, data = data)
}, error = function(e) NULL)

growth_shapiro <- shapiro.test(residuals(growth_anova))
growth_levene  <- car::leveneTest(Growth_cm ~ Species * Lab, data = data)

# pairwise (nonparametric) summary
growth_pairwise <- data %>%
  group_by(Species) %>%
  rstatix::pairwise_wilcox_test(Growth_cm ~ Lab, p.adjust.method = "BH")

# save to .csv
write.csv(as.data.frame(growth_levene),
          file.path(out_dir, "growth_levene_summary.csv"),
          row.names = FALSE)
write.csv(broom::tidy(growth_shapiro),
          file.path(out_dir, "growth_shapiro_summary.csv"),
          row.names = FALSE)
write.csv(growth_pairwise,
          file.path(out_dir, "growth_pairwise_wilcoxon.csv"),
          row.names = FALSE)

# Polyps_change

# re-run tests and save their summaries

polyps_anova <- tryCatch({
  aov(Polyps_change ~ Species * Lab, data = data)
}, error = function(e) NULL)

polyps_shapiro <- shapiro.test(residuals(polyps_anova))
polyps_levene  <- car::leveneTest(Polyps_change ~ Species * Lab, data = data)

# pairwise (nonparametric) summary
polyps_pairwise <- data %>%
  group_by(Species) %>%
  rstatix::pairwise_wilcox_test(Polyps_change ~ Lab, p.adjust.method = "BH")

# Save to csvs
write.csv(as.data.frame(polyps_levene),
          file.path(out_dir, "polyps_levene_summary.csv"),
          row.names = FALSE)
write.csv(broom::tidy(polyps_shapiro),
          file.path(out_dir, "polyps_shapiro_summary.csv"),
          row.names = FALSE)
write.csv(polyps_pairwise,
          file.path(out_dir, "polyps_pairwise_wilcoxon.csv"),
          row.names = FALSE)

# kruskal–wallis summaries
growth_kw <- data %>%
  group_by(Species) %>%
  rstatix::kruskal_test(Growth_cm ~ Lab)
polyps_kw <- data %>%
  group_by(Species) %>%
  rstatix::kruskal_test(Polyps_change ~ Lab)

write.csv(growth_kw,  file.path(out_dir, "growth_kruskal_summary.csv"),  row.names = FALSE)
write.csv(polyps_kw, file.path(out_dir, "polyps_kruskal_summary.csv"), row.names = FALSE)

```

```{r}
lab_colors <- list(
  "Muricea pendula" = c("GAL" = "grey30", "GNV" = "grey60", "HML" = "grey80"),
  "Swiftia exserta" = c("GAL" = "#cc7a00", "GNV" = "#e69f00", "HML" = "#ffb84d"),
  "Thesea nivea"          = c("GAL" = "#6a51a3", "GNV" = "#9467BD", "HML" = "#c5b0d5")
)

make_species_plot <- function(sp) {
  sp_data <- data |> filter(Species == sp)
  p <- ggplot(sp_data, aes(x = Lab, y = Growth_cm, fill = Lab)) +
    geom_boxplot(color = "black", outlier.shape = 16, outlier.size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.6) + 
    scale_fill_manual(values = lab_colors[[sp]]) +
    theme_minimal(base_size = 13) +
    labs(title = bquote(italic(.(sp))), x = "Lab", y = "Growth Rate (cm/year)")
  p
}

p_muricea <- make_species_plot("Muricea pendula")
p_swiftia <- make_species_plot("Swiftia exserta")
p_thesea  <- make_species_plot("Thesea nivea")

p_muricea
p_swiftia
p_thesea

# save each:
ggsave("Muricea_pendula_boxplot.png", p_muricea, width = 6, height = 4, dpi = 300)
ggsave("Swiftia_exserta_boxplot.png", p_swiftia, width = 6, height = 4, dpi = 300)
ggsave("Thesea_boxplot.png", p_thesea, width = 6, height = 4, dpi = 300)
```



#running analyses to better understand significant differences between species

#install packages                       
```{r}
install.packages("tidyr")
install.packages("tibble")
```


```{r}
# read in libs
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(car)
library(tidyr)
library(tibble)
library(purrr)

# load data (same as above- just putting in to be sure)
data <- read_excel("CLGE data.xlsx", sheet = "data")
colnames(data) <- c("Lab", "Species", "Colony_ID", "Growth_cm", "Polyps_change")
data <- data |> mutate(Species = factor(Species), Lab = factor(Lab))

# assign species colors (same as above- just putting in to be sure)
species_colors <- c(
  "Muricea pendula" = "grey60",
  "Swiftia exserta" = "#E69F00", # orange
  "Thesea nivea"          = "#9467BD"  # purple
)

# create stat table with manual y-positions
make_sig_table <- function(df, value_col, x_col = "Species",
                           method = c("auto"), # we decide inside based on assumptions
                           step = 0.08) {
  # df: data.frame with columns x_col and value_col
  # value_col: "Growth_cm" or "Polyps_change"
  # step: fraction of data range to stack bars
  yvals <- df[[value_col]]
  xmax  <- max(yvals, na.rm = TRUE)
  xmin  <- min(yvals, na.rm = TRUE)
  yrng  <- xmax - xmin
  if (!is.finite(yrng) || yrng == 0) yrng <- 1

  # determine test path via assumptions
  form <- as.formula(paste(value_col, "~", x_col))
  aov_fit <- aov(form, data = df)
  sw <- tryCatch(shapiro.test(residuals(aov_fit)), error = function(e) list(p.value = 0))
  lv <- tryCatch(car::leveneTest(form, data = df), error = function(e) data.frame(`Pr(>F)` = 0))

  use_anova <- is.numeric(sw$p.value) && sw$p.value > 0.05 &&
               is.numeric(lv$`Pr(>F)`[1]) && lv$`Pr(>F)`[1] > 0.05

  if (use_anova) {
    # Tukey via rstatix (tidy with group1/group2 to easily understand results)
    pw <- rstatix::tukey_hsd(aov_fit) %>%
      transmute(group1, group2, p.adj,
                p.signif = case_when(
                  p.adj <= 0.001 ~ "***",
                  p.adj <= 0.01  ~ "**",
                  p.adj <= 0.05  ~ "*",
                  TRUE ~ "ns"
                ))
  } else {
    # Kruskal + pairwise Wilcoxon
    pw <- rstatix::pairwise_wilcox_test(df, formula = form, p.adjust.method = "BH") %>%
      transmute(group1, group2, p.adj,
                p.signif = case_when(
                  p.adj <= 0.001 ~ "***",
                  p.adj <= 0.01  ~ "**",
                  p.adj <= 0.05  ~ "*",
                  TRUE ~ "ns"
                ))
  }

  # keep only significant rows for plotting (hide.ns = TRUE)
  pw <- pw %>% filter(p.signif != "ns")

  if (nrow(pw) == 0) {
    # return an empty frame with required cols so stat_pvalue_manual is happy
    return(tibble(group1 = character(), group2 = character(),
                  p.adj = numeric(), p.signif = character(),
                  y.position = numeric()))
  }

  # assign stacked y positions (just above the current max, then step up)
  # ff many comparisons, this will stack them neatly.
  base_y <- xmax + 0.05 * yrng
  pw <- pw %>%
    mutate(y.position = base_y + (row_number() - 1) * (step * yrng))

  return(pw)
}

# looking at growth between species and plotting
growth_sig <- make_sig_table(data, value_col = "Growth_cm")

p_growth_species <- ggplot(data, aes(x = Species, y = Growth_cm, fill = Species)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.6) +
  geom_boxplot(color = "black", outlier.shape = 16, outlier.size = 1.5) +
  scale_fill_manual(values = species_colors) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(face = "italic")) +
  labs(title = "Annual average growth rate/branch (cm)",
       x = "Species", y = "Growth (cm/branch/year)")

# Add sig bars if any
if (nrow(growth_sig) > 0) {
  p_growth_species <- p_growth_species +
    ggpubr::stat_pvalue_manual(
      growth_sig, label = "p.signif", tip.length = 0.01, hide.ns = TRUE
    ) +
    # ensure plot range accommodates bars
    expand_limits(y = max(growth_sig$y.position, na.rm = TRUE) * 1.03)
}
p_growth_species


# looking at change in polyps between species
polyps_sig <- make_sig_table(data, value_col = "Polyps_change")

p_polyps_species <- ggplot(data, aes(x = Species, y = Polyps_change, fill = Species)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "red", linewidth = 0.6) +
  geom_boxplot(color = "black", outlier.shape = 16, outlier.size = 1.5) +
  scale_fill_manual(values = species_colors) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(face = "italic")) +
  labs(title = "%\u0394 Polyps over 6 mo.",
       x = "Species", y = "%\u0394 Polyps over 6 mo.")

if (nrow(polyps_sig) > 0) {
  p_polyps_species <- p_polyps_species +
    ggpubr::stat_pvalue_manual(
      polyps_sig, label = "p.signif", tip.length = 0.01, hide.ns = TRUE
    ) +
    expand_limits(y = max(polyps_sig$y.position, na.rm = TRUE) * 1.03)
}
p_polyps_species

#  save as 300 dpi
ggsave("species_growth_boxplot_vproof.png", p_growth_species, width = 6.5, height = 5, dpi = 300)
ggsave("species_polyps_boxplot_vproof.png", p_polyps_species, width = 6.5, height = 5, dpi = 300)


# export and show summary results of statistical tests

# output to working directory
out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)

# growth data
cat("\n================= GROWTH BETWEEN SPECIES =================\n")

# test assumptions for ANOVA or Kruskal decision
growth_model <- aov(Growth_cm ~ Species, data = data)
shapiro_g <- shapiro.test(residuals(growth_model))
levene_g  <- car::leveneTest(Growth_cm ~ Species, data = data)

# print assumption results
cat("Shapiro-Wilk test for normality:\n"); print(shapiro_g)
cat("\nLevene’s test for homogeneity of variances:\n"); print(levene_g)

if (shapiro_g$p.value > 0.05 & levene_g$`Pr(>F)`[1] > 0.05) {
  cat("\n Running one-way ANOVA (Growth_cm ~ Species)\n")
  growth_anova <- summary(growth_model)
  print(growth_anova)

  cat("\nPost-hoc: Tukey HSD\n")
  growth_tukey <- TukeyHSD(growth_model)
  print(growth_tukey)

  # export
  broom::tidy(growth_model) %>%
    write.csv(file.path(out_dir, "growth_anova_results.csv"), row.names = FALSE)
  broom::tidy(growth_tukey) %>%
    write.csv(file.path(out_dir, "growth_tukey_posthoc.csv"), row.names = FALSE)

} else {
  cat("\n Running non-parametric Kruskal–Wallis + pairwise Wilcoxon\n") #data was nonparametric
  growth_kw <- rstatix::kruskal_test(data, Growth_cm ~ Species)
  print(growth_kw)

  growth_pairwise <- data %>%
    rstatix::pairwise_wilcox_test(Growth_cm ~ Species, p.adjust.method = "BH")
  print(growth_pairwise) #siginifcant differences observed between species

  # export
  write.csv(growth_kw, file.path(out_dir, "growth_kruskal_results.csv"), row.names = FALSE)
  write.csv(growth_pairwise, file.path(out_dir, "growth_pairwise_wilcoxon.csv"), row.names = FALSE)
}


# ----------------- POLYPS -----------------
cat("\n================= POLYPS BETWEEN SPECIES =================\n")

polyps_model <- aov(Polyps_change ~ Species, data = data)
shapiro_p <- shapiro.test(residuals(polyps_model))
levene_p  <- car::leveneTest(Polyps_change ~ Species, data = data)

# print assumption results
cat("Shapiro-Wilk test for normality:\n"); print(shapiro_p)
cat("\nLevene’s test for homogeneity of variances:\n"); print(levene_p)

if (shapiro_p$p.value > 0.05 & levene_p$`Pr(>F)`[1] > 0.05) {
  cat("\n Running one-way ANOVA (Polyps_change ~ Species)\n") #data met assumptions
  polyps_anova <- summary(polyps_model)
  print(polyps_anova) #siginifcant differences observed between species

  cat("\nPost-hoc: Tukey HSD\n")
  polyps_tukey <- TukeyHSD(polyps_model)
  print(polyps_tukey)

  broom::tidy(polyps_model) %>%
    write.csv(file.path(out_dir, "polyps_anova_results.csv"), row.names = FALSE)
  broom::tidy(polyps_tukey) %>%
    write.csv(file.path(out_dir, "polyps_tukey_posthoc.csv"), row.names = FALSE)

} else {
  cat("\n Running non-parametric Kruskal–Wallis + pairwise Wilcoxon\n")
  polyps_kw <- rstatix::kruskal_test(data, Polyps_change ~ Species)
  print(polyps_kw)

  polyps_pairwise <- data %>%
    rstatix::pairwise_wilcox_test(Polyps_change ~ Species, p.adjust.method = "BH")
  print(polyps_pairwise)

  #export 
  
  write.csv(polyps_kw, file.path(out_dir, "polyps_kruskal_results.csv"), row.names = FALSE)
  write.csv(polyps_pairwise, file.path(out_dir, "polyps_pairwise_wilcoxon.csv"), row.names = FALSE)
}

cat("\nStatistical results saved in:", normalizePath(out_dir), "\n")
```
                 
```{r}
# printing summary tables with mean, SD, etc.

library(dplyr)

# output directory
out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)

# function for SE
se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# compute summaries
species_summary <- data %>%
  group_by(Species) %>%
  summarise(
    n_growth  = sum(!is.na(Growth_cm)),
    mean_growth = mean(Growth_cm, na.rm = TRUE),
    se_growth   = se(Growth_cm),
    n_polyps  = sum(!is.na(Polyps_change)),
    mean_polyps = mean(Polyps_change, na.rm = TRUE),
    se_polyps   = se(Polyps_change)
  ) %>%
  mutate(
    Growth_mean_SE  = sprintf("%.3f ± %.3f", mean_growth, se_growth),
    Polyps_mean_SE  = sprintf("%.2f ± %.2f", mean_polyps, se_polyps)
  ) %>%
  select(Species, n_growth, Growth_mean_SE, n_polyps, Polyps_mean_SE)

# print in console
cat("\n================= SUMMARY TABLE: Mean ± SE per Species =================\n")
print(species_summary)

# Export to csv
write.csv(species_summary,
          file.path(out_dir, "species_mean_SE_summary.csv"),
          row.names = FALSE)

cat("\n Summary table saved to:", normalizePath(file.path(out_dir, "species_mean_SE_summary.csv")), "\n")
```
