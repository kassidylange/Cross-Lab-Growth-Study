```{r}
#install packages
install.packages("tidyverse")
install.packages("FactoMineR")
install.packages("factoextra")
```

```{r}
# load packages
library(tidyverse)
library(FactoMineR)
library(factoextra)

# read df
df <- read.csv("CLGE tanks.csv")

# make lab factor
df$Lab <- factor(df$Lab, levels = c("HML", "GNV", "GAL"))

# assign parameters
params <- df %>%
  select(Temp, Sal, pH, NH4, NO2, NO3, PO4, Alk, Ca, Mg)

# run a PCA
pca_res <- PCA(params, scale.unit = TRUE, graph = FALSE)

# assign custom colors for labs
lab_colors <- c(
  "HML" = "purple",
  "GNV" = "blue",
  "GAL" = "green"
)

# generate contribution table
contrib_table <- pca_res$var$contrib
write.csv(contrib_table, "PCA_variable_contributions.csv", row.names = TRUE)

# generate eigenvector table
eigen_table <- pca_res$eig
write.csv(eigen_table, "PCA_eigenvalues.csv", row.names = TRUE)

# make PCA bipolt with labelled eigenvectors
plot_labeled <- fviz_pca_biplot(
  pca_res,
  geom.ind = "point",
  col.ind = df$Lab,
  palette = lab_colors,
  label = "var",         # show labels
  col.var = "black",
  addEllipses = TRUE,
  ellipse.type = "confidence",
  ellipse.fill = TRUE,           
  ellipse.alpha = 0.2,             
  ellipse.level = 0.95,           # 95% confidence ellipse
  repel = TRUE,
  title = "PCA Biplot with Eigenvector Labels"
)

ggsave("PCA_biplot_labeled_300dpi.png", plot_labeled,
       dpi = 300, width = 8, height = 6)

# make PCA bipolt without labelled eigenvectors
plot_unlabeled <- fviz_pca_biplot(
  pca_res,
  geom.ind = "point",
  col.ind = df$Lab,
  palette = lab_colors,
  label = "none",        # remove labels
  col.var = "black",
  addEllipses = TRUE,
    ellipse.type = "confidence",
  ellipse.fill = TRUE,           
  ellipse.alpha = 0.2,              
  ellipse.level = 0.95,           # 95% confidence
  repel = FALSE,
  title = "PCA Biplot without Eigenvector Labels"
)

ggsave("PCA_biplot_unlabeled_300dpi.png", plot_unlabeled,
       dpi = 300, width = 8, height = 6)

# scree plots to understand PCs

# plotting
scree_plot <- fviz_eig(
  pca_res,
  addlabels = TRUE,
  barfill = "gray70",
  barcolor = "black",
  linecolor = "black",
  ggtheme = theme_minimal(),
  main = "Scree Plot of PCA Eigenvalues"
)

ggsave("PCA_scree_plot_300dpi.png",
       scree_plot, dpi = 300, width = 8, height = 6)

```

```{r}
#install packages
install.packages("MVN")
install.packages("car")
install.packages("vegan")
install.packages("mvnTest")
install.packages("heplots")
```
```{r}
# load
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(vegan)

# read df
df <- read.csv("CLGE tanks.csv")
df$Lab <- factor(df$Lab, levels = c("HML", "GNV", "GAL"))

# select parameters
params <- df %>% select(Temp, Sal, pH, NH4, NO2, NO3, PO4, Alk, Ca, Mg)

# pca
pca_res <- PCA(params, scale.unit = TRUE, graph = FALSE)

# extract first four PCs (=71% of variance)
pca_scores <- as.data.frame(pca_res$ind$coord)[, 1:4]
pca_scores$Lab <- df$Lab

# code to fix missing values
if(any(is.na(pca_scores))){
  na_rows <- rowSums(is.na(pca_scores)) == 0
  pca_scores_clean <- pca_scores[na_rows, 1:4]
  labs_clean <- pca_scores$Lab[na_rows]
} else {
  pca_scores_clean <- pca_scores[, 1:4]
  labs_clean <- pca_scores$Lab
}

# run a PERMANOVA
permanova_res <- adonis2(pca_scores_clean ~ labs_clean, method = "euclidean", permutations = 999)
print(permanova_res)

# pairwise
labs <- levels(df$Lab)
pairs <- combn(labs, 2, simplify = FALSE)

pairwise_results <- lapply(pairs, function(p){
  
  # pull out subsets
  subset_rows <- df$Lab %in% p
  pca_subset <- pca_scores_clean[subset_rows, ]
  labs_subset <- df$Lab[subset_rows]
  
  # Run PERMANOVA
  res <- adonis2(pca_subset ~ labs_subset, method = "euclidean", permutations = 999)
  
  # Extract R2 and p-value
  data.frame(
    Lab1 = p[1],
    Lab2 = p[2],
    R2 = res$R2[1],
    F = res$F[1],
    p = res$`Pr(>F)`[1]
  )
})

pairwise_results_df <- bind_rows(pairwise_results) %>%
  mutate(p_adj = p.adjust(p, method = "BH"))  # adjust p-values for multiple comparisons

print(pairwise_results_df)

# bar plot of R2
plot_pairwise <- ggplot(pairwise_results_df, aes(x = comparison, y = R2, fill = R2)) +
  geom_col() +
  geom_text(aes(label = sig_label), vjust = -0.5, size = 5) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(x = "Lab Comparison", y = expression(R^2), title = "Pairwise PERMANOVA R² with Significance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# save
ggsave("Pairwise_PERMANOVA_R2_significance_300dpi.png", plot_pairwise, dpi = 300, width = 8, height = 6)

```
