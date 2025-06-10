# Load required libraries
library(limma)
library(ggplot2)
library(reshape2)
library(readxl)

expr <- read_excel("C:/Users/crtuser/OneDrive - TCDUD.onmicrosoft.com/Documents/PhD/Project/Data/Gearoid_Proteomics/UTS vs TS.xlsx", 
                   sheet = "Sheet3")

# Convert to data frame and set PG.Genes as rownames
expr <- as.data.frame(expr)
expr <- expr[!is.na(expr$PG.Genes), ]
expr$PG.Genes <- make.unique(as.character(expr$PG.Genes))
rownames(expr) <- expr$PG.Genes
expr <- expr[, -1]  # remove PG.Genes column

# Log2-transform if not already done
expr <- log2(expr + 1)
samples <- colnames(expr)

# Extract group (TS vs T) from sample name
group <- ifelse(grepl("^TS", samples), "TS", "T")

# Extract pair ID (e.g., "1", "2", etc.)
pair_id <- sub("^[A-Z]+", "", samples)

# Create metadata
meta <- data.frame(Sample = samples, Group = group, PairID = pair_id)

# Reorder columns of expr to match metadata
expr <- expr[, meta$Sample]

# --- Run paired limma ---

# Create design matrix and factors
group <- factor(meta$Group, levels = c("TS", "T"))
pair <- factor(meta$PairID)
design <- model.matrix(~group)

# Estimate correlation between pairs
dupcor <- duplicateCorrelation(expr, design, block = pair)

# Fit linear model with blocking
fit <- lmFit(expr, design, block = pair, correlation = dupcor$consensus)
fit <- eBayes(fit)

# Extract differential expression results
results <- topTable(fit, coef = "groupT", number = Inf, adjust.method = "fdr")

# View results
View(results)

# Save to file
#write.csv(results, "limma_paired_results_UTS_vs_TS.csv", row.names = TRUE)

# --- Optional: Plot one protein of interest ---
protein <- "LIPA"  # replace with a protein in your data

expr_melt <- data.frame(
  Sample = colnames(expr),
  Expression = as.numeric(unlist(expr[protein, ])),
  Group = meta$Group,
  PairID = meta$PairID
)

ggplot(expr_melt, aes(x = Group, y = Expression, color = Group)) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  geom_jitter(width = 0.1, size = 2) +
  labs(title = paste("Expression of", protein), y = "Log2 Expression") +
  theme_minimal()


##### WHAT IF A PCA DROPPED ON YOUR HEAD
expr_pca <- expr[apply(expr, 1, function(x) sd(x) > 0), ]

# Transpose for PCA (samples as rows, proteins as columns)
expr_t <- t(expr_pca)

# Run PCA
pca <- prcomp(expr_t, scale. = TRUE)

# Create data frame for plotting
pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = meta$Group
)

# Plot PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  labs(title = "UTS vs TS", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),        # Remove gridlines
    axis.line = element_line(),          # Add axis lines
    plot.title = element_text(hjust = 0.5)
  )

# Loadings (rotation matrix): rows = proteins, columns = PCs
loadings <- pca$rotation
head(loadings[, 1:5])  # View first 5 PCs

get_top_loadings <- function(pc = "PC1", n = 10) {
  if (!(pc %in% colnames(loadings))) stop("That PC doesn't exist.")
  sort(abs(loadings[, pc]), decreasing = TRUE)[1:n]
}

# Example: Top 10 for PC2
get_top_loadings("PC2", 20)


