# Load required libraries
library(limma)
library(ggplot2)
library(reshape2)
library(plotly)
library(rgl) # Often used for interactive 3D plots, not directly used in the 2D plots below but good to have if needed.
library(readxl) # For reading Excel files

# --- Data Loading and Preprocessing ---
# Load data from the specified Excel file
# IMPORTANT: Ensure the sheet name in your "UT vs UTS.xlsx" file is correct,
# currently assuming "Low exp removed" based on your previous script.
expr <- read_excel("C:/Users/crtuser/OneDrive - TCDUD.onmicrosoft.com/Documents/PhD/Project/Data/Gearoid_Proteomics/UT vs UTS.xlsx",
                   sheet = "Low exp excluded")

# Convert to data frame and set PG.Genes as rownames
expr <- as.data.frame(expr)
expr <- expr[!is.na(expr$PG.Genes), ] # Remove rows where PG.Genes is NA
expr$PG.Genes <- make.unique(as.character(expr$PG.Genes)) # Make gene names unique
rownames(expr) <- expr$PG.Genes
expr <- expr[, -1]  # remove PG.Genes column as it's now rownames

# Log2-transform if not already done. Add 1 to avoid issues with zero values.
expr <- log2(expr + 1)

# Get sample names from column headers
samples <- colnames(expr)

# Extract group (UTS vs UT) from sample name
# Assumes sample names start with "UTS" for UTS group and "UT" for UT group
group <- ifelse(grepl("^UTS", samples), "UTS", "UT")

# Extract pair ID (e.g., "1", "2", etc.)
pair_id <- sub("^[A-Z]+", "", samples) # Removes leading letters to get the number

# Create metadata
meta <- data.frame(Sample = samples, Group = group, PairID = pair_id)

# Reorder columns of expr to match metadata (ensures consistency)
expr <- expr[, meta$Sample]

# --- Run paired limma analysis ---

# Create design matrix and factors
# Set "UT" as the reference level, so "UTS" will be the comparison group.
group <- factor(meta$Group, levels = c("UT", "UTS"))
pair <- factor(meta$PairID)
design <- model.matrix(~group)

# Estimate correlation between pairs. This is crucial for paired analysis.
dupcor <- duplicateCorrelation(expr, design, block = pair)

# Fit linear model with blocking for paired design
fit <- lmFit(expr, design, block = pair, correlation = dupcor$consensus)
fit <- eBayes(fit) # Apply empirical Bayes smoothing to the standard errors

# Extract differential expression results for the "groupUTS" contrast
results <- topTable(fit, coef = "groupUTS", number = Inf, adjust.method = "fdr")

# View the top results
View(results)

results_ordered <- results[order(results$t, decreasing = TRUE), ]

# View the head (most upregulated in UTS) and tail (most downregulated in UTS) of the ordered results
message("\nTop upregulated proteins (UTS vs UT):")
head(results_ordered)
message("\nTop downregulated proteins (UTS vs UT):")
tail(results_ordered)

# Save results to a CSV file
write.csv(results_ordered, "C:/Users/crtuser/OneDrive - TCDUD.onmicrosoft.com/Documents/PhD/Project/Data/Gearoid_Proteomics/UT_vs_UTS_paired_limma_ordered.csv", row.names = TRUE)

# --- Optional: Plot expression of a single protein ---
# Replace "FOS" with a protein identifier present in your data's PG.Genes column
protein_of_interest <- "FOS"

# Check if the protein exists in the data
if (protein_of_interest %in% rownames(expr)) {
  expr_melt <- data.frame(
    Sample = colnames(expr),
    Expression = as.numeric(unlist(expr[protein_of_interest, ])),
    Group = meta$Group,
    PairID = meta$PairID
  )
  
  ggplot(expr_melt, aes(x = Group, y = Expression, color = Group)) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    geom_jitter(width = 0.1, size = 2) +
    labs(title = paste("Expression of", protein_of_interest), y = "Log2 Expression") +
    theme_minimal() +
    # Add lines connecting paired samples
    geom_line(aes(group = PairID), color = "gray", linetype = "dashed")
} else {
  message(paste("Protein '", protein_of_interest, "' not found in the data. Skipping protein plot."))
}

expr_pca <- expr[apply(expr, 1, function(x) sd(x) > 0), ]

# Transpose the expression matrix for PCA (samples as rows, proteins as columns)
expr_t <- t(expr_pca)

# Run PCA. `scale. = TRUE` standardizes the variables (proteins) to have unit variance.
pca <- prcomp(expr_t, scale. = TRUE)

# Create data frame for PCA plotting
pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1], # Principal Component 1
  PC2 = pca$x[, 2], # Principal Component 2 (changed from PC3 in your original script for standard 2D PCA)
  Group = meta$Group
)

# Plot PCA with PC1 and PC2
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  labs(title = "PCA of UT vs UTS Samples", x = paste0("PC1 (", round(summary(pca)$importance[2, 1]*100, 2), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2, 2]*100, 2), "%)")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), # Remove gridlines
    axis.line = element_line(),   # Add axis lines
    plot.title = element_text(hjust = 0.5) # Center the title
  )