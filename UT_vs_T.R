library(limma)
library(ggplot2)
library(reshape2)

# Load expression data
expr <- read.csv("C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Gearoid_Proteomics\\proteomics_for_limma.csv")
expr$rownames <- make.unique(as.character(expr$rownames))
rownames(expr) <- expr$rownames
expr <- expr[, -1]

# Load metadata
meta <- read.csv("C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Gearoid_Proteomics\\sample_info.csv")

# Match column order
expr <- expr[, meta$Sample]

# Add a PairID column based on matching T/UT samples
# Assumes UT1 <-> T1, UT2 <-> T2, etc.
meta$PairID <- rep(1:(nrow(meta)/2), 2)

# Create design matrix
group <- factor(meta$Group, levels = c("UT", "T"))
pair <- factor(meta$PairID)
design <- model.matrix(~group)

# Fit model with duplicateCorrelation to account for pairing
dupcor <- duplicateCorrelation(expr, design, block = pair)
fit <- lmFit(expr, design, block = pair, correlation = dupcor$consensus)
fit <- eBayes(fit)

# Get results (coefficient for T vs UT)
results <- topTable(fit, coef = "groupT", number = Inf, adjust.method = "fdr")

# Plot one protein of interest
protein <- "CUTA"
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

View(results)

write.csv(results, 'C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Gearoid_Proteomics\\UT_vs_T_paired_limma.csv')

print(
  ggplot(expr_melt, aes(x = Group, y = Expression, color = Group)) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    geom_jitter(width = 0.1, size = 2) +
    labs(title = paste("Expression of", protein), y = "Log2 Expression") +
    theme_minimal()
)

