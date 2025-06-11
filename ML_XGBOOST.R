library(readxl)
library(dplyr)
library(xgboost)
library(caret)
library(Matrix)

# --- Load Data ---
base_path <- "C:/Users/crtuser/OneDrive - TCDUD.onmicrosoft.com/Documents/PhD/Project/Data/Gearoid_Proteomics/"
file1 <- paste0(base_path, "T vs TS.xlsx")
file2 <- paste0(base_path, "UT vs UTS.xlsx")

data_t_ts <- read_excel(file1, sheet = "Low exp removed")
data_ut_uts <- read_excel(file2, sheet = "Low exp excluded")

# --- Clean and Normalize ---
clean_data <- function(df) {
  df %>%
    filter(!is.na(PG.Genes)) %>%
    group_by(PG.Genes) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
}

data1 <- clean_data(data_t_ts)
data2 <- clean_data(data_ut_uts)

combined_expr <- full_join(data1, data2, by = "PG.Genes") %>%
  filter(!is.na(PG.Genes))

expr_matrix <- as.data.frame(combined_expr)
rownames(expr_matrix) <- expr_matrix$PG.Genes
expr_matrix <- expr_matrix[, -1]
expr_matrix <- log2(expr_matrix + 1)
expr_t <- t(expr_matrix)

# --- Create Metadata ---
sample_ids <- rownames(expr_t)
stim_status <- ifelse(grepl("TS|UTS", sample_ids), "Stim", "Unstim")
treat_status <- ifelse(grepl("^T|TS", sample_ids), "Treated", "Untreated")

meta_df <- data.frame(Sample = sample_ids,
                      Stim = factor(stim_status),
                      Treat = factor(treat_status),
                      stringsAsFactors = FALSE)

expr_df <- data.frame(meta_df, expr_t)

# --- Select Top Genes ---
top_n <- 50
gene_vars <- apply(expr_df[, -(1:3)], 2, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:top_n]
expr_df_filtered <- expr_df[, c("Stim", "Treat", top_genes)]

# --- Train with Cross-Validation ---
train_xgb_cv <- function(df, label_col) {
  df[[label_col]] <- factor(df[[label_col]])
  
  fit_control <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"
  )
  
  xgb_grid <- expand.grid(
    nrounds = 50,
    max_depth = 3,
    eta = 0.3,
    gamma = 0,
    colsample_bytree = 0.8,
    min_child_weight = 1,
    subsample = 1
  )
  
  model <- train(
    form = as.formula(paste(label_col, "~ .")),
    data = df,
    method = "xgbTree",
    trControl = fit_control,
    tuneGrid = xgb_grid,
    metric = "ROC"
  )
  
  return(model)
}

# --- Run Models ---
cat("===== Stimulated vs Unstimulated (CV) =====\n")
stim_model <- train_xgb_cv(expr_df_filtered[, c("Stim", top_genes)], "Stim")
print(stim_model)
print(confusionMatrix(stim_model$pred$pred, stim_model$pred$obs))

cat("===== Treated vs Untreated (CV) =====\n")
treat_model <- train_xgb_cv(expr_df_filtered[, c("Treat", top_genes)], "Treat")
print(treat_model)
print(confusionMatrix(treat_model$pred$pred, treat_model$pred$obs))


##### PLOTTING
importance <- varImp(treat_model, scale = TRUE)
print(importance)
plot(importance, top = 20)

importance <- varImp(stim_model, scale = TRUE)
print(importance)
plot(importance, top = 20)


##### MULTICLASS PREDICTION:
# --- Load Libraries ---
library(readxl)
library(dplyr)
library(xgboost)
library(caret)
library(Matrix)

# --- File Paths ---
base_path <- "C:/Users/crtuser/OneDrive - TCDUD.onmicrosoft.com/Documents/PhD/Project/Data/Gearoid_Proteomics/"
file1 <- paste0(base_path, "T vs TS.xlsx")
file2 <- paste0(base_path, "UT vs UTS.xlsx")

# --- Load Excel Sheets ---
data_t_ts <- read_excel(file1, sheet = "Low exp removed")
data_ut_uts <- read_excel(file2, sheet = "Low exp excluded")

# --- Clean and Average Duplicate Genes ---
data_t_ts_clean <- data_t_ts %>%
  filter(!is.na(PG.Genes)) %>%
  group_by(PG.Genes) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

data_ut_uts_clean <- data_ut_uts %>%
  filter(!is.na(PG.Genes)) %>%
  group_by(PG.Genes) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

# --- Combine and Normalize ---
combined_expr <- full_join(data_t_ts_clean, data_ut_uts_clean, by = "PG.Genes") %>%
  filter(!is.na(PG.Genes))

expr_matrix <- as.data.frame(combined_expr)
rownames(expr_matrix) <- expr_matrix$PG.Genes
expr_matrix <- expr_matrix[, -1]
expr_matrix <- log2(expr_matrix + 1)
expr_t <- t(expr_matrix)

# --- Metadata and Multiclass Labels ---
sample_ids <- rownames(expr_t)
class_labels <- gsub("[0-9]", "", sample_ids)  # Get UT, UTS, T, TS
meta_df <- data.frame(Sample = sample_ids, Class = factor(class_labels))

# --- Join metadata and expression ---
expr_df <- data.frame(meta_df, expr_t)

# --- Select Top N Variable Genes ---
top_n <- 50
gene_vars <- apply(expr_df[, -(1:2)], 2, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:top_n]
expr_df_filtered <- expr_df[, c("Class", top_genes)]

# --- Encode Labels as 0-based integers ---
label_numeric <- as.numeric(expr_df_filtered$Class) - 1
dmat <- xgb.DMatrix(data = as.matrix(expr_df_filtered[, -1]), label = label_numeric)

# --- Cross-Validation and Training ---
set.seed(42)
cv_folds <- createFolds(label_numeric, k = 5, list = TRUE, returnTrain = FALSE)

params <- list(
  objective = "multi:softprob",
  num_class = length(unique(label_numeric)),
  eval_metric = "mlogloss"
)

xgb_cv <- xgb.cv(
  params = params,
  data = dmat,
  folds = cv_folds,
  nrounds = 100,
  early_stopping_rounds = 10,
  verbose = 1,
  prediction = TRUE
)

# --- Final model ---
best_nrounds <- xgb_cv$best_iteration
xgb_model <- xgb.train(
  params = params,
  data = dmat,
  nrounds = best_nrounds,
  verbose = 1
)

# --- Predict and Evaluate ---
pred_probs <- predict(xgb_model, as.matrix(expr_df_filtered[, -1]))
pred_matrix <- matrix(pred_probs, ncol = length(unique(label_numeric)), byrow = TRUE)
pred_labels <- max.col(pred_matrix) - 1  # Predicted class

confusion <- confusionMatrix(
  factor(pred_labels, labels = levels(expr_df_filtered$Class)),
  factor(label_numeric, labels = levels(expr_df_filtered$Class))
)

print(confusion)

# --- Top Important Features ---
importance <- xgb.importance(model = xgb_model)
cat("Top genes overall:\n")
print(head(importance, 20))


##### PLOTTING: 
library(ggplot2)

# Example: Plot expression of one gene
gene_to_plot <- 'PRKAR1B'  # Top gene

ggplot(expr_df, aes_string(x = "Class", y = gene_to_plot, fill = "Class")) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Hide default outliers for clarity
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, aes_string(color = "Class")) +  # Add individual points
  labs(title = paste("Expression of", gene_to_plot), y = "Log2 Normalized Expression") +
  theme_minimal()

