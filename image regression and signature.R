# loading and formatting the data -----------------------------------------
KNN_data <- read_csv("/data/home/ha231442/HF/KNN_data.csv")
covar_all <- read_csv("/data/home/ha231442/HF/ckd_cmr_covariates_NA.csv")
image_dates <- read_csv("/data/home/ha231442/HF/image_dates_HF.csv")


KNN_data <-  KNN_data %>%
  as.data.frame() %>%
  mutate(eid = as.character(eid)) %>% 
  select(- ...1)

covar_all <-  covar_all %>%
  as.data.frame() %>%
  rename(eid = f.eid) %>%
  mutate(eid = as.character(eid)) %>% 
  select(- ...1)

image_dates <- image_dates %>%
  as.data.frame() %>%
  rename(eid = feid) %>%
  select(eid,Baseline_visit_date,Imaging_visit_date1)


# inverse-rank normalising the data  --------------------------------------
# Assuming the first column is 'eid' and should not be transformed
eid_column <- KNN_data[, 1, drop = FALSE]  # Drop = FALSE to keep it as a dataframe
data_to_normalize <- KNN_data[, -1]  # All columns except the first one

# Function to inverse-rank normalize a single column
inverse_rank_normalize <- function(column) {
  ranked <- rank(column, ties.method = "average", na.last = "keep")
  n <- length(na.omit(ranked))
  normalized <- qnorm((ranked - 0.5) / n)
  return(normalized)
}

# Apply the function to each column except the first one
normalized_data <- as.data.frame(lapply(data_to_normalize, inverse_rank_normalize))

# Combine the 'eid' column back with the normalized data
KNN_data_normalized <- cbind(eid_column, normalized_data)


####merge proteins and idp
prot_idp_covar <- merge(covar_all , KNN_data_normalized, by = "eid") ##4048 people in total 
prot_idp_covar <- merge(prot_idp_covar, image_dates, by = "eid")
prot_idp_covar$Baseline_visit_date <- as.Date(prot_idp_covar$Baseline_visit_date, format = "%d/%m/%Y")
prot_idp_covar$Imaging_visit_date1 <- as.Date(prot_idp_covar$Imaging_visit_date1, format = "%d/%m/%Y")
prot_idp_covar$time_dif <- interval(prot_idp_covar$Baseline_visit_date, prot_idp_covar$Imaging_visit_date1) / months(1)

write_csv(prot_idp_covar, "prot_idp_covar_simplified.csv")


##selecting variables for imaging regresison 
which(colnames(prot_idp_covar) == "LVEDV") ##from col 55 
which(colnames(prot_idp_covar) == "RVEF") ##to 66 
idp <- prot_idp_covar[, 55:66]
which(colnames(prot_idp_covar) == "a1bg")
prot <- prot_idp_covar[, 180:3090]




#  regression analysis ----------------------------------------------------
# Create a list of prot names and imaging names
prot_list <- colnames(prot) 
img_list <- colnames(idp)

# Create a blank result dataframe
lm_results <- data.frame()

# Model 1
for (i in img_list) {
  print(i)
  
  for (f in prot_list) {
    print(f)
    
    lm_model <- lm(prot_idp_covar[, i] ~ prot_idp_covar[, f] + Age + Sex + Height + Weight + egfr_creat_epi_2021_baseline + time_dif, data = prot_idp_covar)
    
    lm_results[f, paste0("beta.", i)] <- summary(lm_model)$coefficients[2, 1]  # beta value
    lm_results[f, paste0("p.", i)] <- summary(lm_model)$coefficients[2, 4]    # p value
  }
  
  # FDR adjust p value
  lm_results[, paste0("FDR.", i)] <- p.adjust(lm_results[, paste0("p.", i)], method = "BH")
  
  # Bonferroni adjust p value
  lm_results[, paste0("Bonferroni.", i)] <- p.adjust(lm_results[, paste0("p.", i)], method = "bonferroni")
}

write.csv(lm_results, "basic_lm_idp_simplified.csv")

##creating a new df for prots with bonferroni <0.05 at least for one idp
significant_results <- lm_results[apply(lm_results[, grep("Bonferroni", colnames(lm_results))], 1, function(x) any(x < 0.05)), ]
write.csv(significant_results, "significant_lm_idp.csv")



# penalisation ------------------------------------------------------------
library(glmnet)
library(Metrics)

sig_prots <-  significant_results$protnames

##removing NAs
##let's see where the NA in covars are 

subset_cols <- c("Age", "Sex", "Height", "Weight", "egfr_creat_epi_2021_baseline", "time_dif")

# Get columns with missing values and their counts
na_summary <- colSums(is.na(prot_idp_covar[, subset_cols]))
na_summary[na_summary > 0] ##so the only NAs are in eGFR col which has 211 NA. This is a continuous variable 

##mean imputation for eGFR
prot_idp_covar$egfr_creat_epi_2021_baseline <- ifelse(is.na(prot_idp_covar$egfr_creat_epi_2021_baseline), 
                                                      mean(prot_idp_covar$egfr_creat_epi_2021_baseline, na.rm = TRUE), 
                                                      prot_idp_covar$egfr_creat_epi_2021_baseline)


##removing IDPs with Nas 
na_rows <- which(rowSums(is.na(prot_idp_covar[, 55:66])) > 0)
prot_idp_covar <- prot_idp_covar[-na_rows, ]


##formatting for EN
prot_idp_covar$Sex[prot_idp_covar$Sex == "Female"] <- "2" #assigning these levels would have generated NAs if f and m were not character vectors 
prot_idp_covar$Sex[prot_idp_covar$Sex == "Male"] <- "1"
prot_idp_covar$Sex <- as.numeric(prot_idp_covar$Sex)

##keeping only significant proteins 
prot_cols <- colnames(prot_idp_covar[, 180:3090]) ##all the prots 
covar_cols <- prot_idp_covar[, c(1:179, 3093)]
# Extract protein names from the row names of significant_results
protein_names <- rownames(significant_results)
subset_prot_idp_covar <- prot_idp_covar[, intersect(prot_cols, protein_names)] ##only "significat" prots
subset_prot_idp_covar <- cbind(covar_cols, subset_prot_idp_covar)

#split into test and train s
set.seed(123)
index <- sample.int(nrow(prot_idp_covar), 0.7 * nrow(subset_prot_idp_covar))
train_data <- subset_prot_idp_covar[index, ]
test_data <- subset_prot_idp_covar[-index, ]


response_list <- colnames(idp)
protein_names <- rownames(significant_results)

# Initialize empty data frames to store coefficients and metrics
coef_results <- data.frame(protein = character(), coefficient = numeric(), phenotype = character(), stringsAsFactors = FALSE)
metrics_results <- data.frame(phenotype = character(), MSE = numeric(), MAE = numeric(), R_squared = numeric(), stringsAsFactors = FALSE)

set.seed(123)
# Loop over each image-derived phenotype in response_list
for (i in response_list) {
  
  # Extract the response variable (image-derived phenotype)
  response <- train_data[, i]
  
  # Prepare the predictor matrix (proteins and covariates)
  predictors <- as.matrix(train_data[, c(protein_names, "Age", "Sex", "Height", "Weight", "egfr_creat_epi_2021_baseline", "time_dif")])
  
  # Fit the elastic net model using cross-validation to find the best lambda
  cv_model <- cv.glmnet(predictors, response, alpha = 0.5, type.measure = 'mse')
  
  # Extract coefficients at the best lambda (1 standard error rule)
  best_lambda <- cv_model$lambda.min
  en_model <- glmnet(predictors, response, lambda = best_lambda, alpha = 0.5)
  
  # Extract coefficients for proteins in protein_names only
  coefficients <- coef(en_model)
  
  # Ensure coefficients are in a usable format
  coefficients <- as.matrix(coefficients)
  
  # Create a data frame of coefficients with the phenotype identifier
  coef_df <- data.frame(protein = rownames(coefficients), coefficient = as.vector(coefficients), phenotype = i, stringsAsFactors = FALSE)
  
  # Save coefficients to coef_results dataframe
  coef_results <- rbind(coef_results, coef_df)
  
  # Evaluate performance on test set using regression metrics
  test_y <- test_data[, i]
  test_x <- as.matrix(test_data[, c(protein_names, "Age", "Sex", "Height", "Weight", "egfr_creat_epi_2021_baseline", "time_dif")])
  pred <- predict(en_model, test_x, type = "response")
  
  # Calculate regression metrics
  mse_value <- mse(test_y, as.numeric(pred))
  mae_value <- mae(test_y, as.numeric(pred))
  r2_value <- cor(test_y, as.numeric(pred))^2
  
  # Save metrics to metrics_results dataframe
  metrics_df <- data.frame(phenotype = i, MSE = mse_value, MAE = mae_value, R_squared = r2_value)
  metrics_results <- rbind(metrics_results, metrics_df)
}

# Save dataframes to CSV files (optional)
write.csv(coef_results, "BF_coef_results.csv", row.names = FALSE)
write.csv(metrics_results, "BF_metrics_results.csv", row.names = FALSE)

# View the coefficient results
print(coef_results)

# View the metrics results
print(metrics_results)