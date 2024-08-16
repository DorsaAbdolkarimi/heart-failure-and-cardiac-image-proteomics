install.packages("BiocManager")
BiocManager::install("impute")
install.packages("tidyverse")
library(readr)
library(tidyverse)
library(impute)
library(ggplot2)


proteomic <- read_csv("/data/home/ha231442/HF/data.csv")
##filtering 
proteomic_filtered <- proteomic %>%
  select(where(~ mean(is.na(.)) <= 0.2)) ##removed the said proteins 

# Remove rows with more than 20% NAs
proteomic_filtered <- proteomic_filtered[apply(proteomic_filtered, 1, function(row) mean(is.na(row)) <= 0.2), ] ##44191 people and 2911 proteisn 
write.csv(proteomic_filtered, "proteomic_filtered.csv")

## now let's impute with the 2 different methods and compare
set.seed(2)
prot <- proteomic_filtered[, 2:2912] %>%
  as.matrix()
eid <- proteomic_filtered[, 1]
KNN <- impute.knn(prot, k = 10)
KNN_prot <- as.data.frame(KNN$data)
KNN_data <- cbind(eid, KNN_prot)
KNN_data <- read_csv("/data/home/ha231442/HF/KNN_data.csv")


# Calculate the number of NAs in each column and sort them
sorted_na_counts <- sort(colSums(is.na(prot)), decreasing = TRUE) 
write_csv(sorted_na_counts, "sorted_na_counts.csv")

head(sorted_na_counts, 20)
# Get the column names with the highest NA counts
top_columns <- names(head(sorted_na_counts, 100))

# Find the original column numbers in proteomic_filtered
columns_to_compare <- match(top_columns, names(proteomic_filtered))

# Display the column names and their original indices
df <- data.frame(top_columns= top_columns, columns_to_compare = columns_to_compare)


###Density plots 
# Create a directory to store the plots if it doesn't exist
dir.create("density_plots", showWarnings = FALSE)

# Iterate over the column names in df$top_columns
for (top_columns in df$top_columns) {
  if (top_columns %in% colnames(prot) && top_columns %in% colnames(KNN_prot)) {
    # Print the current column being processed (optional)
    cat("Processing column:", top_columns, "\n")
    
    # Open a PDF file to save the plot
    pdf(file = paste0("density_plots/", top_columns, "_density_plot.pdf"))
    
    # Generate the density plot for the original data
    plot(density(prot[[top_columns]], na.rm = TRUE), 
         main = paste("Density Plot of", top_columns), 
         xlab = top_columns, 
         ylab = "Density")
    
    # Add the density plot for the imputed data
    lines(density(KNN_prot[[top_columns]], na.rm = TRUE), col = "red", lty = 3, lwd = 2)
    
    # Add a legend to the plot
    legend("topright", legend = c("Original", "Imputed"), 
           col = c("black", "red"), lty = c(1, 3), lwd = c(1, 2))
    
    # Close the PDF file
    dev.off()
  } else {
    cat("Skipping column:", top_columns, " - not present in both data frames\n")
  }
}


###summary stats 
summary_original <- summary(prot[, top_columns]) %>%
  as.data.frame() %>%
  select(-Var1) %>%
  mutate(Var2 = as.character(Var2)) %>%
  separate(Freq, into = c("Statistic", "Value"), sep = ":") %>%
  mutate(across(c(Statistic, Value), str_trim)) %>%
  pivot_wider(names_from = Statistic, values_from = Value) %>%
  distinct(Var2, .keep_all = TRUE) %>%
  mutate(Var2 = factor(Var2, levels = unique(Var2))) %>%
  column_to_rownames(var = "Var2") %>%
  rownames_to_column(var = "prot")
write_csv(summary_original, "summary_original_100.csv")


summary_KNN <- summary(KNN_prot[, top_columns]) %>%
  as.data.frame() %>%
  mutate(Var2 = as.character(Var2)) %>%
  separate(Freq, into = c("Statistic", "Value"), sep = ":") %>%
  mutate(across(c(Statistic, Value), str_trim)) %>%
  pivot_wider(names_from = Statistic, values_from = Value) %>%
  distinct(Var2, .keep_all = TRUE) %>%
  select(-Var1) %>%
  column_to_rownames(var = "Var2") %>%
  rownames_to_column(var = "prot")
write_csv(summary_KNN, "summary_KNN_100.csv")




####Kolmogorov–Smirnov statistic
# Initialize an empty vector to store the p-values
ks_test_results <- numeric(length(columns_to_compare))

# Loop through each column index in columns_to_compare
for (i in seq_along(columns_to_compare)) {
  col_index <- columns_to_compare[i]
  
  # Perform the KS test and store the p-value
  ks_test_results[i] <- ks.test(prot[, col_index], KNN_prot[, col_index])$p.value
}

# Print the results
test_results <- data.frame(col = columns_to_compare, ks_test_results = ks_test_results)
write_csv(test_results, "ks_test_results_KNN.csv")

####2-Wasserstein Distance
# checking impp quality ---------------------------------------------------
library(transport)

# Function to calculate 2-Wasserstein Distance for a column
wasserstein_distance <- function(prot, KNN_prot) {
  return(wasserstein1d(prot, KNN_prot, p = 2))
}

# Calculate the distance for each column
distances <- sapply(df$top_columns, function(col) {
  wasserstein_distance(prot[[col]], KNN_prot[[col]])
})

distances <- distances %>%
  as.data.frame() %>%
  rownames_to_column()
write_csv(distances, "wasserstein_distance_KNN.csv")

####3- KL divergence 
install.packages("entropy")
library(entropy)

# Function to calculate KL Divergence for a column
kl_divergence <- function(original, imputed) {
  # Remove NA values
  original <- na.omit(original)
  imputed <- na.omit(imputed)
  
  # Estimate densities
  original_density <- density(original)
  imputed_density <- density(imputed)
  
  # Interpolate densities
  interp_density <- approx(original_density$x, original_density$y, xout = imputed_density$x)$y
  
  # Ensure no zero densities to avoid log(0)
  interp_density[interp_density == 0] <- .Machine$double.eps
  imputed_density$y[imputed_density$y == 0] <- .Machine$double.eps
  
  # Calculate KL Divergence
  kl <- KL.empirical(interp_density, imputed_density$y)
  return(kl)
}

# Calculate the divergence for each column
divergences <- sapply(df$top_columns, function(col) {
  kl_divergence(prot[[col]], KNN_prot[[col]])
})
# Calculate the divergence for each column
divergences <- sapply(df$top_columns, function(col) {
  kl_divergence(prot[[col]], KNN_prot[[col]])
})

divergences <- divergences %>%
  as.data.frame() %>%
  rownames_to_column()
write_csv(divergences, "KL_divergences_KNN.csv")