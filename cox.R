library(survival)
KNN_data <- read_csv("/data/home/ha231442/HF/KNN_data.csv")
covar_all <- read_rds("/data/home/ha231442/HF/HF/ckd_ukb_all_covariates_NA.rds")


KNN_data <-  KNN_data %>%
  as.data.frame() %>%
  mutate(eid = as.character(eid)) %>% 
  select(- ...1)

covar_all <-  df %>%
  as.data.frame() %>%
  rename(eid = f.eid) %>%
  mutate(eid = as.character(eid)) 

cox_data <- merge(covar_all, KNN_data, by = "eid")

#remove prevalent HF
cox_data <- cox_data %>% filter(hf_prevalent != "Yes")
write.csv(cox_data, "cox_data.csv")


cox_data$ethnicity_binary <- as.factor(cox_data$ethnicity_binary)
cox_data$Sex <- as.factor(cox_data$Sex)

# simple cox ------------------------------------------------------------
##prepring for cox 
which(colnames(cox_data) == "a1bg") ##from col 63
which(colnames(cox_data) == "zpr1") ##to 2973
prot <- colnames(cox_data[, 63:2973])
matrix_res <- matrix(data = NA, nrow = length(prot), ncol = 4)
rownames(matrix_res) <- prot
colnames(matrix_res) <- paste(c(paste("coef"), paste("lower .95"), paste("upper .95"), paste("Pr(>|z|)")), sep = "")


## no protein adjustments to loook at apoe effect 
cox_model <- coxph(Surv(time_hf_event_date , event_hf_event_date) ~ age + factor(inflam_data$sex) + apoe + PLAUR ,data = cox_data, na.action = na.exclude)
cox_zph <- cox.zph(cox_model) ## to check redisuals 
cox_zph


for (i in prot) {
  str(i)
  cox_model <- coxph(Surv(time_hf_event_date, event_hf_event_date) ~ cox_data[, i] + Age + Sex + egfr_creat_epi_2021_baseline, data = cox_data, na.action = na.exclude)
  sum_gwa <- summary(cox_model) # Save summaries table
  coef_gwa <- sum_gwa$coefficients # Save coef table
  conf_gwa <- confint(cox_model, level = 0.95) # Save conf int table
  b_gwa <- coef_gwa[1, 2]  # Choose coef you want
  L95_gwa <- exp(conf_gwa[1, 1])
  U95_gwa <- exp(conf_gwa[1, 2])
  p_gwa <- coef_gwa[1, 5]
  # Add them in table
  matrix_res[i, 1] <- b_gwa
  matrix_res[i, 2] <- L95_gwa
  matrix_res[i, 3] <- U95_gwa
  matrix_res[i, 4] <- p_gwa
}
data_res <- as.data.frame(matrix_res)
data_res$fdr <- p.adjust(data_res$`Pr(>|z|)`, method = "fdr")
write.csv(data_res, "basic_cox.csv", row.names = TRUE)

data_res_filtered <- data_res %>% filter(fdr < 0.05)
write.csv(data_res_filtered, "basic_cox_FDR.csv", row.names = TRUE)

install.packages("synapser", repos=c("http://ran.synapse.org", "https://cloud.r-project.org"))
install.packages("synapser")

##model 2 
cox_data$smoking_binary <- as.factor(cox_data$smoking_binary)
cox_data$ethnicity_binary <- as.factor(cox_data$ethnicity_binary)
cox_data$ethnicity_binary <- as.factor(cox_data$ethnicity_binary)

for (i in prot) {
  str(i)
  cox_model <- coxph(Surv(time_hf_event_date, event_hf_event_date) ~ cox_data[, i] + Age + Sex + egfr_creat_epi_2021_baseline + smoking_binary + BMI + ethnicity_binary, data = cox_data, na.action = na.exclude)
  sum_gwa <- summary(cox_model) # Save summaries table
  coef_gwa <- sum_gwa$coefficients # Save coef table
  conf_gwa <- confint(cox_model, level = 0.95) # Save conf int table
  b_gwa <- coef_gwa[1, 2]  # Choose coef you want
  L95_gwa <- exp(conf_gwa[1, 1])
  U95_gwa <- exp(conf_gwa[1, 2])
  p_gwa <- coef_gwa[1, 5]
  # Add them in table
  matrix_res[i, 1] <- b_gwa
  matrix_res[i, 2] <- L95_gwa
  matrix_res[i, 3] <- U95_gwa
  matrix_res[i, 4] <- p_gwa
}
data_res <- as.data.frame(matrix_res)
data_res$fdr <- p.adjust(data_res$`Pr(>|z|)`, method = "fdr")
write.csv(data_res, "model2_cox.csv", row.names = TRUE)

data_res_filtered <- data_res %>% filter(fdr < 0.05)
write.csv(data_res_filtered, "model2_cox_FDR.csv", row.names = TRUE)

# forest plots ------------------------------------------------------------
top_20_fdr <- read_excel("FP_basic_cox.xlsx")
top_20_fdr$lower <- top_20_fdr$`lower .95` 
top_20_fdr$upper <- top_20_fdr$`upper .95` 
top_20_fdr$prot <- toupper(top_20_fdr$prot)

top_20_fdr$prot <- factor(top_20_fdr$prot, levels = rev(unique(top_20_fdr$prot)))

basic_model_top20 <- ggplot(data=top_20_fdr, aes(x=prot, y=coef, ymin= lower, ymax= upper)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Protein") + ylab("HR (95% CI)") +
  theme_bw() +  # use a white background
  ggtitle("Basic Cox PH") +  # Add your title here
  theme(axis.title.x = element_text(size = 18),  # adjust x-axis label size
        axis.title.y = element_text(size = 18),  # adjust y-axis label size
        axis.text.x = element_text(size = 18),   # adjust x-axis text size
        axis.text.y = element_text(size = 18),   # adjust y-axis text size
        plot.title = element_text(size = 20))  # adjust title size
basic_model_top20
ggsave("basic_model_top20.png")

top_20_fdr <- read_excel("FP_m2_cox.xlsx")
top_20_fdr$lower <- top_20_fdr$`lower .95` 
top_20_fdr$upper <- top_20_fdr$`upper .95` 
top_20_fdr$prot <- toupper(top_20_fdr$prot)
top_20_fdr$prot <- factor(top_20_fdr$prot, levels = rev(unique(top_20_fdr$prot)))

model2_top20 <- ggplot(data=top_20_fdr, aes(x=prot, y=coef, ymin= lower, ymax= upper)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Protein") + ylab("HR (95% CI)") +
  theme_bw() +  # use a white background
  ggtitle("Adjusted Cox PH") +  # Add your title here
  theme(axis.title.x = element_text(size = 18),  # adjust x-axis label size
        axis.title.y = element_text(size = 18),  # adjust y-axis label size
        axis.text.x = element_text(size = 18),   # adjust x-axis text size
        axis.text.y = element_text(size = 18),   # adjust y-axis text size
        plot.title = element_text(size = 20))  # adjust title size

model2_top20
ggsave("model2_top20.png")


# volcano plots -----------------------------------------------------------
###basic cox 
cox_basic <- read_excel("basic_cox.xlsx")
cox_basic <- cox_basic[order(cox_basic$bf), ]
cox_basic$log_fdr <- -log10(cox_basic$bf)
cox_basic$prot <- toupper(cox_basic$prot)
top_significant <- head(cox_basic[cox_basic$bf < 0.05, ], 10)


p <- ggplot(cox_basic, aes(x = coef, y = log_fdr, color = bf < 0.05)) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c("darkgray", "darkslateblue"), 
    name = "p-value", 
    labels = c("p-value < 0.05", "p-value ≥ 0.05")  # Adjust legend labels
  ) +
  geom_text_repel(data = top_significant, aes(label = prot), size = 5, show.legend = FALSE) +
  geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +  # Add vertical line at x = 1
  labs(x = "HR", y = "-log10(p-value)",  # Adjust axis labels
       title = "Basic Cox PH") +
  theme_minimal() +
  theme(
    plot.margin = margin(1, 1, 1, 1, "cm"),  # Adjust plot margins
    axis.title = element_text(size = 16),  # Increase axis label size
    axis.text = element_text(size = 16), 
    plot.title = element_text(size = 18),  # Increase title size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),   # Increase legend text size
    legend.key.size = unit(1.5, "lines")     # Increase legend key size
  ) + 
  coord_cartesian(ylim = c(min(cox_basic$log_fdr), max(cox_basic$log_fdr) + 1))  # Expand y-axis
# Save plot with increased dimensions
ggsave("volcano_plot.png", plot = p, width = 10, height = 8, units = "in")

###adjusted cox
model2_cox <- read_excel("adjusted_cox.xlsx")
model2_cox <- model2_cox[order(model2_cox$bf), ]
model2_cox$log_fdr <- -log10(model2_cox$bf)
model2_cox$prot <- toupper(model2_cox$prot)
top_significant <- head(model2_cox[model2_cox$bf < 0.05, ], 10)


p <- ggplot(model2_cox, aes(x = coef, y = log_fdr, color = bf < 0.05)) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c("darkgray", "darkslateblue"), 
    name = "p-value", 
    labels = c("p-value < 0.05", "p-value ≥ 0.05")  # Adjust legend labels
  ) +
  geom_text_repel(data = top_significant, aes(label = prot), size = 5, show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +  # Add vertical line at x = 1
  labs(x = "HR", y = "-log10(p-value)",  # Adjust axis labels
       title = "Adjusted Cox PH") +
  theme_minimal() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),  # Adjust plot margins
        axis.title = element_text(size = 16),  # Increase axis label size
        axis.text = element_text(size = 16), 
        plot.title = element_text(size = 18)) +  # Increase axis text size
  coord_cartesian(ylim = c(min(model2_cox$log_fdr), max(model2_cox$log_fdr) + 1))  # Expand y-axis
# Save plot with increased dimensions
ggsave("volcano_plot_model2_cox.png", plot = p, width = 10, height = 8, units = "in")

