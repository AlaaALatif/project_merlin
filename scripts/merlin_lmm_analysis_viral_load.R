# Install required libraries 
required_packages <- c("lmerTest", "rms", "emmeans", "pbkrtest")

# Install packages that are not yet installed
installed_packages <- rownames(installed.packages())
for(pkg in required_packages) {
  if(!(pkg %in% installed_packages)) {
    install.packages(pkg)
  }
}
# Load required libraries
library(ggplot2)
library(lmerTest)    # For lmer() with p-values
library(rms)         # For restricted cubic spline (rcs)
library(emmeans)     # For estimated marginal means and pairwise contrasts

# Ensure Treatment is a factor and Time is numeric.
# (Assume your data is stored in the data frame 'df'.)
df = read.csv("/Users/alaa/Documents/ucsf/data/suliman/merlin/meta/experiment2_viral_load_lmm_data.csv")
head(df)
df$treatment_group <- factor(df$treatment_group, levels = c("Immediate", "Deferred"))
df$timepoint <- as.numeric(df$timepoint)  # e.g., 3, 6, 12, 24

# Fit a linear mixed model with a random intercept for PID and RCS for Time.
# Use rcs(Time, 5) for 5 knots (3 internal + 2 boundary knots).
model <- lmer(log10_viral_load ~ treatment_group * factor(timepoint) + (1 | patient_id), data = df)

# Show model summary to inspect fixed effects (including interactions)
summary(model)

# Extract the summary of the model
model_summary <- summary(model)

# Extract fixed effects coefficients as a data frame
fixed_effects <- as.data.frame(coef(model_summary))

# Save the fixed effects table to CSV
write.csv(fixed_effects, "/Users/alaa/Documents/ucsf/data/suliman/merlin/phenotype/experiment2_viral_load_lmm_fixed_effects_summary_v2.csv", row.names = TRUE)

# Extract random effects variance components as a data frame
random_effects <- as.data.frame(VarCorr(model))

# Save random effects table to CSV
write.csv(random_effects, "/Users/alaa/Documents/ucsf/data/suliman/merlin/phenotype/experiment2_viral_load_lmm_random_effects_summary_v2.csv", row.names = FALSE)

# Use ANOVA (with Satterthwaite approximation) to test the overall significance
# of the interaction between Treatment and Time.
anova(model)
# Run ANOVA on the model
anova_results <- anova(model)

# Convert the ANOVA results to a data frame
anova_df <- as.data.frame(anova_results)

# Save the ANOVA table to a CSV file (including row names for clarity)
write.csv(anova_df, "/Users/alaa/Documents/ucsf/data/suliman/merlin/phenotype/experiment2_viral_load_lmm_anova_results_v2.csv", row.names = TRUE)

# For plotting, we want time to be numeric.
# Here we assume that the factor levels (e.g., "0", "1", "2", "3") represent the timepoints.
df$timepoint_num <- as.numeric(as.character(df$timepoint))

# Compute estimated marginal means.
# Since the model was fit with factor(timepoint), use the same levels.
emm_df <- as.data.frame(emmeans(model, ~ treatment_group | timepoint))
# Convert timepoint from factor to numeric for plotting
emm_df$timepoint_num <- as.numeric(as.character(emm_df$timepoint))

# Create the spaghetti plot
p <- ggplot() +
  # Plot individual trajectories (spaghetti lines)
  geom_line(data = df, 
            aes(x = timepoint_num, y = log10_viral_load, 
                group = patient_id, color = treatment_group),
            alpha = 0.3) +
  # Overlay the estimated marginal means as points
  geom_point(data = emm_df, 
             aes(x = timepoint_num, y = emmean, color = treatment_group),
             size = 3) +
  # Connect the estimated marginal means with lines
  geom_line(data = emm_df, 
            aes(x = timepoint_num, y = emmean, color = treatment_group, group = treatment_group),
            size = 1.5) +
  labs(x = "Timepoint", y = "Log Viral Load",
       title = "Spaghetti Plot with Estimated Marginal Means",
       color = "Treatment Group") +
  theme_minimal()

# Print the plot
print(p)
