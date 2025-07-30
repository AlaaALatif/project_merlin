library(tibble)
library(ggplot2)
library(lme4)
library(lmerTest)   # lme4 + p‑values
library(emmeans)
library(dplyr)

# ------------------------------------------------------------------------------
#  RISK‑6 signature scorer – official method
# ------------------------------------------------------------------------------

risk6_score <- function(expr,
                        input = c("rna", "ct"),   # use "ct" when the matrix holds raw Ct values
                        pseudocount = 1) {
  
  input <- match.arg(input)
  
  # 1.  Canonical gene list (HGNC symbols; FCGR1B was renamed → FCGR1BP)
  up   <- c("GBP2", "FCGR1BP", "SERPING1")     # up‑regulated in TB progressors
  down <- c("TUBGCP6", "TRMT2A", "SDR39U1")    # down‑regulated in progressors
  genes <- c(up, down)
  
  # 2.  Sanity checks ----------------------------------------------------------
  if (!all(genes %in% rownames(expr))) {
    stop("Missing gene(s): ",
         paste(setdiff(genes, rownames(expr)), collapse = ", "),
         call. = FALSE)
  }
  
  # 3.  Subset and coerce to matrix -------------------------------------------
  m <- as.matrix(expr[genes, , drop = FALSE])
  
  # 4.  Pre‑processing ---------------------------------------------------------
  #     • RNA mode: if values look like raw counts, apply log2(x + 1)
  #     • Ct mode: use values as‑is
  if (input == "rna" && max(m, na.rm = TRUE) > 100) {
    m <- log2(m + pseudocount)
  }
  
  # 5.  Core calculation -------------------------------------------------------
  col_gmean <- function(x) exp(colMeans(log(x), na.rm = TRUE))
  
  if (input == "rna") {
    # RISK6_geo  (log2‑ratio of geometric means)
    gm_up   <- col_gmean(m[up,   , drop = FALSE])
    gm_down <- col_gmean(m[down, , drop = FALSE])
    score   <- log2(gm_up) - log2(gm_down)      # positive ⇒ higher TB risk
  } else {  # Ct data (lower Ct = higher expression ⇒ sign flips)
    ct_up_mean   <- colMeans(m[up,   , drop = FALSE], na.rm = TRUE)
    ct_down_mean <- colMeans(m[down, , drop = FALSE], na.rm = TRUE)
    score <- ct_down_mean - ct_up_mean          # positive ⇒ higher TB risk
  }
  
  # 6.  Return tidy tibble -----------------------------------------------------
  tibble::tibble(sample_id = colnames(expr),
                 RISK6      = unname(score))
}


#  1. Read your gene × sample table (TSV with 'gene_symbol' first column)
expr_file <- "~/Downloads/formatted_counts_by_gene_2024-12-01.tsv"

expr_mat <- readr::read_tsv(expr_file, show_col_types = FALSE) %>% 
  tibble::column_to_rownames("gene_symbol") %>% 
  as.matrix()

# 2. Score each sample (auto‑detects raw vs. logged RNA‑seq values)
risk6_tbl <- risk6_score(expr_mat, input = "rna")
print(risk6_tbl, n = 5)

# Read metadata
meta_path  <- "~/Documents/ucsf/data/suliman/merlin/meta/metadata_round1_2025-07-21.csv"      # adjust as needed
meta_df    <- read_csv(meta_path, show_col_types = FALSE)
# Filter metadata for the desired timepoints
desired_times <- c("BL", "X-1", "DX", "M3", "M6", "M12")
filtered_df <- meta_df %>% 
  filter(timepoint %in% desired_times)
results_df <- risk6_tbl %>%
  inner_join(filtered_df, by = "sample_id")
# Create new field for patient ID
results_df <- results_df %>%
  mutate(
    patient_id = sub("_.*$", "", sample_id)   # drop everything from first "_" to end
  )
# Create new time and treatment fields with appropriate factorization
results_df <- results_df %>%
  mutate(
    time = factor(timepoint,
                       levels = c("BL","X-1","DX","M3","M6","M12")),  # BL = reference
    treatment = factor(art_treatment_group, levels = c("Immediate","Deferred"))
  )
# Flag anything that is not a finite number
results_df <- results_df %>%
  mutate(
    RISK6 = ifelse(!is.finite(RISK6), NA_real_, RISK6)   # turn NaN / ±Inf into NA
  )
# Drop samples with non-finite risk scores
results_df <- results_df %>%
  filter(!is.na(RISK6))
# Basic visualization of RISK6 scores
# pre‑calculate the dodge & jitter position objects (makes the code tidier)
pos_dodge <- position_dodge(width = 0.75)
pos_jit   <- position_jitterdodge(
  jitter.width = 0.20,   # horizontal jitter for individual points
  jitter.height = 0,     # no vertical jitter
  dodge.width  = 0.75)   # align with the box‑plot dodge
# generate box plot
ggplot(results_df,
       aes(time, RISK6,
           fill   = treatment,
           colour = treatment)) +
  
  geom_boxplot(
    outlier.shape = NA,
    alpha         = 0.6,
    width         = 0.65,
    position      = pos_dodge) +
  
  geom_jitter(
    size      = 1.1,
    alpha     = 0.35,
    position  = pos_jit) +          # ← keep only this positional spec
  #   (no width/height arguments here)
  
  labs(
    x     = "Visit",
    y     = "RISK6 score (raw)",
    fill  = "ART start",
    colour= "ART start",
    title = "Raw RISK6 score distribution across timepoints") +
  
  theme_classic(base_size = 12) +
  theme(legend.position = "top",
        plot.title.position = "plot")


# Fit a linear fixed effects model
risk_mod <- lmer(
  RISK6 ~ time * treatment + (1 | patient_id),
  data = results_df
)
anova(risk_mod, type = 3)        # Type III tests with Satterthwaite df
summary(risk_mod)$coefficients   # fixed‑effect estimates
# Rankings per timepoint
tp_means <- emmeans(risk_mod, ~ time)
pairs(tp_means, adjust = "tukey")   # all pairwise comparisons
summary(tp_means, infer = TRUE)     # CIs & ranking
# Contrasting each timepoint to baseline within each treatment group
tp_by_trt <- emmeans(risk_mod, ~ time | treatment)
bl_vs_later <- contrast(
  tp_by_trt,
  method = "trt.vs.ctrl",
  ref = "BL",
  adjust = "BH"    # FDR control, optional
)
bl_vs_later
# Comparing treatment groups at each timepoint
# Immediate – Deferred at each timepoint
tp_trt_emm <- emmeans(risk_mod, ~ treatment | time)
trt_vs_trt <- contrast(tp_trt_emm,
                       method   = "pairwise",
                       by       = "time",
                       adjust   = "BH")   # BH or "none" if you prefer raw p's
trt_vs_trt
# Visualizing fitted trajectories
# 2.1  Prediction grid
newdat <- expand.grid(
  time = levels(results_df$time),
  treatment = levels(results_df$treatment)
)

# keep factor structure identical to model
newdat <- newdat %>% mutate(
  time = factor(time, levels = levels(results_df$time)),
  treatment = factor(treatment, levels = levels(results_df$treatment))
)

# 2.2  Point predictions (no random effects)
newdat$pred <- predict(risk_mod, newdata = newdat, re.form = NA)

# 2.3  Standard errors via the fixed‑effect variance–covariance matrix
X   <- model.matrix(~ time * treatment, newdat)
vc  <- vcov(risk_mod)                # fixed‑effects V‑C matrix
newdat$se <- sqrt(diag(X %*% vc %*% t(X)))
newdat <- newdat %>%
  mutate(lwr = pred - 1.96 * se,
         upr = pred + 1.96 * se)

# 2.4  Plot
ggplot(newdat,
       aes(time, pred,
           colour = treatment, group = treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                width = 0.2) +
  labs(x = "Timepoint",
       y = "Predicted RISK6 (population level)",
       colour = "ART start") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")