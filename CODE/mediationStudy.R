##################################################################################################
# Script Mediation-by-time (GLMM + LMM) with subject-level bootstrap
# Author: Claudio Aleman Morillo    email: caleman@us.es
# Version: 1.0
# -----------------------------------------------------------
# Purpose
#     1) fit two mixed-effects models (Poisson GLMM for a symptom outcome; LMM for a brain-area mediator)
#     2) perform subject-level nonparametric bootstrap using identical resamples across both models
#     3) compute time-dependent direct, mediated, and total effects at three representative points: EARLY, MID, LATE
#     4) export bootstrap draws and summary tables to ./outputs/
#
#
# Reproducibility
#   Tested with: R >= 4.2, lme4, lmerTest, dplyr, tidyr
# -----------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)       # glmer/lmer engines
  library(lmerTest)   # p-values for lmer (not required for glmer)
})

options(stringsAsFactors = FALSE)

# -----------------------------
# 0) Configuration
# -----------------------------
# INPUT: A SINGLE CSV with all needed columns.
# EXPECTED COLUMNS (standardized):
#   SUBJECT_ID                : subject identifier (character)
#   Assessment                : visit index (integer)
#   Sex                       : "Male"/"Female" (factor with two levels)
#   Age_MRI                   : age at MRI / assessment (numeric)
#   Age_Inclusion             : baseline age (numeric)
#   eTIV                      : estimated intracranial volume (numeric)
#   CPZ_equivalent            : medication dose (chlorpromazine equivalents; numeric)
#   Symptom                   : symptom outcome (count-like; numeric, e.g., BPRS total)
#   AreaValue                 : mediator value for the brain region of interest (numeric)
#
# You can override CONFIG via environment variables if desired.
CONFIG <- list(
  INPUT_DATA       = Sys.getenv("INPUT_DATA", unset = "data/analysis_data.csv"),
  SYMPTOM_COL      = Sys.getenv("SYMPTOM_COL", unset = "Symptom"),
  AREA_COL         = Sys.getenv("AREA_COL",    unset = "AreaValue"),
  SEX_COL          = Sys.getenv("SEX_COL",     unset = "Sex"),
  AGE_MRI_COL      = Sys.getenv("AGE_MRI_COL", unset = "Age_MRI"),
  AGE_INC_COL      = Sys.getenv("AGE_INC_COL", unset = "Age_Inclusion"),
  ETIV_COL         = Sys.getenv("ETIV_COL",    unset = "eTIV"),
  CPZ_COL          = Sys.getenv("CPZ_COL",     unset = "CPZ_equivalent"),
  ID_COL           = Sys.getenv("ID_COL",      unset = "SUBJECT_ID"),
  ASSESS_COL       = Sys.getenv("ASSESS_COL",  unset = "Assessment"),
  
  # analysis limits and bootstrap settings
  MAX_ASSESSMENT   = as.numeric(Sys.getenv("MAX_ASSESSMENT", unset = 10)),  # keep assessments <= this value
  B_BOOT           = as.integer(Sys.getenv("B_BOOT",         unset = 500)), # number of bootstrap replicates
  MAX_RETRIES      = as.integer(Sys.getenv("MAX_RETRIES",    unset = 200)),
  RNG_SEED         = as.integer(Sys.getenv("RNG_SEED",       unset = 123)),
  
  # Output directory
  OUT_DIR          = Sys.getenv("OUT_DIR", unset = "outputs")
)

# Make sure output dir exists
if (!dir.exists(CONFIG$OUT_DIR)) dir.create(CONFIG$OUT_DIR, recursive = TRUE)

# -----------------------------
# 1) Load and preprocess data
# -----------------------------
stopifnot(file.exists(CONFIG$INPUT_DATA))
D <- read.csv(CONFIG$INPUT_DATA)

# Keep only relevant columns and drop incomplete rows
vars_needed <- c(CONFIG$SYMPTOM_COL, CONFIG$AGE_INC_COL, CONFIG$SEX_COL, CONFIG$ETIV_COL,
                 CONFIG$AREA_COL, "Treatment_Time", CONFIG$CPZ_COL, CONFIG$ID_COL,
                 CONFIG$AGE_MRI_COL, CONFIG$ASSESS_COL)

# If Treatment_Time is not present, compute it from Age_MRI - Age_Inclusion
if (!"Treatment_Time" %in% names(D)) {
  stopifnot(all(c(CONFIG$AGE_MRI_COL, CONFIG$AGE_INC_COL) %in% names(D)))
  D$Treatment_Time <- D[[CONFIG$AGE_MRI_COL]] - D[[CONFIG$AGE_INC_COL]]
}

# Filter by maximum assessment index (if present)
if (CONFIG$ASSESS_COL %in% names(D)) {
  D <- D %>% filter(.data[[CONFIG$ASSESS_COL]] <= CONFIG$MAX_ASSESSMENT)
}

# Type casting and scaling if it is necessary
D <- D %>%
  transmute(
    Subject        = as.factor(.data[[CONFIG$ID_COL]]),
    Assessment     = if (CONFIG$ASSESS_COL %in% names(.)) as.numeric(.data[[CONFIG$ASSESS_COL]]) else NA_real_,
    sex            = as.factor(.data[[CONFIG$SEX_COL]]),
    Age_inclusion  = as.numeric(.data[[CONFIG$AGE_INC_COL]]),
    Age_MRI        = as.numeric(.data[[CONFIG$AGE_MRI_COL]]),
    residual_etiv  = as.numeric(scale(residuals(lm(
      as.numeric(.data[[CONFIG$ETIV_COL]]) ~ as.numeric(.data[[CONFIG$AGE_MRI_COL]]) + as.factor(.data[[CONFIG$SEX_COL]])
    )))),
    Treatment_Time = as.numeric(.data[["Treatment_Time"]]),
    CPZ_equivalent = as.numeric(.data[[CONFIG$CPZ_COL]]),
    symptom        = as.numeric(.data[[CONFIG$SYMPTOM_COL]]),
    area           = as.numeric(.data[[CONFIG$AREA_COL]])
  ) %>%
  filter(complete.cases(.)) %>%
  mutate(
    # Standardize continuous covariates used in models
    Age_inclusion  = as.numeric(scale(Age_inclusion)),
    Age_MRI        = as.numeric(scale(Age_MRI)),
    Treatment_Time = as.numeric(scale(Treatment_Time)),
    CPZ_equivalent = as.numeric(scale(CPZ_equivalent)),
    # Normalize sex levels to ensure consistent ordering
    sex            = factor(sex, levels = c("Male", "Female"))
  ) %>%
  as.data.frame()

# -----------------------------
# 2) Model formulas
# -----------------------------
# Outcome (symptom) model: Poisson GLMM with log link
form_Y <- as.formula(
  "symptom ~ 1 + Age_inclusion + sex + residual_etiv + area + Treatment_Time +Treatment_Time:CPZ_equivalent + CPZ_equivalent + area:Treatment_Time + area:CPZ_equivalent + (1|Subject)"
)

# Mediator (brain area) model: LMM
form_M <- as.formula(
  "area ~ Treatment_Time + residual_etiv + Age_inclusion + sex + CPZ_equivalent +CPZ_equivalent:Treatment_Time + (1|Subject)"
)

# Fit base models (sanity check)
fit_Y <- glmer(form_Y, data = D, family = poisson(link = "log"),
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
fit_M <- lmer(form_M,  data = D,
              control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

# -----------------------------
# 3) Subject-level bootstrap
# -----------------------------
set.seed(CONFIG$RNG_SEED)
B      <- CONFIG$B_BOOT
ids    <- levels(D$Subject)
K      <- length(ids)

# Helper checks
is_rank_deficient <- function(fm, data) {
  X <- model.matrix(lme4::nobars(fm), data)
  qr(X)$rank < ncol(X)
}

has_aliased_fixef <- function(mod) any(is.na(lme4::fixef(mod)))

ok_mermod <- function(mod) {
  conv_ok <- isTRUE(is.null(mod@optinfo$conv$lme4$messages)) ||
    !any(grepl("failed|converge", mod@optinfo$conv$lme4$messages, ignore.case = TRUE))
  !lme4::isSingular(mod, tol = 1e-6) && conv_ok && !has_aliased_fixef(mod)
}

# Resample-by-subject builder
build_boot_dataset <- function(id_sample, base = D) {
  out <- dplyr::bind_rows(lapply(seq_along(id_sample), function(j) {
    id <- id_sample[j]
    block <- base[base$Subject == id, , drop = FALSE]
    if (nrow(block) == 0L) return(NULL)
    block$Subject <- factor(paste0(id, "_b", j)) # unique cluster per draw
    block
  }))
  droplevels(out)
}

pY <- length(fixef(fit_Y))
boot_fixef_Y <- matrix(NA_real_, nrow = B, ncol = pY,
                       dimnames = list(paste0("rep", seq_len(B)), names(fixef(fit_Y))))

pM <- length(fixef(fit_M))
boot_fixef_M <- matrix(NA_real_, nrow = B, ncol = pM,
                       dimnames = list(paste0("rep", seq_len(B)), names(fixef(fit_M))))

conv_ok <- logical(B)

b <- 1L
while (b <= B) {
  tries <- 0L
  repeat {
    tries <- tries + 1L
    ids_b <- sample(ids, size = K, replace = TRUE)
    dat_b <- build_boot_dataset(ids_b)
    
    if (is_rank_deficient(form_Y, dat_b) || is_rank_deficient(form_M, dat_b)) {
      if (tries >= CONFIG$MAX_RETRIES) break else next
    }
    
    mY <- try(glmer(form_Y, data = dat_b, family = poisson(link = "log"),
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))), silent = TRUE)
    if (inherits(mY, "try-error")) { if (tries >= CONFIG$MAX_RETRIES) break else next }
    
    mM <- try(lmer(form_M, data = dat_b,
                   control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))), silent = TRUE)
    if (inherits(mM, "try-error")) { if (tries >= CONFIG$MAX_RETRIES) break else next }
    
    if (!ok_mermod(mY) || !ok_mermod(mM)) { if (tries >= CONFIG$MAX_RETRIES) break else next }
    
    # Align and store
    cfY <- fixef(mY); cfM <- fixef(mM)
    boot_fixef_Y[b, names(cfY)] <- cfY
    boot_fixef_M[b, names(cfM)] <- cfM
    
    conv_ok[b] <- TRUE
    message(sprintf("Bootstrap replicate %d succeeded in %d attempt(s).", b, tries))
    break
  }
  if (!conv_ok[b]) warning(sprintf("No valid replicate for b=%d after %d attempts; leaving NA.", b, CONFIG$MAX_RETRIES))
  b <- b + 1L
}

# Save raw bootstrap draws
sym_name  <- CONFIG$SYMPTOM_COL
area_name <- CONFIG$AREA_COL

write.csv(boot_fixef_Y, file.path(CONFIG$OUT_DIR, sprintf("bootstrap_fixef_Y_%s_%s.csv", sym_name, area_name)), row.names = TRUE)
write.csv(boot_fixef_M, file.path(CONFIG$OUT_DIR, sprintf("bootstrap_fixef_M_%s_%s.csv", sym_name, area_name)), row.names = TRUE)

# -----------------------------
# 4) Summaries and mediated-effect decomposition
# -----------------------------
summary_boot <- function(boot_mat, ok = NULL, probs = c(.025, .5, .975)) {
  if (!is.null(ok)) boot_mat <- boot_mat[ok, , drop = FALSE]
  boot_mat <- boot_mat[rowSums(is.na(boot_mat)) < ncol(boot_mat), , drop = FALSE]
  est_mean   <- apply(boot_mat, 2, function(x) mean(x, na.rm = TRUE))
  est_median <- apply(boot_mat, 2, function(x) median(x, na.rm = TRUE))
  est_se     <- apply(boot_mat, 2, function(x) sd(x, na.rm = TRUE))
  ci_lo      <- apply(boot_mat, 2, function(x) quantile(x, probs = probs[1], na.rm = TRUE, names = FALSE))
  ci_hi      <- apply(boot_mat, 2, function(x) quantile(x, probs = probs[3], na.rm = TRUE, names = FALSE))
  n_eff      <- apply(boot_mat, 2, function(x) sum(is.finite(x)))
  data.frame(coef = colnames(boot_mat), mean = est_mean, median = est_median,
             se = est_se, q2.5 = ci_lo, q97.5 = ci_hi, n_boot = n_eff,
             row.names = NULL, check.names = FALSE)
}

summ_stats <- function(x) c(mean = mean(x), l95 = quantile(x, .025), u95 = quantile(x, .975))

# Time points (labels only). Values are derived from scaled Treatment_Time: EARLY=min, MID=0, LATE=max.
TT_vals <- c(
  EARLY = min(D$Treatment_Time, na.rm = TRUE),
  MID   = 0,                          # since Treatment_Time is scaled
  LATE  = max(D$Treatment_Time, na.rm = TRUE)
)

calc_time_tables <- function(bootY, bootM, tt_vals) {
  bootY <- as.data.frame(bootY); bootM <- as.data.frame(bootM)
  
  A1 <- bootY$area
  A2 <- bootY$`area:Treatment_Time`
  A3 <- bootY$`Treatment_Time:CPZ_equivalent`
  A4 <- bootY$CPZ_equivalent
  
  
  B1 <- bootM$CPZ_equivalent
  B2 <- bootM$`Treatment_Time:CPZ_equivalent`
  
  res <- list()
  for (lbl in names(tt_vals)) {
    t  <- tt_vals[[lbl]]
    # Medication -> Area (path a)
    pathA <- B1 + B2 * t
    # Area   -> Symptom given time (path b)
    pathB <- A1 + A2 * t
    # Direct effect of medication on symptom at time t (path c')
    direct <- A4 + A3 * t
    # Indirect (mediated) effect via Area
    mediated <- (pathA) * (pathB)
    total <- mediated + direct
    
    df <- data.frame(
      time_label = lbl,
      rbind(
        data.frame(effect = "PathA_Medication_to_Area",    t(summ_stats(pathA))),
        data.frame(effect = "PathB_Area_to_Symptom",       t(summ_stats(pathB))),
        data.frame(effect = "Direct",                       t(summ_stats(direct))),
        data.frame(effect = "Mediated",                     t(summ_stats(mediated))),
        data.frame(effect = "Total",                        t(summ_stats(total)))
      ), row.names = NULL
    )
    res[[lbl]] <- df
  }
  do.call(rbind, res)
}

summary_Y <- summary_boot(boot_fixef_Y, ok = conv_ok)
summary_M <- summary_boot(boot_fixef_M, ok = conv_ok)

med_table <- calc_time_tables(boot_fixef_Y[conv_ok, , drop = FALSE],
                              boot_fixef_M[conv_ok, , drop = FALSE],
                              TT_vals)

# Save summaries
write.csv(summary_Y, file.path(CONFIG$OUT_DIR, sprintf("summary_fixef_Y_%s_%s.csv", sym_name, area_name)), row.names = FALSE)
write.csv(summary_M, file.path(CONFIG$OUT_DIR, sprintf("summary_fixef_M_%s_%s.csv", sym_name, area_name)), row.names = FALSE)
write.csv(med_table, file.path(CONFIG$OUT_DIR, sprintf("mediation_summary_%s_%s.csv", sym_name, area_name)), row.names = FALSE)
