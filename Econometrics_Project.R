# 0. Setup --------------------------------------------------------------------
# All library() calls live here so the script has a single dependency surface.
# install.packages(c("sandwich", "lmtest", "MASS", "mfx", "aod", "dplyr",
#                    "logistf", "stargazer", "car", "corrplot"))

library(sandwich)    # robust / cluster-robust variance-covariance matrices
library(lmtest)      # coeftest(), lrtest()
library(MASS)
library(mfx)         # logitmfx(): marginal effects for logit
library(aod)         # wald.test()
library(dplyr)
library(logistf)     # Firth's bias-reduced logit (for quasi-separation)
library(stargazer)   # publication-quality summary / regression tables
library(car)         # vif()
library(corrplot)    # correlation heatmaps
library(pROC)        # roc(), auc(), coords() for ROC curve


# 1. Data Loading & Cleaning --------------------------------------------------
# Read the raw Mass Mobilization file, keep only protest events, and build the
# analytic features used by the model (size, duration, demand dummies, lags).

## 1.1 Load & filter ----------------------------------------------------------
data <- read.csv("dataverse_files/mmALL_073120_csv.csv", stringsAsFactors = FALSE)
str(data)

# Non-protest rows have no meaningful DV (protesterviolence is undefined there).
data <- data[data$protest == 1, ]

data$country <- factor(data$country)
data$region  <- factor(data$region)
# Region codes (left as numeric labels):
# 1 South America, 2 Central America, 3 North America, 4 Europe,
# 5 Asia, 6 MENA, 7 Africa

## 1.2 Dates & duration -------------------------------------------------------
# Build proper Date objects and compute event length in days.
data$start_date <- as.Date(paste(data$startyear, data$startmonth, data$startday, sep = "-"),
                           format = "%Y-%m-%d")
data$end_date   <- as.Date(paste(data$endyear,   data$endmonth,   data$endday,   sep = "-"),
                           format = "%Y-%m-%d")

data$duration_days <- as.numeric(data$end_date - data$start_date)
data$duration_days[data$duration_days < 0] <- 0  # guard against bad rows

## 1.3 Protest size -----------------------------------------------------------
# Two source columns exist: a bucketed `participants_category` and a free-text
# `participants` field. We harmonise both into a single numeric `protest_size`,
# log-transform it (heavy right skew), and impute missing values by country
# median while flagging the imputation.

category_midpoints <- c(
  "50-99"      = 75,
  "100-999"    = 550,
  "1000-1999"  = 1500,
  "2000-4999"  = 3500,
  "5000-10000" = 7500,
  ">10000"     = 15000
)
data$size_from_category <- category_midpoints[data$participants_category]

# Parse the free-text column: handles plain numbers and "1000s" / "100s" suffixes.
data$size_from_raw <- sapply(data$participants, function(x) {
  if (is.na(x) || x == "") return(NA)
  x <- trimws(as.character(x))
  if (grepl("^[0-9]+s$", x, ignore.case = TRUE)) {
    return(as.numeric(gsub("s$", "", x, ignore.case = TRUE)))
  } else if (grepl("^[0-9]+$", x)) {
    return(as.numeric(x))
  } else {
    return(NA)
  }
})

# Prefer the bucketed estimate; fall back to the parsed raw number.
data$protest_size <- ifelse(!is.na(data$size_from_category),
                            data$size_from_category,
                            data$size_from_raw)

# Imputation: country-specific median, global median fallback, with indicator.
data$protest_size_missing <- ifelse(is.na(data$protest_size), 1, 0)

country_medians <- data %>%
  group_by(country) %>%
  summarise(med_size = median(protest_size, na.rm = TRUE), .groups = "drop")

data <- data %>%
  left_join(country_medians, by = "country") %>%
  mutate(
    protest_size = ifelse(is.na(protest_size), med_size, protest_size),
    protest_size = ifelse(is.na(protest_size),
                          median(protest_size, na.rm = TRUE),
                          protest_size)
  ) %>%
  select(-med_size)

data$log_protest_size <- log(data$protest_size + 1)
data$large_protest    <- ifelse(data$protest_size >= 1000, 1, 0)

# Provenance — how many observations came from each source.
cat("Category available:", sum(!is.na(data$size_from_category)), "\n")
cat("Raw filled gaps:",   sum(is.na(data$size_from_category) & !is.na(data$size_from_raw)), "\n")
cat("Imputed (missing indicator = 1):", sum(data$protest_size_missing), "\n")

data$size_from_category <- NULL
data$size_from_raw      <- NULL

## 1.4 Demand features --------------------------------------------------------
# The raw data lists up to 4 free-text demands per event. We collapse them into
# binary topic dummies via keyword matching, plus a count of total demands.

demand_cols <- c("protesterdemand1", "protesterdemand2",
                 "protesterdemand3", "protesterdemand4")

make_demand_dummy <- function(pattern) {
  apply(data[, demand_cols], 1, function(row) {
    as.integer(any(grepl(pattern, row, ignore.case = TRUE)))
  })
}

data$demand_labor      <- make_demand_dummy("labor|wage")
data$demand_land       <- make_demand_dummy("land|farm")
data$demand_police     <- make_demand_dummy("police|brutal")
data$demand_political  <- make_demand_dummy("political|process")
data$demand_prices     <- make_demand_dummy("price|tax")
data$demand_corruption <- make_demand_dummy("remov|corrupt")
data$demand_social     <- make_demand_dummy("social")

data$n_demands <- apply(data[, demand_cols], 1, function(row) {
  sum(row != "" & !is.na(row))
})

## 1.5 State-violence lags ----------------------------------------------------
# Past state repression is a strong theoretical predictor of escalation.
# We build (a) a 1-period lag and (b) the cumulative running mean of state
# violence, both within country and computed BEFORE the current event so they
# can serve as predictors without leakage.

response_cols <- c("stateresponse1", "stateresponse2", "stateresponse3",
                   "stateresponse4", "stateresponse5", "stateresponse6",
                   "stateresponse7")

# Temporary per-event flag: did the state respond with physical force?
data$violent_response <- apply(data[, response_cols], 1, function(row) {
  as.integer(any(grepl("kill|shoot|beat", row, ignore.case = TRUE)))
})

# Order chronologically within each country before computing lags.
data <- data[order(data$country, data$start_date), ]

data <- data %>%
  group_by(country) %>%
  mutate(
    lag1_state_violence = lag(violent_response, 1),
    cum_state_violence  = lag(cummean(violent_response), 1)
  ) %>%
  ungroup()

# First event per country has no prior history -> assume no prior repression.
data$lag1_state_violence[is.na(data$lag1_state_violence)] <- 0
data$cum_state_violence[is.na(data$cum_state_violence)]   <- 0

data$violent_response <- NULL  # only needed to compute the lags above

## 1.6 Finalize analytic sample ----------------------------------------------
# Year as a factor enables clean year fixed effects in the regressions.
data$year_factor <- factor(data$year)

cat("Final dataset:", nrow(data), "protest events,",
    length(unique(data$country)), "countries,",
    paste(range(data$year), collapse = "-"), "\n")

# Drop raw columns whose information is now captured by engineered features.
drop_cols <- c("id", "ccode", "protest",
               "startday", "startmonth", "startyear",
               "endday", "endmonth", "endyear",
               "start_date", "end_date",
               "location", "participants_category", "participants",
               "protesteridentity",
               "protesterdemand1", "protesterdemand2",
               "protesterdemand3", "protesterdemand4",
               "stateresponse1", "stateresponse2", "stateresponse3",
               "stateresponse4", "stateresponse5", "stateresponse6",
               "stateresponse7",
               "sources", "notes")
data <- data[, !(names(data) %in% drop_cols)]

View(data)


# 2. Exploratory Data Analysis ------------------------------------------------
# Descriptive view of the analytic sample before any modeling: a summary table,
# a missingness guard, outcome distribution by region / year / regressor, the
# size distribution that motivates the log transform, and a multicollinearity
# diagnostic (correlation matrix + VIFs).

## 2.1 Summary statistics -----------------------------------------------------
summary_vars <- c("protesterviolence",
                  "protest_size", "log_protest_size", "protest_size_missing",
                  "large_protest", "duration_days", "n_demands", "protestnumber",
                  "demand_labor", "demand_land", "demand_police",
                  "demand_political", "demand_prices", "demand_corruption",
                  "demand_social",
                  "lag1_state_violence", "cum_state_violence")

stargazer(as.data.frame(data[, summary_vars]),
          type = "text",
          title = "Descriptive statistics — analytic sample",
          digits = 3,
          summary.stat = c("n", "mean", "sd", "min", "p25", "median", "p75", "max"))

## 2.2 Missingness sanity check -----------------------------------------------
# All NAs are resolved by construction (size imputed, dummies are 0/1, lags
# filled with 0 for first event). Guards against silent regressions if the
# cleaning logic changes later.
stopifnot(sum(sapply(data, function(x) sum(is.na(x)))) == 0)
cat("Missingness check passed: 0 NAs across", ncol(data), "columns.\n")

## 2.3 Outcome distribution: violent vs peaceful ------------------------------
cat("\nOverall share of violent protests:",
    round(mean(data$protesterviolence, na.rm = TRUE), 3),
    "  (", sum(data$protesterviolence, na.rm = TRUE), "events)\n")

# Region-level shares + sample counts.
viol_by_region <- aggregate(protesterviolence ~ region, data = data,
                            FUN = function(x) c(n = length(x),
                                                share = mean(x, na.rm = TRUE)))
print(viol_by_region)

# Year-level trend.
viol_by_year <- aggregate(protesterviolence ~ year, data = data, FUN = mean)
plot(viol_by_year$year, viol_by_year$protesterviolence, type = "l",
     lwd = 2, col = "firebrick",
     xlab = "Year", ylab = "Share of violent protests",
     main = "Share of violent protests over time")
abline(h = mean(data$protesterviolence, na.rm = TRUE),
       col = "grey50", lty = 2)

## 2.4 Protest size distribution (raw vs log) ---------------------------------
# Motivation for the log transform applied during cleaning.
par(mfrow = c(1, 2))
hist(data$protest_size, breaks = 50,
     col = "steelblue", border = "white",
     main = "Protest Size (Raw)", xlab = "Participants")
abline(v = median(data$protest_size, na.rm = TRUE),
       col = "red", lwd = 2, lty = 2)
hist(data$log_protest_size, breaks = 30,
     col = "darkgreen", border = "white",
     main = "Protest Size (Log)", xlab = "log(Participants + 1)")
par(mfrow = c(1, 1))

## 2.5 Outcome vs key regressors ----------------------------------------------
# Quick visual evidence on whether candidate features actually shift P(violent).
par(mfrow = c(2, 2))

barplot(tapply(data$protesterviolence, data$region, mean),
        col = "steelblue", border = "white",
        main = "Violent share by region", ylab = "P(violent)")

barplot(tapply(data$protesterviolence, data$large_protest, mean),
        names.arg = c("Small (<1000)", "Large (>=1000)"),
        col = "darkorange", border = "white",
        main = "Violent share by protest size", ylab = "P(violent)")

barplot(tapply(data$protesterviolence, data$lag1_state_violence, mean),
        names.arg = c("No prior state violence", "Prior state violence"),
        col = "darkred", border = "white",
        main = "Violent share by lagged state violence", ylab = "P(violent)")

plot(density(data$log_protest_size[data$protesterviolence == 0]),
     col = "steelblue", lwd = 2,
     main = "log(size) density: violent vs peaceful",
     xlab = "log(protest_size + 1)")
lines(density(data$log_protest_size[data$protesterviolence == 1]),
      col = "firebrick", lwd = 2)
legend("topright", legend = c("Peaceful", "Violent"),
       col = c("steelblue", "firebrick"), lwd = 2, bty = "n")

par(mfrow = c(1, 1))

# Per-demand-dummy violent share (helps decide which topics matter).
demand_dummies <- c("demand_labor", "demand_land", "demand_police",
                    "demand_political", "demand_prices",
                    "demand_corruption", "demand_social")
demand_shares <- sapply(demand_dummies, function(v) {
  c(share_when_1 = mean(data$protesterviolence[data[[v]] == 1], na.rm = TRUE),
    share_when_0 = mean(data$protesterviolence[data[[v]] == 0], na.rm = TRUE),
    n_when_1     = sum(data[[v]] == 1, na.rm = TRUE))
})
print(round(t(demand_shares), 3))

## 2.6 Correlation matrix & VIF -----------------------------------------------
# Pre-modeling multicollinearity diagnostic.
num_vars <- c("log_protest_size", "duration_days", "n_demands", "protestnumber",
              "large_protest", "protest_size_missing",
              "demand_labor", "demand_land", "demand_police",
              "demand_political", "demand_prices", "demand_corruption",
              "demand_social",
              "lag1_state_violence", "cum_state_violence")

cor_mat <- cor(data[, num_vars], use = "pairwise.complete.obs")
corrplot(cor_mat, method = "color", type = "upper",
         tl.cex = 0.7, tl.col = "black",
         addCoef.col = "black", number.cex = 0.55,
         title = "Correlation matrix of candidate regressors",
         mar = c(0, 0, 1.5, 0))

vif_lpm <- lm(protesterviolence ~ . - country - year - year_factor,
              data = data[, c("protesterviolence", num_vars,
                              "region", "country", "year", "year_factor")])
vif_vals <- car::vif(vif_lpm)
cat("\nVariance Inflation Factors (GVIF^(1/(2*Df)) for factors):\n")
print(round(vif_vals, 3))


# 3. General-to-Specific Feature Selection ------------------------------------
# Following the LSE (Hendry) G2S approach demonstrated in Lab 08.
# Protocol (mirrors AE_Lab_08_Gets_Better_codes.R):
#   Start from the GUM containing every candidate regressor.
#   At each step:
#     1. Read cluster-robust p-values (coeftest + vcovCL by country).
#     2. Identify the scalar regressor with the highest p-value.
#     3. Test that ALL variables dropped so far (including this new one) are
#        jointly zero in the GUM via an LR test (lrtest(gum, restricted)).
#        - If p > 0.05 for the LR test (cannot reject H0) -> drop the variable.
#        - If p <= 0.05 (would significantly worsen the GUM) -> stop.
#   Stop when every remaining scalar coefficient is significant at alpha = 0.05.
#   The surviving model is saved as `model_final`.
#
# NOTE: year_factor (multi-level FE) and region (factor) are treated as blocks:
# they are never dropped individually — only as a group via a Wald/LR block test
# at the end if all their constituent dummies are jointly insignificant.

## 3.1 General Unrestricted Model (GUM) ---------------------------------------
# All substantive candidates plus region and year fixed effects.
gum <- glm(protesterviolence ~
             log_protest_size +
             large_protest +
             protest_size_missing +
             duration_days +
             demand_political +
             demand_labor +
             demand_land +
             demand_police +
             demand_prices +
             demand_corruption +
             demand_social +
             n_demands +
             protestnumber +
             lag1_state_violence +
             cum_state_violence +
             region +
             year_factor,
           data = data,
           family = binomial(link = "logit"))

# Cluster-robust SEs (clustered by country, same as Lab 08 / standard in panel
# and cross-country data to account for within-country error correlation).
ct_gum <- coeftest(gum, vcov = vcovCL(gum, cluster = data$country))
print(ct_gum)

## 3.2 Step 1: candidate drop — demand_social (p = 0.58) ----------------------
# Most insignificant scalar regressor in the GUM. Test whether {demand_social=0}
# is jointly acceptable in the GUM via LR test (Lab 08 pattern: always test
# the cumulative dropped set against the GUM, not against the current model).
step1 <- glm(protesterviolence ~
               log_protest_size +
               large_protest +
               protest_size_missing +
               duration_days +
               demand_political +
               demand_labor +
               demand_land +
               demand_police +
               demand_prices +
               demand_corruption +
               n_demands +
               protestnumber +
               lag1_state_violence +
               cum_state_violence +
               region +
               year_factor,
             data = data,
             family = binomial(link = "logit"))

# LR test: H0: beta_demand_social = 0 in the GUM
# p = 0.4894 -> cannot reject H0 -> demand_social can be dropped.
lrtest(gum, step1)

ct_step1 <- coeftest(step1, vcov = vcovCL(step1, cluster = data$country))
print(ct_step1)

## 3.3 Step 2: candidate drop — protestnumber (p = 0.20 in Step 1) ------------
# Next most insignificant scalar. Test whether {demand_social=0, protestnumber=0}
# are jointly zero in the GUM.
step2_candidate <- glm(protesterviolence ~
                         log_protest_size +
                         large_protest +
                         protest_size_missing +
                         duration_days +
                         demand_political +
                         demand_labor +
                         demand_land +
                         demand_police +
                         demand_prices +
                         demand_corruption +
                         n_demands +
                         lag1_state_violence +
                         cum_state_violence +
                         region +
                         year_factor,
                       data = data,
                       family = binomial(link = "logit"))

# LR test: H0: beta_demand_social = beta_protestnumber = 0 in the GUM
# p = 0.0063 -> REJECT H0 -> cannot drop protestnumber -> STOP elimination.
lrtest(gum, step2_candidate)

## 3.4 Final model ------------------------------------------------------------
# Elimination stops after Step 1. The final model is the GUM minus demand_social.
# `protestnumber` remains despite its p = 0.20 because the joint LR test
# confirms it contributes significantly when tested against the full GUM.
model_final <- step1

ct_final <- coeftest(model_final, vcov = vcovCL(model_final, cluster = data$country))
cat("\n=== Final model: cluster-robust standard errors ===\n")
print(ct_final)

## 3.5 Block tests: year_factor and region ------------------------------------
# Individual levels of a factor block are never dropped one at a time in G2S.
# Instead, we test each block as a whole: can ALL year dummies simultaneously
# be set to zero without significantly worsening model_final?

# year_factor block: model_final vs model_final without year_factor
model_no_year <- glm(protesterviolence ~
                       log_protest_size +
                       large_protest +
                       protest_size_missing +
                       duration_days +
                       demand_political +
                       demand_labor +
                       demand_land +
                       demand_police +
                       demand_prices +
                       demand_corruption +
                       n_demands +
                       protestnumber +
                       lag1_state_violence +
                       cum_state_violence +
                       region,
                     data = data,
                     family = binomial(link = "logit"))

# LR test: H0: all year dummies = 0 jointly
# Chi2(30) = 103, p = 6.2e-10 -> REJECT H0 -> year FE are jointly significant,
# keep year_factor in the final model.
lrtest(model_final, model_no_year)

# region block: model_final vs model_final without region
model_no_region <- glm(protesterviolence ~
                         log_protest_size +
                         large_protest +
                         protest_size_missing +
                         duration_days +
                         demand_political +
                         demand_labor +
                         demand_land +
                         demand_police +
                         demand_prices +
                         demand_corruption +
                         n_demands +
                         protestnumber +
                         lag1_state_violence +
                         cum_state_violence +
                         year_factor,
                       data = data,
                       family = binomial(link = "logit"))

# LR test: H0: all region dummies = 0 jointly
# Chi2(7) = 129, p < 2.2e-16 -> REJECT H0 -> region FE are jointly significant,
# keep region in the final model.
lrtest(model_final, model_no_region)

## 3.6 GUM vs final model comparison table ------------------------------------
# Displays both models side by side for the report; shows what was dropped and
# what changed in significance.
stargazer(gum, model_final,
          type = "text",
          title = "G2S selection: GUM vs final model",
          column.labels = c("GUM", "Final"),
          omit = "year_factor",
          omit.labels = "Year FE",
          add.lines = list(c("Year FE", "Yes", "Yes")),
          digits = 3)

## 3.7 Interaction term extension --------------------------------------------
# Theoretical motivation: prior state violence may have a stronger effect on
# LARGE protests, because the combination of crowd dynamics and state
# escalation could be multiplicative (escalation amplifies with crowd size).
# We test whether adding large_protest * lag1_state_violence improves fit.
model_int <- glm(update(model_final$formula,
                        . ~ . + large_protest:lag1_state_violence),
                 data = data, family = binomial(link = "logit"))

ct_int <- coeftest(model_int, vcov = vcovCL(model_int, cluster = data$country))
cat("\nInteraction term (large_protest x lag1_state_violence):\n")
print(ct_int[grep("large_protest|lag1_state", rownames(ct_int)), ])

# LR test: H0: interaction coefficient = 0
# Chi2(1) = 0.07, p = 0.79 -> cannot reject H0 -> interaction does NOT improve
# fit. The marginal effect of prior state violence is statistically the same
# for large and small protests. We retain model_final without the interaction
# but report this test as required diagnostic content.
lrtest(model_final, model_int)


# 4. Modeling and Interpretation ----------------------------------------------
# All tests operate on `model_final` (logit, GUM minus demand_social, from Part 3).
# Pattern follows Lab 03: LR vs null -> Wald block -> quasi-separation check ->
# linktest -> GoF (HL / OR / Stukel) -> pseudo-R2 -> ROC -> confusion matrix ->
# marginal effects -> LPM benchmark -> probit robustness -> stargazer table.

## 4.1 Joint significance: LR test vs null model -----------------------------
# H0: all slope coefficients = 0 simultaneously.
# Chi2(51) = 1162, p < 2.2e-16 -> model is globally highly significant.
null_model <- glm(protesterviolence ~ 1, data = data,
                  family = binomial(link = "logit"))
lrtest(null_model, model_final)

## 4.2 Wald test: all demand dummies jointly = 0 -----------------------------
# Demand-topic dummies were added for theoretical reasons; test them as a block.
# Chi2(6) = 150, p ≈ 0 -> demand type matters jointly, keep all six dummies.
demand_idx <- grep("^demand_", names(coef(model_final)))
H_demand   <- diag(length(coef(model_final)))[demand_idx, ]
wald.test(b = coef(model_final), Sigma = vcov(model_final), L = H_demand)

## 4.3 Quasi-separation check -------------------------------------------------
# Fitted probs outside (0.001, 0.999) would indicate near-perfect separation
# and require Firth's bias-reduced logit (logistf).
# Min = 0.047, Max = 0.990 -> no separation; standard logit is appropriate.
fp <- model_final$fitted.values
cat("Fitted prob range: [", round(min(fp), 4), ",", round(max(fp), 4), "]\n")
cat("Extreme obs (p < 0.001 or > 0.999):", sum(fp < 0.001 | fp > 0.999), "\n")

## 4.4 Linktest (specification check) ----------------------------------------
# yhat should be significant, yhat2 should NOT be (Lab 03, linktest.R).
# Result: yhat2 p = 4.65e-07 -> significant -> suggests possible misspecification
# (likely missing nonlinear or interaction terms; noted as a caveat).
source("class_materials/Lab 03 (2026-03-06)-20260526/linktest.R")
linktest_result <- linktest(model_final)

## 4.5 Goodness-of-fit tests --------------------------------------------------
# Three complementary GoF tests from AllGOFTests.R (Lab 03).
source("class_materials/Lab 03 (2026-03-06)-20260526/AllGOFTests.R")

# Hosmer-Lemeshow (g = 10 groups): X2 = 16.3, p = 0.039
# -> marginal rejection at 5%; model calibration could be improved.
HLTest(model_final, g = 10)

# Osius-Rojek: z = 1.69, p = 0.091
# -> cannot reject H0 at 5%; passes this test.
o.r.test(model_final)

# Stukel: stat = 39.1, p = 3.2e-09
# -> strongly rejects; confirms the linktest signal that the logit link
# function may not capture the full nonlinearity in the tails.
stukel.test(model_final)

## 4.6 Pseudo-R2 statistics ---------------------------------------------------
# None of these is R2 in the OLS sense; McFadden > 0.2 is conventionally
# considered a good fit for binary logit; 0.066 is modest but typical for
# aggregate protest data with significant unobserved heterogeneity.
ll_full <- as.numeric(logLik(model_final))
ll_null <- as.numeric(logLik(null_model))
n       <- nrow(data)

R2_mcfadden  <- 1 - ll_full / ll_null
R2_coxsnell  <- 1 - exp((2 / n) * (ll_null - ll_full))
R2_nagelkerke <- R2_coxsnell / (1 - exp((2 / n) * ll_null))

cat(sprintf("McFadden R2:   %.4f\n", R2_mcfadden))    # 0.0660
cat(sprintf("Cox-Snell R2:  %.4f\n", R2_coxsnell))    # 0.0734
cat(sprintf("Nagelkerke R2: %.4f\n", R2_nagelkerke))  # 0.1072

# McKelvey-Zavoina R2: based on the latent-variable interpretation of logit.
# Formula: var(linear index) / (var(linear index) + pi^2/3) for logit.
# Best comparable to OLS R2 because it operates on the latent y*.
# Lab 03 (AE_Lab_03_Logit.R, line ~91) uses BaylorEdPsych::PseudoR2() which
# returns this same statistic; we compute it manually to avoid the archived
# CRAN dependency.
yhat_idx <- predict(model_final, type = "link")
R2_MZ    <- var(yhat_idx) / (var(yhat_idx) + pi^2 / 3)
cat(sprintf("McKelvey-Zavoina R2: %.4f\n", R2_MZ))    # 0.1113

# Count R2: classification accuracy at the 0.5 cutoff.
# Naive measure; can look high even for a poor model if the outcome is rare.
# Lab 03 obtains this from BaylorEdPsych::PseudoR2() as well.
pred05    <- as.integer(model_final$fitted.values >= 0.5)
correct   <- sum(pred05 == data$protesterviolence)
R2_count  <- correct / n
cat(sprintf("Count R2:            %.4f  (%d / %d correct)\n",
            R2_count, correct, n))                     # 0.7417

# Adjusted count R2: improvement over the majority-class baseline (Long 1997).
# Negative or near-zero -> model barely beats predicting the modal outcome.
# Lab 03 reference: BaylorEdPsych::PseudoR2() "Count Adj" entry.
n_majority   <- max(table(data$protesterviolence))
R2_count_adj <- (correct - n_majority) / (n - n_majority)
cat(sprintf("Adjusted count R2:   %.4f  (majority baseline = %d)\n",
            R2_count_adj, n_majority))                 # 0.0245

# Verdict: the model's added value over "always predict peaceful" is real
# but small in absolute terms (only ~2.5pp gain in classification),
# matching the AUC=0.677 evidence of modest discrimination.

## 4.7 ROC curve & AUC --------------------------------------------------------
# AUC = 0.677: model discriminates meaningfully above chance (0.5) but has
# room for improvement — consistent with the unobserved-heterogeneity caveat.
g_roc <- roc(data$protesterviolence, model_final$fitted.values, quiet = TRUE)
plot(g_roc, col = "firebrick", lwd = 2,
     main = paste("ROC curve — Final logit model  (AUC =",
                  round(auc(g_roc), 3), ")"))
abline(a = 0, b = 1, col = "grey60", lty = 2)

## 4.8 Confusion matrix at optimal threshold ----------------------------------
# Youden-optimal threshold = 0.264 (balances sensitivity and specificity).
# Accuracy 0.646, Sensitivity 0.614, Specificity 0.658.
best_coords <- coords(g_roc, "best",
                      ret = c("threshold", "sensitivity", "specificity", "accuracy"))
print(round(best_coords, 4))

thresh     <- as.numeric(best_coords["threshold"])
pred_class <- ifelse(model_final$fitted.values >= thresh, 1, 0)
cm_tbl     <- table(Predicted = pred_class, Actual = data$protesterviolence)
print(cm_tbl)
cat("Accuracy:   ", round(mean(pred_class == data$protesterviolence), 4), "\n")
cat("Sensitivity:", round(cm_tbl[2, 2] / sum(cm_tbl[, 2]), 4), "\n")
cat("Specificity:", round(cm_tbl[1, 1] / sum(cm_tbl[, 1]), 4), "\n")

## 4.9 Marginal effects at the mean (cluster-robust) --------------------------
# Average partial effects evaluated at the sample mean of all regressors.
# year_factor omitted from the mfx formula (many dummies; region is kept).
# Interpretation examples:
#   demand_police: +26pp — protests demanding police accountability are the
#     most likely to turn violent (consistent with confrontational framing).
#   cum_state_violence: +35pp — countries with a higher historical rate of
#     state repression see dramatically more violent protests (escalation trap).
#   log_protest_size: -2pp per log-unit — larger crowds marginally reduce
#     per-event violence risk (diffusion / policing effect).
fml_no_yf <- protesterviolence ~
  log_protest_size + large_protest + protest_size_missing + duration_days +
  demand_political + demand_labor + demand_land + demand_police +
  demand_prices + demand_corruption +
  n_demands + protestnumber + lag1_state_violence + cum_state_violence + region

logitmfx(fml_no_yf, data = data, atmean = TRUE,
         robust = TRUE, clustervar1 = "country")

## 4.10 LPM benchmark ---------------------------------------------------------
# The linear probability model is estimated for comparison.
# BP test: p < 2.2e-16 -> strong heteroscedasticity (expected for LPM with a
# binary DV), which is why we prefer the logit.
# Coefficient signs and significance are directionally consistent with logit.
lpm <- lm(update(fml_no_yf, . ~ . + year_factor), data = data)
bptest(lpm)  # Breusch-Pagan: BP = 866, p < 2.2e-16

ct_lpm <- coeftest(lpm, vcov = vcovHC(lpm, type = "HC"))
scalar_rows <- !grepl("^year_factor|^region", rownames(ct_lpm))
cat("LPM HC-robust scalar coefficients:\n")
print(round(ct_lpm[scalar_rows, ], 4))

## 4.11 Probit robustness check -----------------------------------------------
# Estimated with the same specification as model_final.
# All signs identical to logit; same variables significant at 5%.
# Confirms results are not an artefact of the logistic functional form.
model_probit <- glm(model_final$formula, data = data,
                    family = binomial(link = "probit"))
ct_probit <- coeftest(model_probit,
                      vcov = vcovCL(model_probit, cluster = data$country))
scalar_p <- !grepl("^year_factor|^region", rownames(ct_probit))
cat("Probit cluster-robust scalar coefficients:\n")
print(round(ct_probit[scalar_p, ], 4))

## 4.12 Three-model comparison table (logit / probit / LPM) ------------------
stargazer(model_final, model_probit, lpm,
          type = "text",
          title = "Final model comparison: logit, probit, LPM",
          column.labels = c("Logit", "Probit", "LPM"),
          omit = "year_factor",
          omit.labels = "Year FE",
          add.lines = list(c("Year FE", "Yes", "Yes", "Yes"),
                           c("SE type", "Cluster", "Cluster", "HC")),
          digits = 3)


# 5. Conclusion ---------------------------------------------------------------
cat("
================================================================================
CONCLUSION
================================================================================

This study examines what makes mass demonstrations turn violent using 15,239
protest events from the Mass Mobilization dataset (1990-2020, 166 countries).
A binary logit model was selected via LSE general-to-specific elimination and
estimated with country-clustered standard errors throughout.

KEY FINDINGS
------------
1. Cumulative state repression (+35pp) is the dominant predictor of protest
   violence. Countries with a long history of violent state responses trap
   subsequent protests in an escalation cycle, regardless of the immediate
   cause of mobilisation.

2. Demand type is the second major driver. Protests demanding police
   accountability (+26pp), price/tax relief (+24pp), and land rights (+15pp)
   are the most likely to turn violent, reflecting the inherently confrontational
   framing of these grievances. Protests combining more demand types are
   paradoxically less violent (-8pp per additional demand), possibly because
   broader coalitions are harder to provoke.

3. Protest size has a small negative effect (-2pp per log-unit), suggesting
   that larger, more visible mobilisations face heavier policing or benefit
   from crowd diffusion effects that dampen individual-level violence.

4. Short-run state violence (1-period lag, +6pp) reinforces the long-run
   repression effect, confirming that escalation is path-dependent.

5. Regional effects are substantial: protests in Europe (-10pp) and MENA
   (-9pp) are less likely to be violent than the African baseline, conditional
   on all other covariates, pointing to unobserved institutional differences.

MODEL FIT AND CAVEATS
---------------------
Fit is modest (McFadden R2 = 0.066, AUC = 0.677), as expected for
aggregate event-level data where unobserved factors (police tactics on the
day, protest leadership, media coverage) play a large role. The Stukel test
and linktest both flag potential nonlinearity in the tails of the index
function, suggesting that interaction terms or higher-order size effects
could improve calibration. All qualitative conclusions are robust to
replacing the logit with a probit or a heteroscedasticity-corrected LPM.

POLICY IMPLICATION
------------------
The single most actionable lever is the state's own behaviour: reducing
habitual repressive responses lowers the baseline probability of violence
at future protests. Demand-side interventions (addressing police brutality
complaints, price pressures) also reduce violence risk by defusing the
most confrontational mobilisation types.
================================================================================
")

