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


