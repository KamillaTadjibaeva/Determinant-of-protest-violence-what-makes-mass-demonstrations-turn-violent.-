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


