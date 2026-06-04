# 0. Setup
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


# 1. Data Loading & Cleaning
# We read the raw Mass Mobilization file, keep only protest events, and build the
# analytic features we use in the model (size, duration, demand dummies, lags).

# Load & filter
data <- read.csv("dataverse_files/mmALL_073120_csv.csv", stringsAsFactors = FALSE)
str(data)

# We drop non-protest rows since protesterviolence is undefined for them.
data <- data[data$protest == 1, ]

data$country <- factor(data$country)
data$region  <- factor(data$region)
# Region codes (left as numeric labels):
# 1 South America, 2 Central America, 3 North America, 4 Europe,
# 5 Asia, 6 MENA, 7 Africa

# Dates & duration
# We build proper Date objects and compute event length in days.
data$start_date <- as.Date(paste(data$startyear, data$startmonth, data$startday, sep = "-"),
                           format = "%Y-%m-%d")
data$end_date   <- as.Date(paste(data$endyear,   data$endmonth,   data$endday,   sep = "-"),
                           format = "%Y-%m-%d")

data$duration_days <- as.numeric(data$end_date - data$start_date)
data$duration_days[data$duration_days < 0] <- 0  # guard against bad rows

# Protest size
# Two source columns exist: a bucketed `participants_category` and a free-text
# `participants` field. We harmonise both into a single numeric `protest_size`,
# log-transform it (heavy right skew), and impute missing values by country
# median while flagging the imputation with a dummy.

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

# Demand features
# The raw data lists up to 4 free-text demands per event. We collapse them into
# binary topic dummies via keyword matching, plus we count total demands.

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

# State-violence lags
# Past state repression is a strong theoretical predictor of escalation.
# We build (a) a 1-period lag and (b) the cumulative running mean of state
# violence, both within country and computed BEFORE the current event so we
# avoid information leakage.

response_cols <- c("stateresponse1", "stateresponse2", "stateresponse3",
                   "stateresponse4", "stateresponse5", "stateresponse6",
                   "stateresponse7")

# We create a temporary per-event flag: did the state respond with physical force?
data$violent_response <- apply(data[, response_cols], 1, function(row) {
  as.integer(any(grepl("kill|shoot|beat", row, ignore.case = TRUE)))
})

# We order chronologically within each country before computing lags.
data <- data[order(data$country, data$start_date), ]

data <- data %>%
  group_by(country) %>%
  mutate(
    lag1_state_violence = lag(violent_response, 1),
    cum_state_violence  = lag(cummean(violent_response), 1)
  ) %>%
  ungroup()

# First event per country has no prior history -> we assume no prior repression.
data$lag1_state_violence[is.na(data$lag1_state_violence)] <- 0
data$cum_state_violence[is.na(data$cum_state_violence)]   <- 0

data$violent_response <- NULL  # only needed to compute the lags above

# Finalize analytic sample
# Year as a factor enables clean year fixed effects in the regressions.
data$year_factor <- factor(data$year)

# Region dummies (baseline = Africa)
data$region_southam  <- as.integer(data$region == "South America")
data$region_centam   <- as.integer(data$region == "Central America")
data$region_northam  <- as.integer(data$region == "North America")
data$region_europe   <- as.integer(data$region == "Europe")
data$region_asia     <- as.integer(data$region == "Asia")
data$region_mena     <- as.integer(data$region == "MENA")
data$region_oceania  <- as.integer(data$region == "Oceania")

# Year dummies (baseline = 1990)
for (yr in 1991:2020) {
  data[[paste0("year_", yr)]] <- as.integer(data$year == yr)
}

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


# 2. Exploratory Data Analysis
# We present a descriptive view of the analytic sample before modeling: summary
# table, missingness check, outcome distribution by region/year/regressor,
# size distribution motivating the log transform, and multicollinearity
# diagnostics (correlation matrix + VIFs).

# Summary statistics
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

# Missingness sanity check
# All NAs should be resolved by construction (size imputed, dummies are 0/1,
# lags filled with 0 for first event). We verify this here.
stopifnot(sum(sapply(data, function(x) sum(is.na(x)))) == 0)
cat("Missingness check passed: 0 NAs across", ncol(data), "columns.\n")

# Outcome distribution: violent vs peaceful
cat("\nOverall share of violent protests:",
    round(mean(data$protesterviolence, na.rm = TRUE), 3),
    "  (", sum(data$protesterviolence, na.rm = TRUE), "events)\n")

# Region-level shares + sample counts.
viol_by_region <- aggregate(protesterviolence ~ region, data = data,
                            FUN = function(x) c(n = length(x),
                                                share = mean(x, na.rm = TRUE)))
print(viol_by_region)

# Country-level shares + sample counts, sorted highest to lowest.
viol_by_country <- aggregate(protesterviolence ~ country, data = data,
                             FUN = function(x) c(n = length(x),
                                                 share = mean(x, na.rm = TRUE)))
viol_by_country <- viol_by_country[order(-viol_by_country$protesterviolence[, "share"]), ]
print(viol_by_country)

# Year-level trend.
viol_by_year <- aggregate(protesterviolence ~ year, data = data, FUN = mean)
plot(viol_by_year$year, viol_by_year$protesterviolence, type = "l",
     lwd = 2, col = "firebrick",
     xlab = "Year", ylab = "Share of violent protests",
     main = "Share of violent protests over time")
abline(h = mean(data$protesterviolence, na.rm = TRUE),
       col = "grey50", lty = 2)

# Protest size distribution (raw vs log)
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

# Outcome vs key regressors
# We check visually whether candidate features actually shift P(violent).
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
  c(pct_violent_with_demand    = mean(data$protesterviolence[data[[v]] == 1], na.rm = TRUE),
    pct_violent_without_demand = mean(data$protesterviolence[data[[v]] == 0], na.rm = TRUE),
    n_protests_with_demand     = sum(data[[v]] == 1, na.rm = TRUE))
})
print(round(t(demand_shares), 3))

# Correlation matrix & VIF
# We check for multicollinearity before modeling.
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


# 3. General-to-Specific Feature Selection

# LOGIT VS PROBIT via AIC/BIC
gum_probit <- glm(gum$formula, data = data, family = binomial(link = "probit"))


#GENERAL TO SPECIFIC
# General Unrestricted Model (GUM)
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
             region_asia + region_centam + region_europe + region_mena +
             region_northam + region_oceania + region_southam +
             year_1991 + year_1992 + year_1993 + year_1994 + year_1995 +
             year_1996 + year_1997 + year_1998 + year_1999 + year_2000 +
             year_2001 + year_2002 + year_2003 + year_2004 + year_2005 +
             year_2006 + year_2007 + year_2008 + year_2009 + year_2010 +
             year_2011 + year_2012 + year_2013 + year_2014 + year_2015 +
             year_2016 + year_2017 + year_2018 + year_2019 + year_2020,
           data = data,
           family = binomial(link = "logit"))
summary(gum)


#First testing if insignificant variables are jointly insignificant

#model without ALL insignificant variables at once: 
#demand_social, region_oceania, year_1993, year_1994, year_1995,
#year_1998, year_1999, year_2001, year_2006, year_2009, 
#year_2014, year_2019
gum1 <- glm(protesterviolence ~
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
             region_asia + region_centam + region_europe + region_mena +
             region_northam + region_southam +
             year_1991 + year_1992 +
             year_1996 + year_1997 + year_2000 + year_2002 + year_2003 + 
              year_2004 + year_2005 +
              + year_2007 + year_2008 + year_2010 +
             year_2011 + year_2012 + year_2013 + year_2015 +
             year_2016 + year_2017 + year_2018 + year_2020,
           data = data,
           family = binomial(link = "logit"))
lrtest(gum, gum1)
#Using lrtest for logit
# H0: they are jointly insignificant
#Fail to reject H0, so they are jointly insignificant, we can drop them
summary(gum1)

#2nd model without ALL insignificant variables at once: 
#year_1991, year_1992, year_2020
gum2 <- glm(protesterviolence ~
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
              region_asia + region_centam + region_europe + region_mena +
              region_northam + region_southam +
              year_1996 + year_1997 + year_2000 + year_2002 + year_2003 + 
              year_2004 + year_2005 +
              + year_2007 + year_2008 + year_2010 +
              year_2011 + year_2012 + year_2013 + year_2015 +
              year_2016 + year_2017 + year_2018,
            data = data,
            family = binomial(link = "logit"))

lrtest(gum, gum2)
#H0: dropped variables are jointly insignificant
#Fail to reject, so dropped variables are jointly insignificant
summary(gum2)

#3rd model without ALL insignificant variables at once: 
#year_1996, year_2008
gum3 <- glm(protesterviolence ~
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
              region_asia + region_centam + region_europe + region_mena +
              region_northam + region_southam +
              year_1997 + year_2000 + year_2002 + year_2003 + 
              year_2004 + year_2005 +
              + year_2007 + year_2010 +
              year_2011 + year_2012 + year_2013 + year_2015 +
              year_2016 + year_2017 + year_2018,
            data = data,
            family = binomial(link = "logit"))

lrtest(gum, gum3)
#H0: dropped variables are jointly insignificant
#Fail to reject, so dropped variables are jointly insignificant
summary(gum3)

#4th model without ALL insignificant variables at once: 
#year_2007, year_2011, year_2013, year_2018
gum4 <- glm(protesterviolence ~
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
              region_asia + region_centam + region_europe + region_mena +
              region_northam + region_southam +
              year_1997 + year_2000 + year_2002 + year_2003 + 
              year_2004 + year_2005 +
              + year_2010 +
              year_2012 + year_2015 +
              year_2016 + year_2017,
            data = data,
            family = binomial(link = "logit"))

lrtest(gum, gum4)
#H0: dropped variables are jointly insignificant
#Reject, so dropped variables are jointly SIGNIFICANT, so we have to remove them one by one
summary (gum3)

#Now we will remove the insignificant variables one by one
#starting with the most insignificant year_2013

gum5 <- glm(protesterviolence ~
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
              region_asia + region_centam + region_europe + region_mena +
              region_northam + region_southam +
              year_1997 + year_2000 + year_2002 + year_2003 + 
              year_2004 + year_2005 +
              + year_2007 + year_2010 +
              year_2011 + year_2012 + year_2015 +
              year_2016 + year_2017 + year_2018,
            data = data,
            family = binomial(link = "logit"))

summary(gum5)
#The new model still has some insignificant variables
#So we drop the most insignificant out of those which is year_2011

linearHypothesis(gum, c(
  "demand_social=0",
  "region_oceania=0",
  "year_1993=0", "year_1994=0", "year_1995=0",
  "year_1998=0", "year_1999=0", "year_2001=0",
  "year_2006=0", "year_2009=0", "year_2014=0", "year_2019=0",
  "year_1991=0", "year_1992=0", "year_2020=0",
  "year_1996=0", "year_2008=0",
  "year_2013=0", "year_2011=0"
), test = "Chisq")

#H0: they are indeed all equal to 0
#Fail to reject , so they are indeed all equal to 0
#So we can drop year_2011 form the model

gum6 <- glm(protesterviolence ~
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
              region_asia + region_centam + region_europe + region_mena +
              region_northam + region_southam +
              year_1997 + year_2000 + year_2002 + year_2003 + 
              year_2004 + year_2005 +
              + year_2007 + year_2010 +
              year_2012 + year_2015 +
              year_2016 + year_2017 + year_2018,
            data = data,
            family = binomial(link = "logit"))

summary(gum6)
#The new model still has some insignificant variables
#So we drop the most insignificant out of those which is year_2007

linearHypothesis(gum, c(
  "demand_social=0",
  "region_oceania=0",
  "year_1993=0", "year_1994=0", "year_1995=0",
  "year_1998=0", "year_1999=0", "year_2001=0",
  "year_2006=0", "year_2009=0", "year_2014=0", "year_2019=0",
  "year_1991=0", "year_1992=0", "year_2020=0",
  "year_1996=0", "year_2008=0",
  "year_2013=0", "year_2011=0", "year_2007=0"
), test = "Chisq")

#H0: they are indeed all equal to 0
#Fail to reject , so they are indeed all equal to 0
#So we can drop year_2007 form the model

gum7 <- glm(protesterviolence ~
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
              region_asia + region_centam + region_europe + region_mena +
              region_northam + region_southam +
              year_1997 + year_2000 + year_2002 + year_2003 + 
              year_2004 + year_2005 +
              + year_2010 +
              year_2012 + year_2015 +
              year_2016 + year_2017 + year_2018,
            data = data,
            family = binomial(link = "logit"))

summary(gum7)
#The new model still has some insignificant variables
#So we drop the most insignificant out of those which is year_2018

linearHypothesis(gum, c(
  "demand_social=0",
  "region_oceania=0",
  "year_1993=0", "year_1994=0", "year_1995=0",
  "year_1998=0", "year_1999=0", "year_2001=0",
  "year_2006=0", "year_2009=0", "year_2014=0", "year_2019=0",
  "year_1991=0", "year_1992=0", "year_2020=0",
  "year_1996=0", "year_2008=0",
  "year_2013=0", "year_2011=0", "year_2007=0", "year_2018=0"
), test = "Chisq")

#H0: they are indeed all equal to 0
#Reject , so they are NOT all equal to 0
#So we cannot drop year_2018

#Now we have to stop the process, this is the best model we can get. 
model_final <- gum7



# GUM vs final model comparison table
stargazer(gum, model_final,
          type = "text",
          title = "G2S selection: GUM vs final model",
          column.labels = c("GUM", "Final"),
          omit = "year_factor",
          omit.labels = "Year FE",
          add.lines = list(c("Year FE", "Yes", "Yes")),
          digits = 3)


# Models Diagnostic tests

# Linktest (specification check)
# We run linktest: yhat should be significant, yhat2 should NOT be.
# Result: yhat2 is significant -> suggests possible misspecification,
# likely missing nonlinear or interaction terms.
source("class_materials/Lab 03 (2026-03-06)-20260526/linktest.R")
linktest_result <- linktest(model_final)

# Linktest rejects (yhat2 significant) -> misspecification detected.
# We try adding theoretically motivated interaction variables.

# Interaction 1: large_protest × lag1_state_violence
model_int1 <- glm(update(model_final$formula,
                         . ~ . + large_protest:lag1_state_violence),
                  data = data, family = binomial(link = "logit"))
lrtest(model_final, model_int1)
# p = 0.81 -> NOT significant. Large protests don't react differently to
# prior state violence. Discard this interaction.

# Interaction 2: cum_state_violence × demand_police
# We expect that police-brutality protests in chronically repressive countries
# escalate disproportionately (compounding grievance + opportunity).
model_int2 <- glm(update(model_final$formula,
                         . ~ . + cum_state_violence:demand_police),
                  data = data, family = binomial(link = "logit"))
summary(model_int2)
lrtest(model_final, model_int2)
#interaction is significant, so we will update the final model with the new interaction
model_final <- model_int2


# Interaction 3: log_protest_size × duration_days
# We expect that longer protests with more participants have compounding
# escalation dynamics (sustained crowds build frustration over time).
model_int3 <- glm(update(model_final$formula,
                         . ~ . + log_protest_size:duration_days),
                  data = data, family = binomial(link = "logit"))
summary(model_int3)
lrtest(model_final, model_int3)
#interaction is significant, so we will update the final model with the new interaction
model_final <- model_int3
summary(model_final)

# Re-run linktest on updated model
linktest(model_final)

#adding interaction variables did not help to make the yhat2 insignificant, 
# but it did help to make it LESS significant

source("class_materials/Lab 03 (2026-03-06)-20260526/AllGOFTests.R")

# Hosmer-Lemeshow (g = 10 groups): X2 = 16.3, p = 0.039
#H0: model is ok
HLTest(model_final, g = 10)
#p-value = 0.01235 , Reject, so model is not ok

# Osius-Rojek: z = 1.69, p = 0.091
#H0: Test is ok
o.r.test(model_final)
#p-value =  0.8285437 , so Fail to reject, so the model is ok

# Stukel test
#H0: model is ok
stukel.test(model_final)
#p-value =  7.23038e-07 , so reject, so model is not ok


# Pseudo-R2 statistics
library(BaylorEdPsych)
PseudoR2(model_final)
# Count:  74% of all observations predicted correctly
# McKelvey.Zavoina:  if hidden variable was observed, then it would explain 11.4% of variation
# Adj.Count 3.2% of all observations correctly ONLY due to variation of X variable

# HYPOTHESIS TESTING

# LR test vs null model
# H0: all slope coefficients = 0 simultaneously.
null_model <- glm(protesterviolence ~ 1, data = data,
                  family = binomial(link = "logit"))
lrtest(null_model, model_final)
# We reject H0 -> our model coefficients are jointly significant, so the model
# has predictive power beyond the intercept-only (naive) model.

# Wald test: all demand dummies jointly = 0
# H0: all demand coefficients (including the interaction) are jointly zero.
demand_idx <- grep("^demand_", names(coef(model_final)))
H_demand   <- diag(length(coef(model_final)))[demand_idx, ]
wald.test(b = coef(model_final), Sigma = vcov(model_final), L = H_demand)
# We reject H0 -> demand type significantly affects the probability of protest violence.

# Wald test: all region dummies jointly = 0
# H0: region does not affect protest violence (all region coefficients = 0).
region_idx <- grep("^region_", names(coef(model_final)))
H_region   <- diag(length(coef(model_final)))[region_idx, ]
wald.test(b = coef(model_final), Sigma = vcov(model_final), L = H_region)
# We reject H0 -> region significantly affects protest violence.

# Wald test: state violence variables jointly = 0
# H0: past state repression (lag + cumulative + interaction) has no effect.
sv_idx <- grep("state_violence|cum_state", names(coef(model_final)))
H_sv   <- diag(length(coef(model_final)))[sv_idx, ]
wald.test(b = coef(model_final), Sigma = vcov(model_final), L = H_sv)
# We reject H0 -> state violence history significantly predicts protest violence.

# Wald test: all year dummies jointly = 0
# H0: no time trend (all year coefficients = 0).
year_idx <- grep("^year_", names(coef(model_final)))
H_year   <- diag(length(coef(model_final)))[year_idx, ]
wald.test(b = coef(model_final), Sigma = vcov(model_final), L = H_year)
# We reject H0 -> year fixed effects are jointly significant, time trends matter.

# Wald test: protest size variables jointly = 0
# H0: protest size (log_protest_size + large_protest + interaction) has no effect.
size_idx <- grep("protest_size|large_protest", names(coef(model_final)))
H_size   <- diag(length(coef(model_final)))[size_idx, ]
wald.test(b = coef(model_final), Sigma = vcov(model_final), L = H_size)
# We reject H0 -> protest size significantly affects violence probability.

#MARGINAL EFFECTS
# Marginal effects at the mean
logitmfx(model_final$formula, data = data, atmean = TRUE,
         robust = TRUE, clustervar1 = "country")

# Average marginal effects (averaged over all observations)
logitmfx(model_final$formula, data = data, atmean = FALSE,
         robust = TRUE, clustervar1 = "country")


