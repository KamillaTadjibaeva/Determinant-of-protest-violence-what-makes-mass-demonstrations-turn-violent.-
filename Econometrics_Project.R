# install.packages("mfx")
# install.packages(c('Hmisc', 'polspline'))
# install.packages('biglm')
# install.packages(c('rms', 'statmod', 'speedglm'))
# install.packages("countrycode")

library("sandwich")
library("lmtest")
library("MASS")
library("mfx")
library("aod")
library(dplyr)
library("logistf")
library(countrycode)

data <- read.csv("dataverse_files/mmALL_073120_csv.csv", stringsAsFactors = FALSE)
str(data)

# Keep only protest events (non-protest rows have no meaningful DV)
data <- data[data$protest == 1, ]

data$country <- factor(data$country)
# Region factor with labels
data$region <- factor(data$region)
                      # levels = 1:7,
                      # labels = c("South America", "Central America",
                      #            "North America", "Europe", "Asia",
                      #            "MENA", "Africa"))

# Convert start/end dates
data$start_date <- as.Date(paste(data$startyear, data$startmonth, data$startday, sep = "-"),
                           format = "%Y-%m-%d")
data$end_date <- as.Date(paste(data$endyear, data$endmonth, data$endday, sep = "-"),
                         format = "%Y-%m-%d")

# Protest duration in days
data$duration_days <- as.numeric(data$end_date - data$start_date)
data$duration_days[data$duration_days < 0] <- 0

# Midpoint encoding for participants_category
category_midpoints <- c(
  "50-99"      = 75,
  "100-999"    = 550,
  "1000-1999"  = 1500,
  "2000-4999"  = 3500,
  "5000-10000" = 7500,
  ">10000"     = 15000
)

# encode participants_category via midpoints
data$size_from_category <- category_midpoints[data$participants_category]

# clean raw participants column (handles "1000s", "100s", plain numbers)
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

# unified column — category first, raw as backup
data$protest_size <- ifelse(!is.na(data$size_from_category),
                            data$size_from_category,
                            data$size_from_raw)

par(mfrow = c(1, 2))

# skew problem
hist(data$protest_size[!is.na(data$protest_size)], breaks = 50,
     col = "steelblue", border = "white",
     main = "Protest Size (Raw)", xlab = "Participants")
abline(v = median(data$protest_size, na.rm = TRUE), col = "red", lwd = 2, lty = 2)

# Log — why we transform
hist(log(data$protest_size[!is.na(data$protest_size)] + 1), breaks = 30,
     col = "darkgreen", border = "white",
     main = "Protest Size (Log)", xlab = "log(Participants + 1)")

par(mfrow = c(1, 1))

# Impute missing protest_size with country-specific median + add missing indicator
data$protest_size_missing <- ifelse(is.na(data$protest_size), 1, 0)

country_medians <- data %>%
  group_by(country) %>%
  summarise(med_size = median(protest_size, na.rm = TRUE), .groups = "drop")

data <- data %>%
  left_join(country_medians, by = "country") %>%
  mutate(
    protest_size = ifelse(is.na(protest_size), med_size, protest_size),
    # If a country has NO observed sizes at all, fall back to global median
    protest_size = ifelse(is.na(protest_size), median(protest_size, na.rm = TRUE), protest_size)
  ) %>%
  select(-med_size)

# Log-transform for model input
data$log_protest_size <- log(data$protest_size + 1)

# Check coverage
cat("Category available:", sum(!is.na(data$size_from_category)), "\n")
cat("Raw filled gaps:", sum(is.na(data$size_from_category) & !is.na(data$size_from_raw)), "\n")
cat("Imputed (missing indicator = 1):", sum(data$protest_size_missing), "\n")

data$size_from_category <- NULL
data$size_from_raw <- NULL

# Creating demand dummies by searching text across ALL 4 demand columns
demand_cols <- c("protesterdemand1", "protesterdemand2", "protesterdemand3", "protesterdemand4")

data$demand_labor <- apply(data[, demand_cols], 1, function(row) {
  as.integer(any(grepl("labor|wage", row, ignore.case = TRUE)))
})

data$demand_land <- apply(data[, demand_cols], 1, function(row) {
  as.integer(any(grepl("land|farm", row, ignore.case = TRUE)))
})

data$demand_police <- apply(data[, demand_cols], 1, function(row) {
  as.integer(any(grepl("police|brutal", row, ignore.case = TRUE)))
})

data$demand_political <- apply(data[, demand_cols], 1, function(row) {
  as.integer(any(grepl("political|process", row, ignore.case = TRUE)))
})

data$demand_prices <- apply(data[, demand_cols], 1, function(row) {
  as.integer(any(grepl("price|tax", row, ignore.case = TRUE)))
})

data$demand_corruption <- apply(data[, demand_cols], 1, function(row) {
  as.integer(any(grepl("remov|corrupt", row, ignore.case = TRUE)))
})

data$demand_social <- apply(data[, demand_cols], 1, function(row) {
  as.integer(any(grepl("social", row, ignore.case = TRUE)))
})

# Number of demands per protest (count non-empty demand columns)
data$n_demands <- apply(data[, demand_cols], 1, function(row) {
  sum(row != "" & !is.na(row))
})

# Large protest indicator (>= 1000 participants)
data$large_protest <- ifelse(data$protest_size >= 1000, 1, 0)

# State used violence at this event? (needed ONLY to compute lags)
response_cols <- c("stateresponse1", "stateresponse2", "stateresponse3",
                   "stateresponse4", "stateresponse5", "stateresponse6", "stateresponse7")

data$violent_response <- apply(data[, response_cols], 1, function(row) {
  as.integer(any(grepl("kill|shoot|beat", row, ignore.case = TRUE)))
})

# Compute lags (within country, ordered by date)
data <- data[order(data$country, data$start_date), ]

data <- data %>%
  group_by(country) %>%
  mutate(
    lag1_state_violence = lag(violent_response, 1),
    cum_state_violence = lag(cummean(violent_response), 1)
  ) %>%
  ungroup()

# First protest in each country has no history — fill with 0 (no prior info)
data$lag1_state_violence[is.na(data$lag1_state_violence)] <- 0
data$cum_state_violence[is.na(data$cum_state_violence)] <- 0

# Remove the temporary column
data$violent_response <- NULL

cat("Lag1 state violence rate:", mean(data$lag1_state_violence), "\n")
cat("Cumulative state violence mean:", mean(data$cum_state_violence), "\n")

# Year as factor for fixed effects
data$year_factor <- factor(data$year)

cat("Final dataset:", nrow(data), "protest events,",
    length(unique(data$country)), "countries,",
    paste(range(data$year), collapse = "-"), "\n")
cat("Violent protests:", sum(data$protesterviolence, na.rm = TRUE),
    "(", round(mean(data$protesterviolence, na.rm = TRUE) * 100, 1), "%)\n")

# Drop columns no longer needed (already encoded or unused)
drop_cols <- c("id", "ccode", "protest",
               "startday", "startmonth", "startyear",
               "endday", "endmonth", "endyear",
               "start_date", "end_date",
               "location", "participants_category", "participants",
               "protesteridentity",
               "protesterdemand1", "protesterdemand2", "protesterdemand3", "protesterdemand4",
               "stateresponse1", "stateresponse2", "stateresponse3",
               "stateresponse4", "stateresponse5", "stateresponse6", "stateresponse7",
               "sources", "notes")
data <- data[, !(names(data) %in% drop_cols)]

View(data)

# Model 1: Baseline — what predicts protester violence?
model1 <- glm(protesterviolence ~
                log_protest_size +
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

summary(model1)

# Marginal effects 
logitmfx(protesterviolence ~
           log_protest_size +
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
         data = data, robust = TRUE, clustervar1 = "country")

# Model 2: Adding large protest indicator
model2 <- glm(protesterviolence ~
                log_protest_size +
                large_protest +
                protest_size_missing +
                duration_days +
                demand_political +
                demand_labor +
                demand_land +
                demand_corruption +
                n_demands +
                protestnumber +
                lag1_state_violence +
                cum_state_violence +
                region +
                year_factor,
              data = data,
              family = binomial(link = "logit"))

summary(model2)
coeftest(model2, vcov = vcovCL(model2, cluster = data$country))
