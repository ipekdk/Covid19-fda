library(tidyverse) # programming tools
library(fda) # functional data analysis
library(curl) # for downloading data
library(functional) # functional programming tools
library(countrycode) # country code schemes (NOTE: PACKAGE IS NEW)

### GET DATA AND REMOVE COUNTRIES WHERE POPULATION DATA IS NOT PRESENT ###
country_data <- read.csv(curl("https://datahub.io/core/covid-19/r/countries-aggregated.csv"), stringsAsFactors = F) # data
population_data <- world_bank_pop # population data from the world bank (dating to 2017)
population_totals <- subset(population_data, indicator == "SP.POP.TOTL") # get population totals
population_2017 <- select(population_totals, c("country", "2017")) # select most recent year
country_data$Country <- countrycode(country_data$Country, "country.name", "wb") # convert country names to world bank format where possible
country_data <- subset(country_data, is.na(country_data$Country) == F) # drop rows where conversion was not possible
country_data <- subset(country_data, country_data$Country %in% population_2017$country) # drop rows where population data is missing
countries <- sort(unique(country_data$Country)) # these are the countries for which we could convert format and have population data on

### EXTRACT CONFIRMED, DEATHS, AND RECOVERY DATA ###
confirmed <- list()
deaths <- list()
recovered <- list()

for(place in countries){
  confirmed[[place]] <- subset(country_data, country_data$Country == place)$Confirmed
  deaths[[place]] <- subset(country_data, country_data$Country == place)$Deaths
  recovered[[place]] <- subset(country_data, country_data$Country == place)$Recovered
}

### SCALE CONFIRMED AND DEATHS BY POPULATION SIZE ###
confirmed_scaled <- list()
deaths_scaled <- list()

for(place in countries){
  confirmed_scaled[[place]] <- (10^6) * (confirmed[[place]]/subset(population_2017, population_2017$country == place)$`2017`)
  deaths_scaled[[place]] <- (10^6) * (deaths[[place]]/subset(population_2017, population_2017$country == place)$`2017`)
}

### COMPUTE RATIOS ###
case_fatality_ratios <- list()
case_recovery_ratios <- list()

for(place in countries){
  case_fatality_ratios[[place]] <- 100 * (deaths[[place]]/max(confirmed[[place]], 1))
  case_recovery_ratios[[place]] <- 100 * (recovered[[place]]/max(confirmed[[place]], 1))
}

### CONSTRUCT BASIS, SMOOTH DATA INTO BASIS, AND EXTRACT FD OBJECTS ###
days <- 1:length(confirmed$AFG) # days from January 22nd (inclusive)
spline_basis <- create.bspline.basis(rangeval = range(days), nbasis = 10) # B-spline basis with 10 basis functions
fd_builder <- function(x){as.fd(smooth.basis(argvals = days, y = x, spline_basis))}

confirmed_functions_fd <- lapply(confirmed_scaled, fd_builder)
death_functions_fd <- lapply(deaths_scaled, fd_builder)
cfr_functions_fd <- lapply(case_fatality_ratios, fd_builder)
crr_functions_fd <- lapply(case_recovery_ratios, fd_builder)

### SUBSET SELECTED COUNTRIES FOR CONSIDERATION ###

places <- sort(c("China", "South Korea", "United Kingdom", "France", "United States", "Germany", "Spain", "Italy"))
places_codes <- countrycode(places, "country.name", "wb")

# NB: If you want more than 8 countries, you will need to expand R's colour scheme

selected_confirmed <- list()
selected_deaths <- list()
selected_cfr <- list()
selected_crr <- list()

for(place in places_codes){
  selected_confirmed[[place]] <- confirmed_functions_fd[[place]]
  selected_deaths[[place]] <- death_functions_fd[[place]]
  selected_cfr[[place]] <- cfr_functions_fd[[place]]
  selected_crr[[place]] <- crr_functions_fd[[place]]
}

# OBTAIN MAXIMUMS AND MINIMUMS TO HAVE A DATA-ADAPTIVE Y-AXIS WHEN PLOTTING
getmax <- function(derivative, x){max(eval.fd(days, deriv(x, Lfdobj = derivative)))}
getmax_deriv0 <- Curry(getmax, 0)
getmax_deriv1 <- Curry(getmax, 1)
getmax_deriv2 <- Curry(getmax, 2)

getmin <- function(derivative, x){min(eval.fd(days, deriv(x, Lfdobj = derivative)))}
getmin_deriv0 <- Curry(getmin, 0)
getmin_deriv1 <- Curry(getmin, 1)
getmin_deriv2 <- Curry(getmin, 2)

max_confirmed_deriv0 <- max(unlist(lapply(selected_confirmed, getmax_deriv0), use.names = F))
max_confirmed_deriv1 <- max(unlist(lapply(selected_confirmed, getmax_deriv1), use.names = F))
max_confirmed_deriv2 <- max(unlist(lapply(selected_confirmed, getmax_deriv2), use.names = F))
max_deaths_deriv0 <- max(unlist(lapply(selected_deaths, getmax_deriv0), use.names = F))
max_deaths_deriv1 <- max(unlist(lapply(selected_deaths, getmax_deriv1), use.names = F))
max_deaths_deriv2 <- max(unlist(lapply(selected_deaths, getmax_deriv2), use.names = F))
max_cfr <- max(unlist(lapply(selected_cfr, getmax_deriv0), use.names = F))
max_crr <- max(unlist(lapply(selected_crr, getmax_deriv0), use.names = F))

min_confirmed_deriv0 <- min(unlist(lapply(selected_confirmed, getmin_deriv0), use.names = F))
min_confirmed_deriv1 <- min(unlist(lapply(selected_confirmed, getmin_deriv1), use.names = F))
min_confirmed_deriv2 <- min(unlist(lapply(selected_confirmed, getmin_deriv2), use.names = F))
min_deaths_deriv0 <- min(unlist(lapply(selected_deaths, getmin_deriv0), use.names = F))
min_deaths_deriv1 <- min(unlist(lapply(selected_deaths, getmin_deriv1), use.names = F))
min_deaths_deriv2 <- min(unlist(lapply(selected_deaths, getmin_deriv2), use.names = F))
min_cfr <- min(unlist(lapply(selected_cfr, getmin_deriv0), use.names = F))
min_crr <- min(unlist(lapply(selected_crr, getmin_deriv0), use.names = F))

### PLOTS ###
colours <- palette() # R's default colour scheme

plot.fd(selected_confirmed[[1]], col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Confirmed cases (per million)", main = "Confirmed cases of Covid-19 (per million)",
        ylim = 1.05 * c(min_confirmed_deriv0, max_confirmed_deriv0))
for(i in 2:length(places)){lines(selected_confirmed[[i]], col = colours[i])}
legend("topleft", legend = countrycode(places_codes, "wb", "country.name"), col = colours, lty = 1, cex = 0.7)

plot.fd(deriv.fd(selected_confirmed[[1]]), col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Value of derivative", main = "Derivative of confirmed cases (per million) curve",
        ylim = 1.05 * c(min_confirmed_deriv1, max_confirmed_deriv1))
for(i in 2:length(places)){lines(deriv.fd(selected_confirmed[[i]]), col = colours[i])}
legend("topleft", legend = countrycode(places_codes, "wb", "country.name"), col = colours, lty = 1, cex = 0.7)

plot.fd(deriv.fd(selected_confirmed[[1]], Lfdobj = 2), col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Value of second derivative", main = "Second derivative of confirmed cases (per million) curve",
        ylim = 1.05 * c(min_confirmed_deriv2, max_confirmed_deriv2))
for(i in 2:length(places)){lines(deriv.fd(selected_confirmed[[i]], Lfdobj = 2), col = colours[i])}
legend("topleft", legend = countrycode(places_codes, "wb", "country.name"), col = colours, lty = 1, cex = 0.7)

# DEATHS

plot.fd(selected_deaths[[1]], col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Deaths (per million)", main = "Deaths from Covid-19 (per million)",
        ylim = 1.05 * c(min_deaths_deriv0, max_deaths_deriv0))
for(i in 2:length(places)){lines(selected_deaths[[i]], col = colours[i])}
legend("topleft", legend = countrycode(places_codes, "wb", "country.name"), col = colours, lty = 1, cex = 0.7)

plot.fd(deriv.fd(selected_deaths[[1]]), col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Value of derivative", main = "Derivative of deaths (per million) curve",
        ylim = 1.05 * c(min_deaths_deriv1, max_deaths_deriv1))
for(i in 2:length(places)){lines(deriv.fd(selected_deaths[[i]]), col = colours[i])}
legend("topleft", legend = countrycode(places_codes, "wb", "country.name"), col = colours, lty = 1, cex = 0.7)

plot.fd(deriv.fd(selected_deaths[[1]], Lfdobj = 2), col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Value of second derivative", main = "Second derivative of deaths (per million) curve",
        ylim = 1.05 * c(min_deaths_deriv2, max_deaths_deriv2))
for(i in 2:length(places)){lines(deriv.fd(selected_deaths[[i]], Lfdobj = 2), col = colours[i])}
legend("topleft", legend = countrycode(places_codes, "wb", "country.name"), col = colours, lty = 1, cex = 0.7)

# RATIOS

plot.fd(selected_cfr[[1]], col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Case fatality ratio in %", main = "Case fatality ratio",
        ylim = 1.05 * c(min_cfr, max_cfr))
for(i in 2:length(places)){lines(selected_cfr[[i]], col = colours[i])}
legend("topleft", legend = countrycode(places_codes, "wb", "country.name"), col = colours, lty = 1, cex = 0.7)

plot.fd(selected_crr[[1]], col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Case recovery ratio in %", main = "Case recovery ratio",
        ylim = 1.05 * c(min_crr, max_crr))
for(i in 2:length(places)){lines(selected_crr[[i]], col = colours[i])}
legend("topleft", legend = countrycode(places_codes, "wb", "country.name"), col = colours, lty = 1, cex = 0.7)
