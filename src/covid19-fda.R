library(tidyverse)
library(fda) # functional data analysis
library(curl) # for downloading data

# https://data.humdata.org/dataset/5dff64bc-a671-48da-aa87-2ca40d7abf02 (Source of Case Data)

### GET DATA ###

confirmed <- read_csv(curl("https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_confirmed_global.csv&filename=time_series_covid19_confirmed_global.csv"))

deaths <- read_csv(curl("https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_deaths_global.csv&filename=time_series_covid19_deaths_global.csv"))

### DROP LATITUDE AND LONGITUDE ###

confirmed <- subset(confirmed, select = -c(Lat, Long))
deaths <- subset(deaths, select = -c(Lat, Long))

### RENAME COLUMNS FOR CONVENIENCE ###
days <- 1:(ncol(confirmed) - 2)
names(confirmed) <- c("province", "country", days)
names(deaths) <- c("province", "country", days)

### COUNTRIES TO BE COMPARED ###

countries <- sort(c("US", "United Kingdom", "Italy", "Germany", "China", "Korea, South", "France", "Spain"))

confirmed <- subset(confirmed, confirmed$country %in% countries)
deaths <- subset(deaths, deaths$country %in% countries)

### DATA PREPROCESSING ###

# Each country should have its own column
confirmed_processed <- list()
deaths_processed <- list()

for(place in countries){
  confirmed_processed[[place]] <- NA
  deaths_processed[[place]] <- NA
}

for(place in countries){
  confirmed_subsetted <- subset(confirmed, confirmed$country == place)
  confirmed_index <- match(NA, confirmed_subsetted$province)
  deaths_subsetted <- subset(deaths, deaths$country == place)
  deaths_index <- match(NA, deaths_subsetted$province)
  if(!is.na(confirmed_index)){
    # Ignore any overseas terrorities
    confirmed_processed[[place]] <- as.vector(confirmed_subsetted[confirmed_index, 3:(max(days) + 2)], mode = "numeric")
    deaths_processed[[place]] <- as.vector(deaths_subsetted[deaths_index, 3:(max(days) + 2)], mode = "numeric")
  }else{
    # Aggregate across provinces (this approach is fine owing to how the original data has been structured)
    confirmed_processed[[place]] <- as.vector(apply(confirmed_subsetted[, 3:(max(days) + 2)], 2, sum), mode = "numeric")
    deaths_processed[[place]] <- as.vector(apply(deaths_subsetted[, 3:(max(days) + 2)], 2, sum), mode = "numeric")
    }
}

### Calculate Ratios
ratios <- list()
for(place in countries){
  ratios[[place]] <- 100 * (deaths_processed[[place]]/max(confirmed_processed[[place]], 1))
}

### CREATE SPLINE BASIS ###
spline_basis <- create.bspline.basis(rangeval = range(days), nbasis = 10) ### B-spline basis with 10 basis functions

### COMPUTE LOGARITHMS (ADDING 1 TO AVOID LOG(0))
plus_one_log <- function(x){log(x + 1)}

log_confirmed <- lapply(confirmed_processed, plus_one_log)
log_deaths <- lapply(deaths_processed, plus_one_log)

### SMOOTH DATA INTO THE BASIS ###
smoother <- function(x){smooth.basis(argvals = days, y = x, spline_basis)}

confirmed_functions <- lapply(confirmed_processed, smoother)
death_functions <- lapply(deaths_processed, smoother)
ratio_functions <- lapply(ratios, smoother)

log_confirmed_functions <- lapply(log_confirmed, smoother)
log_death_functions <- lapply(log_deaths, smoother)

### EXTRACT THE FUNCTIONAL DATA OBJECTS ###
extracter <- function(x){x$fd}

confirmed_functions_fd <- lapply(confirmed_functions, extracter)
death_functions_fd <- lapply(death_functions, extracter)
ratios_fd <- lapply(ratio_functions, extracter)
log_confirmed_functions_fd <- lapply(log_confirmed_functions, extracter)
log_death_functions_fd <- lapply(log_death_functions, extracter)

### PLOTS ###
colours <- c('black', 'blue', 'green', 'orange', 'purple', 'red', 'violet', 'yellow')

# CONFIRMED

plot.fd(confirmed_functions_fd[[1]], col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Confirmed cases", main = "Confirmed cases of Covid-19", ylim = c(0, 100000))
for(i in 2:length(countries)){lines(confirmed_functions_fd[[i]], col = colours[i])}
legend("topleft", legend = countries, col = colours, lty = 1, cex = 0.7)

plot.fd(deriv.fd(confirmed_functions_fd[[1]]), col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Value of derivative", main = "Derivative of confirmed cases curve", ylim = c(-100, 20000))
for(i in 2:length(countries)){lines(deriv.fd(confirmed_functions_fd[[i]]), col = colours[i])}
legend("topleft", legend = countries, col = colours, lty = 1, cex = 0.7)

plot.fd(deriv.fd(confirmed_functions_fd[[1]], Lfdobj = 2), col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Value of second derivative", main = "Second derivative of confirmed cases curve", ylim = c(-600, 3000))
for(i in 2:length(countries)){lines(deriv.fd(confirmed_functions_fd[[i]], Lfdobj = 2), col = colours[i])}
legend("topleft", legend = countries, col = colours, lty = 1, cex = 0.7)

# DEATHS

plot.fd(death_functions_fd[[1]], col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Deaths", main = "Deaths from Covid-19", ylim = c(0, 9000))
for(i in 2:length(countries)){lines(death_functions_fd[[i]], col = colours[i])}
legend("topleft", legend = countries, col = colours, lty = 1, cex = 0.7)

plot.fd(deriv.fd(death_functions_fd[[1]]), col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Value of derivative", main = "Derivative of death curve", ylim = c(0, 1000))
for(i in 2:length(countries)){lines(deriv.fd(death_functions_fd[[i]]), col = colours[i])}
legend("topleft", legend = countries, col = colours, lty = 1, cex = 0.7)

plot.fd(deriv.fd(death_functions_fd[[1]], Lfdobj = 2), col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Value of second derivative", main = "Second derivative of death curve", ylim = c(-10, 150))
for(i in 2:length(countries)){lines(deriv.fd(death_functions_fd[[i]], Lfdobj = 2), col = colours[i])}
legend("topleft", legend = countries, col = colours, lty = 1, cex = 0.7)


# LOG CONFIRMED

plot.fd(log_confirmed_functions_fd[[1]], col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Log confirmed cases", main = "Log confirmed cases of Covid-19", ylim = c(0, 12))
for(i in 2:length(countries)){lines(log_confirmed_functions_fd[[i]], col = colours[i])}
legend("topleft", legend = countries, col = colours, lty = 1, cex = 0.7)

plot.fd(deriv.fd(log_confirmed_functions_fd[[1]]), col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Value of derivative", main = "Derivative of log confirmed cases curve", ylim = c(-0.5, 0.7))
for(i in 2:length(countries)){lines(deriv.fd(log_confirmed_functions_fd[[i]]), col = colours[i])}
legend("bottomright", legend = countries, col = colours, lty = 1, cex = 0.7)

plot.fd(deriv.fd(log_confirmed_functions_fd[[1]], Lfdobj = 2), col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Value of second derivative", main = "Second derivative of log confirmed cases curve", ylim = c(-0.2, 0.2))
for(i in 2:length(countries)){lines(deriv.fd(log_confirmed_functions_fd[[i]], Lfdobj = 2), col = colours[i])}
legend("bottomright", legend = countries, col = colours, lty = 1, cex = 0.7)

# LOG DEATHS

plot.fd(log_death_functions_fd[[1]], col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Log deaths", main = "Log deaths from Covid-19", ylim = c(0, 10))
for(i in 2:length(countries)){lines(log_death_functions_fd[[i]], col = colours[i])}
legend("topleft", legend = countries, col = colours, lty = 1, cex = 0.7)

plot.fd(deriv.fd(log_death_functions_fd[[1]]), col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Value of derivative", main = "Derivative of log deaths curve", ylim = c(-0.5, 0.7))
for(i in 2:length(countries)){lines(deriv.fd(log_death_functions_fd[[i]]), col = colours[i])}
legend("topleft", legend = countries, col = colours, lty = 1, cex = 0.7)

plot.fd(deriv.fd(log_death_functions_fd[[1]], Lfdobj = 2), col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Value of second derivative", main = "Second derivative of log deaths curve", ylim = c(-0.1, 0.1))
for(i in 2:length(countries)){lines(deriv.fd(log_death_functions_fd[[i]], Lfdobj = 2), col = colours[i])}
legend("topleft", legend = countries, col = colours, lty = 1, cex = 0.7)

### RATIOS
plot.fd(ratios_fd[[1]], col = colours[1], xlab = "Days from January 22nd (Inclusive)",
        ylab = "Deaths/Confirmed Cases in %", main = "Case fatality ratio", ylim = c(0, 11))
for(i in 2:length(countries)){lines(ratios_fd[[i]], col = colours[i])}
legend("topleft", legend = countries, col = colours, lty = 1, cex = 0.7)
