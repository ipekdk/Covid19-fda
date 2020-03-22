library(tidyverse)
library(fda) # functional data analysis

# https://data.humdata.org/dataset/5dff64bc-a671-48da-aa87-2ca40d7abf02 (Source of Data)

### GET DATA ### [CHANGE AS NEEDED]
confirmed_cases <- read_csv('/home/mrlj/Desktop/covid-19/fda-covid19/data/time_series_2019-ncov-Confirmed.csv')
death_cases <- read_csv('/home/mrlj/Desktop/covid-19/fda-covid19/data/time_series_2019-ncov-Deaths.csv')
recovered_cases <- read_csv('/home/mrlj/Desktop/covid-19/fda-covid19/data/time_series_2019-ncov-Recovered.csv')

### COUNTRIES TO BE COMPARED ###

countries <- c("US", "United Kingdom", "Italy", "Germany", "China", "Korea, South", "France")

### CREATE SPLINE BASIS ###
days <- seq(1, nrow(confirmed_cases), 1) # days from January 22nd
spline_basis <- create.bspline.basis(rangeval = range(days), nbasis = 10) ### B-spline basis with 10 basis functions

### COMPUTE LOGARITHMS (ADDING 1 TO AVOID LOG(0))
log_confirmed_cases <- list()
log_death_cases <- list()
log_recovered_cases <- list()

for(country in countries){
  log_confirmed_cases[[country]] <- log(confirmed_cases[[country]] + 1)
  log_death_cases[[country]] <- log(death_cases[[country]] + 1)
  log_recovered_cases[[country]] <- log(recovered_cases[[country]] + 1)
}

### SMOOTH DATA INTO THE BASIS ###
confirmed_functions <- list()
death_functions <- list()
recovered_functions <- list()
log_confirmed_functions <- list()
log_death_functions <- list()
log_recovered_functions <- list()

for(country in countries){
  confirmed_functions[[country]] <- smooth.basis(argvals = days, y = confirmed_cases[[country]], spline_basis)
  death_functions[[country]] <- smooth.basis(argvals = days, y = death_cases[[country]], spline_basis)
  recovered_functions[[country]] <- smooth.basis(argvals = days, y = recovered_cases[[country]], spline_basis)
  log_confirmed_functions[[country]] <- smooth.basis(argvals = days, y = log_confirmed_cases[[country]], spline_basis)
  log_death_functions[[country]] <- smooth.basis(argvals = days, y = log_death_cases[[country]], spline_basis)
  log_recovered_functions[[country]] <- smooth.basis(argvals = days, y = log_recovered_cases[[country]], spline_basis)
}

### EXTRACT THE FUNCTIONAL DATA OBJECTS ###
confirmed_functions_fd <- list()
death_functions_fd <- list()
recovered_functions_fd <- list()
log_confirmed_functions_fd <- list()
log_death_functions_fd <- list()
log_recovered_functions_fd <- list()

for(country in countries){
  confirmed_functions_fd[[country]] <- confirmed_functions[[country]]$fd
  death_functions_fd[[country]] <- death_functions[[country]]$fd
  recovered_functions_fd[[country]] <- recovered_functions[[country]]$fd
  log_confirmed_functions_fd[[country]] <- log_confirmed_functions[[country]]$fd
  log_death_functions_fd[[country]] <- log_death_functions[[country]]$fd
  log_recovered_functions_fd[[country]] <- log_recovered_functions[[country]]$fd
}

### PLOTS ###
# CONFIRMED
plot.fd(confirmed_functions_fd$`United Kingdom`, col = 'black', xlab = "Days from January 22nd", ylab = "Confirmed cases",
        ylim = c(0, 80000), main = "Confirmed cases of Covid-19")
lines(confirmed_functions_fd$`Korea, South`, col = 'blue')
lines(confirmed_functions_fd$China, col = 'green')
lines(confirmed_functions_fd$US, col = 'orange')
lines(confirmed_functions_fd$Italy, col = 'purple')
lines(confirmed_functions_fd$Germany, col = 'red')
lines(confirmed_functions_fd$France, col = 'violet')
legend(1, 60000, legend=c("UK", "South Korea", "China", "US", "Italy", "Germany", "France"),
       col=c("black", "blue", "green", "orange", "purple", "red", "violet"), lty=1, cex=0.8)

plot.fd(confirmed_functions_fd$`United Kingdom`, col = 'black', Lfdobj = 1, xlab = "Days from January 22nd",
        ylab = "Value of derivative", ylim = c(-100, 6000), main = "Derivative of confirmed cases curve")
lines(deriv.fd(confirmed_functions_fd$`Korea, South`), col = 'blue')
lines(deriv.fd(confirmed_functions_fd$China), col = 'green')
lines(deriv.fd(confirmed_functions_fd$US), col = 'orange')
lines(deriv.fd(confirmed_functions_fd$Italy), col = 'purple')
lines(deriv.fd(confirmed_functions_fd$Germany), col = 'red')
lines(deriv.fd(confirmed_functions_fd$France), col = 'violet')
legend(1, 6000, legend=c("UK", "South Korea", "China", "US", "Italy", "Germany", "France"),
       col=c("black", "blue", "green", "orange", "purple", "red", "violet"), lty=1, cex=0.8)

plot.fd(confirmed_functions_fd$`United Kingdom`, col = 'black', Lfdobj = 2, xlab = "Days from January 22nd",
        ylab = "Value of second derivative", ylim = c(-1000, 1000), main = "Second derivative of confirmed cases curve")
lines(deriv.fd(confirmed_functions_fd$`Korea, South`, Lfdobj = 2), col = 'blue')
lines(deriv.fd(confirmed_functions_fd$China, Lfdobj = 2), col = 'green')
lines(deriv.fd(confirmed_functions_fd$US, Lfdobj = 2), col = 'orange')
lines(deriv.fd(confirmed_functions_fd$Italy, Lfdobj = 2), col = 'purple')
lines(deriv.fd(confirmed_functions_fd$Germany, Lfdobj = 2), col = 'red')
lines(deriv.fd(confirmed_functions_fd$France, Lfdobj = 2), col = 'violet')
legend(1, -100, legend=c("UK", "South Korea", "China", "US", "Italy", "Germany", "France"),
       col=c("black", "blue", "green", "orange", "purple", "red", "violet"), lty=1, cex=0.8)

# DEATHS
plot.fd(death_functions_fd$`United Kingdom`, col = 'black', xlab = "Days from January 22nd", ylab = "Deaths", ylim = c(0, 5000),
        main = "Deaths from Covid-19")
lines(death_functions_fd$`Korea, South`, col = 'blue')
lines(death_functions_fd$China, col = 'green')
lines(death_functions_fd$US, col = 'orange')
lines(death_functions_fd$Italy, col = 'purple')
lines(death_functions_fd$Germany, col = 'red')
lines(death_functions_fd$France, col = 'violet')
legend(1, 5000, legend=c("UK", "South Korea", "China", "US", "Italy", "Germany", "France"),
       col=c("black", "blue", "green", "orange", "purple", "red", "violet"), lty=1, cex=0.8)

plot.fd(death_functions_fd$`United Kingdom`, col = 'black', Lfdobj = 1, xlab = "Days from January 22nd",
        ylab = "Value of derivative", ylim = c(0, 600), main = "Derivative of death curve")
lines(deriv.fd(death_functions_fd$`Korea, South`), col = 'blue')
lines(deriv.fd(death_functions_fd$China), col = 'green')
lines(deriv.fd(death_functions_fd$US), col = 'orange')
lines(deriv.fd(death_functions_fd$Italy), col = 'purple')
lines(deriv.fd(death_functions_fd$Germany), col = 'red')
lines(deriv.fd(death_functions_fd$France), col = 'violet')
legend(1, 500, legend=c("UK", "South Korea", "China", "US", "Italy", "Germany", "France"),
       col=c("black", "blue", "green", "orange", "purple", "red", "violet"), lty=1, cex=0.8)

plot.fd(death_functions_fd$`United Kingdom`, col = 'black', Lfdobj = 2, xlab = "Days from January 22nd",
        ylab = "Value of second derivative", ylim = c(-20, 70), main = "Second derivative of death curve")
lines(deriv.fd(death_functions_fd$`Korea, South`, Lfdobj = 2), col = 'blue')
lines(deriv.fd(death_functions_fd$China, Lfdobj = 2), col = 'green')
lines(deriv.fd(death_functions_fd$US, Lfdobj = 2), col = 'orange')
lines(deriv.fd(death_functions_fd$Italy, Lfdobj = 2), col = 'purple')
lines(deriv.fd(death_functions_fd$Germany, Lfdobj = 2), col = 'red')
lines(deriv.fd(death_functions_fd$France, Lfdobj = 2), col = 'violet')
legend(1, 50, legend=c("UK", "South Korea", "China", "US", "Italy", "Germany", "France"),
       col=c("black", "blue", "green", "orange", "purple", "red", "violet"), lty=1, cex=0.8)

# LOG CONFIRMED
plot.fd(log_confirmed_functions_fd$`United Kingdom`, col = 'black', xlab = "Days from January 22nd", ylab = "Log confirmed cases",
        ylim = c(0, 12), main = "Log confirmed cases of Covid-19")
lines(log_confirmed_functions_fd$`Korea, South`, col = 'blue')
lines(log_confirmed_functions_fd$China, col = 'green')
lines(log_confirmed_functions_fd$US, col = 'orange')
lines(log_confirmed_functions_fd$Italy, col = 'purple')
lines(log_confirmed_functions_fd$Germany, col = 'red')
lines(log_confirmed_functions_fd$France, col = 'violet')
legend(45, 4.5, legend=c("UK", "South Korea", "China", "US", "Italy", "Germany", "France"),
       col=c("black", "blue", "green", "orange", "purple", "red", "violet"), lty=1, cex=0.7)

plot.fd(log_confirmed_functions_fd$`United Kingdom`, col = 'black', Lfdobj = 1, xlab = "Days from January 22nd",
        ylab = "Value of derivative", ylim = c(-0.5, 0.7), main = "Derivative of log confirmed cases curve")
lines(deriv.fd(log_confirmed_functions_fd$`Korea, South`), col = 'blue')
lines(deriv.fd(log_confirmed_functions_fd$China), col = 'green')
lines(deriv.fd(log_confirmed_functions_fd$US), col = 'orange')
lines(deriv.fd(log_confirmed_functions_fd$Italy), col = 'purple')
lines(deriv.fd(log_confirmed_functions_fd$Germany), col = 'red')
lines(deriv.fd(log_confirmed_functions_fd$France), col = 'violet')
legend(40, -0.05, legend=c("UK", "South Korea", "China", "US", "Italy", "Germany", "France"),
       col=c("black", "blue", "green", "orange", "purple", "red", "violet"), lty=1, cex=0.8)

plot.fd(log_confirmed_functions_fd$`United Kingdom`, col = 'black', Lfdobj = 2, xlab = "Days from January 22nd",
        ylab = "Value of second derivative", ylim = c(-0.2, 0.2), main = "Second derivative of log confirmed cases curve")
lines(deriv.fd(log_confirmed_functions_fd$`Korea, South`, Lfdobj = 2), col = 'blue')
lines(deriv.fd(log_confirmed_functions_fd$China, Lfdobj = 2), col = 'green')
lines(deriv.fd(log_confirmed_functions_fd$US, Lfdobj = 2), col = 'orange')
lines(deriv.fd(log_confirmed_functions_fd$Italy, Lfdobj = 2), col = 'purple')
lines(deriv.fd(log_confirmed_functions_fd$Germany, Lfdobj = 2), col = 'red')
lines(deriv.fd(log_confirmed_functions_fd$France, Lfdobj = 2), col = 'violet')
legend(40, -0.05, legend=c("UK", "South Korea", "China", "US", "Italy", "Germany", "France"),
       col=c("black", "blue", "green", "orange", "purple", "red", "violet"), lty=1, cex=0.8)

# LOG DEATHS
plot.fd(log_death_functions_fd$`United Kingdom`, col = 'black', xlab = "Days from January 22nd", ylab = "Log death cases",
        ylim = c(0, 9), main = "Log death cases of Covid-19")
lines(log_death_functions_fd$`Korea, South`, col = 'blue')
lines(log_death_functions_fd$China, col = 'green')
lines(log_death_functions_fd$US, col = 'orange')
lines(log_death_functions_fd$Italy, col = 'purple')
lines(log_death_functions_fd$Germany, col = 'red')
lines(log_death_functions_fd$France, col = 'violet')
legend(1, 4.5, legend=c("UK", "South Korea", "China", "US", "Italy", "Germany", "France"),
       col=c("black", "blue", "green", "orange", "purple", "red", "violet"), lty=1, cex=0.7)

plot.fd(log_death_functions_fd$`United Kingdom`, col = 'black', Lfdobj = 1, xlab = "Days from January 22nd",
        ylab = "Value of derivative", ylim = c(-0.5, 0.7), main = "Derivative of log death cases curve")
lines(deriv.fd(log_death_functions_fd$`Korea, South`), col = 'blue')
lines(deriv.fd(log_death_functions_fd$China), col = 'green')
lines(deriv.fd(log_death_functions_fd$US), col = 'orange')
lines(deriv.fd(log_death_functions_fd$Italy), col = 'purple')
lines(deriv.fd(log_death_functions_fd$Germany), col = 'red')
lines(deriv.fd(log_death_functions_fd$France), col = 'violet')
legend(40, -0.05, legend=c("UK", "South Korea", "China", "US", "Italy", "Germany", "France"),
       col=c("black", "blue", "green", "orange", "purple", "red", "violet"), lty=1, cex=0.8)

plot.fd(log_death_functions_fd$`United Kingdom`, col = 'black', Lfdobj = 2, xlab = "Days from January 22nd",
        ylab = "Value of second derivative", ylim = c(-0.2, 0.2), main = "Second derivative of log death cases curve")
lines(deriv.fd(log_death_functions_fd$`Korea, South`, Lfdobj = 2), col = 'blue')
lines(deriv.fd(log_death_functions_fd$China, Lfdobj = 2), col = 'green')
lines(deriv.fd(log_death_functions_fd$US, Lfdobj = 2), col = 'orange')
lines(deriv.fd(log_death_functions_fd$Italy, Lfdobj = 2), col = 'purple')
lines(deriv.fd(log_death_functions_fd$Germany, Lfdobj = 2), col = 'red')
lines(deriv.fd(log_death_functions_fd$France, Lfdobj = 2), col = 'violet')
legend(40, -0.05, legend=c("UK", "South Korea", "China", "US", "Italy", "Germany", "France"),
       col=c("black", "blue", "green", "orange", "purple", "red", "violet"), lty=1, cex=0.8)
