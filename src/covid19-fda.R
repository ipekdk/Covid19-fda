library(fda) # functional data analysis

# https://data.humdata.org/dataset/5dff64bc-a671-48da-aa87-2ca40d7abf02 (Source of Data)

### GET DATA ### [CHANGE AS NEEDED]
confirmed_cases <- read_csv('/home/mrlj/Desktop/covid-19/fda-covid19/data/time_series_2019-ncov-Confirmed.csv')
death_cases <- read_csv('/home/mrlj/Desktop/covid-19/fda-covid19/data/time_series_2019-ncov-Deaths.csv')
recovered_cases <- read_csv('/home/mrlj/Desktop/covid-19/fda-covid19/data/time_series_2019-ncov-Recovered.csv')

### COUNTRIES TO BE COMPARED ###

countries <- c("US", "United Kingdom", "Italy", "Germany", "China", "Korea, South", "France")

### CREATE SPLINE BASIS ###
days <- seq(1, nrow(confirmed_cases), 1) # days from January 20th
spline_basis <- create.bspline.basis(rangeval = range(days), nbasis = 10) ### B-spline basis with 10 basis functions

### SMOOTH DATA INTO THE BASIS ###
confirmed_functions <- list()
death_functions <- list()
recovered_functions <- list()
cd_ratio_functions <- list()

for(country in countries){
  confirmed_functions[[country]] <- smooth.basis(argvals = days, y = confirmed_cases[[country]], spline_basis)
  death_functions[[country]] <- smooth.basis(argvals = days, y = death_cases[[country]], spline_basis)
  recovered_functions[[country]] <- smooth.basis(argvals = days, y = recovered_cases[[country]], spline_basis)}

### EXTRACT THE FUNCTIONAL DATA OBJECTS ###
confirmed_functions_fd <- list()
death_functions_fd <- list()
recovered_functions_fd <- list()

for(country in countries){
  confirmed_functions_fd[[country]] <- confirmed_functions[[country]]$fd
  death_functions_fd[[country]] <- death_functions[[country]]$fd
  recovered_functions_fd[[country]] <- recovered_functions[[country]]$fd
}

### PLOTS ###
plot.fd(death_functions_fd$`United Kingdom`, col = 'black', xlab = "Days from January 20th", ylab = "Deaths", ylim = c(0, 5000),
        main = "Deaths from Covid-19")
lines(death_functions_fd$`Korea, South`, col = 'blue')
lines(death_functions_fd$China, col = 'green')
lines(death_functions_fd$US, col = 'orange')
lines(death_functions_fd$Italy, col = 'purple')
lines(death_functions_fd$Germany, col = 'red')
lines(death_functions_fd$France, col = 'violet')
legend(1, 5000, legend=c("UK", "South Korea", "China", "US", "Italy", "Germany", "France"),
       col=c("black", "blue", "green", "orange", "purple", "red", "violet"), lty=1, cex=0.8)

plot.fd(death_functions_fd$`United Kingdom`, col = 'black', Lfdobj = 1, xlab = "Days from January 20th",
        ylab = "Value of derivative", ylim = c(0, 150), main = "Derivative of death curve")
lines(deriv.fd(death_functions_fd$`Korea, South`), col = 'blue')
lines(deriv.fd(death_functions_fd$China), col = 'green')
lines(deriv.fd(death_functions_fd$US), col = 'orange')
lines(deriv.fd(death_functions_fd$Italy), col = 'purple')
lines(deriv.fd(death_functions_fd$Germany), col = 'red')
lines(deriv.fd(death_functions_fd$France), col = 'violet')
legend(1, 150, legend=c("UK", "South Korea", "China", "US", "Italy", "Germany", "France"),
       col=c("black", "blue", "green", "orange", "purple", "red", "violet"), lty=1, cex=0.8)

plot.fd(death_functions_fd$`United Kingdom`, col = 'black', Lfdobj = 2, xlab = "Days from January 20th",
        ylab = "Value of second derivative", ylim = c(-20, 70), main = "Second derivative of death curve")
lines(deriv.fd(death_functions_fd$`Korea, South`, Lfdobj = 2), col = 'blue')
lines(deriv.fd(death_functions_fd$China, Lfdobj = 2), col = 'green')
lines(deriv.fd(death_functions_fd$US, Lfdobj = 2), col = 'orange')
lines(deriv.fd(death_functions_fd$Italy, Lfdobj = 2), col = 'purple')
lines(deriv.fd(death_functions_fd$Germany, Lfdobj = 2), col = 'red')
lines(deriv.fd(death_functions_fd$France, Lfdobj = 2), col = 'violet')
legend(1, 50, legend=c("UK", "South Korea", "China", "US", "Italy", "Germany", "France"),
       col=c("black", "blue", "green", "orange", "purple", "red", "violet"), lty=1, cex=0.8)
