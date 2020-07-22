# Modeling the SARS-CoV-2 epidemic in the canton of Zurich, Switzerland
# Christian L. Althaus, 21 July 2020

# Load libraries
library(deSolve)
library(RColorBrewer)
library(lubridate)
library(bbmle)
library(mvtnorm)
library(plotrix)

# Set random number generator
set.seed(174673)

# Load data
# Source: https://github.com/openZH/covid_19
# Retrieved: 18 July 2020
reported <- read.csv("swiss_covid_epidemic_ZH.csv")
reported$date <- ymd(reported$date)
firstcase <- ymd(20200227)
lockdown <- ymd(20200317)
lastpoint <- ymd(20200531)
reported <- subset(reported, date <= lastpoint)

# Seroprevalence data from Emmenegger et al. (2020, medRxiv)
date <- c(ymd(20200315), ymd(20200415), ymd(20200515))
sample <- c(1098, 1469, 1640)
pos <- c(1, 17, 26)
prev <- pos/sample
lower <- numeric()
upper <- numeric()
for(i in 1:length(date)) {
    x <- binom.test(pos[i], sample[i])
    lower[i] <- x$conf.int[1]
    upper[i] <- x$conf.int[2]
}
sero <- data.frame(date, sample, pos, prev, lower, upper)

# SARS-CoV-2 transmission model
model <- function(t, x, parms) {
	with(as.list(c(parms, x)), {
	    tt <- t - seed - control - shift
	    beta <- (1 - (1 - kappa)/(1 + exp(- slope*tt)))*R_0*gamma/popsize
	    epsilon1 <- epsilon1/(epsilon2*epsilon4 + (1 - epsilon2)*epsilon3)
	    dS <- - beta*S*I
	    dE <- beta*S*I - sigma*E
	    dI <- sigma*E - gamma*I
	    dP <- gamma*I - omega1*P
	    dH1 <- epsilon1*omega1*P - omega2*H1
	    dH2 <- (1 - epsilon2)*omega2*H1 - omega3*H2
	    dH3 <- (1 - epsilon3)*omega3*H2 + (1 - epsilon4)*omega3*U - omega4*H3
	    dU <- epsilon2*omega2*H1 - omega3*U
        dR <- (1 - epsilon1)*omega1*P + omega4*H3
        dD <- epsilon3*omega3*H2 + epsilon4*omega3*U
	    dC1 <- beta*S*I
	    dC2 <- omega1*P
		dH_in <- epsilon1*omega1*P
		dH_out <- omega4*H3 + epsilon3*omega3*H2 + epsilon4*omega3*U
		dU_in <- epsilon2*omega2*H1
		dU_out <- omega3*U
		der <- c(dS, dE, dI, dP, dH1, dH2, dH3, dU, dR, dD, dC1, dC2, dH_in, dH_out, dU_in, dU_out)
		list(der)
	})
}

# Negative log-likelihood
nll <- function(seed, R_0, sigma, gamma, omega1, omega2, omega3, omega4, epsilon1, epsilon2, epsilon3, epsilon4, control, kappa, slope, shift) {
    pars <- c(seed = seed, R_0 = R_0, sigma = sigma, gamma = gamma, omega1 = omega1, omega2 = omega2, omega3 = omega3, omega4 = omega4, epsilon1 = epsilon1, epsilon2 = epsilon2, epsilon3 = epsilon3, epsilon4 = epsilon4, control = control, kappa = kappa, slope = slope, shift = shift)
    pars <- trans(pars)
    times <- c(0, pars["seed"] + reported$date - min(reported$date))
    times <- c(times, max(times + 1))
    simulation <- as.data.frame(ode(inits, times, model, pars))
    simulation <- simulation[-1, ]
    w_H <- which(is.na(diff(c(0, reported$hospitalized))))
    w_HT <- which(is.na(reported$hospitalized))
    w_U <- which(is.na(diff(c(0, reported$icu))))
    w_UT <- which(is.na(reported$icu))
    ll <- sum(dpois(reported$deceased, diff(simulation$D), log = TRUE)) + 
        sum(dnbinom(reported$hospitalized[- w_HT], size = 0.5*(simulation$H1 + simulation$H2 + simulation$H3 + simulation$U)[- w_HT], mu = (simulation$H1 + simulation$H2 + simulation$H3 + simulation$U)[- w_HT], log = TRUE)) + 
        sum(dnbinom(reported$icu[- w_UT], size = 0.5*simulation$U[- w_UT], mu = simulation$U[- w_UT], log = TRUE)) +
        sum(dbinom(sero$pos, sero$sample, simulation$C2[sero$date - firstcase]/popsize, log = TRUE))
    return(-ll)
}

# Parameter transformations
trans <- function(pars) {
    pars["R_0"] <- exp(pars["R_0"])
    pars["seed"] <- exp(pars["seed"])
    pars["kappa"] <- plogis(pars["kappa"])
    pars["epsilon1"] <- plogis(pars["epsilon1"])
    pars["epsilon2"] <- plogis(pars["epsilon2"])
    pars["epsilon3"] <- plogis(pars["epsilon3"])
    return(pars)
}

# Fit the model to incidence of deaths
popsize <- 1536000
times <- 0:200			
inits <- c(S = popsize - 1,
           E = 0,
           I = 1,
           P = 0,
           H1 = 0,
           H2 = 0,
           H3 = 0,
           U = 0,
           R = 0,
           D = 0,
           C1 = 1,
           C2 = 0,
           H_in = 0,
           H_out = 0,
           U_in = 0,
           U_out = 0)
fixed <- c(sigma = 1/2.6, # Incubation period
           gamma = 1/2.6, # Infectious period
           omega1 = 1/6, # Duration from symptom onset to hospitalization/recovery
           omega2 = 1/2, # Duration from hospitalization to ICU admission
           omega3 = 1/8, # Additional duration in ICU/hospital
           omega4 = 1/5, # Duration of recovery in hospital
           control = as.numeric(lockdown - firstcase),
           slope = 0.5, # Slope of logistic function
           shift = 0, # Shift of logistic function
           epsilon4 = 0.23) # Probability of death in ICU
free <- c(seed = log(25), # Time at which model is initialized
          R_0 = log(2), # Basic reproduction number
          kappa = qlogis(0.30), # Relative transmission after implementation of NPIs
          epsilon1 = qlogis(0.006), # Probability of hospital admission (re-estimated to obtain an infection fatality ratio of 0.75%)
          epsilon2 = qlogis(0.5), # Probability of ICU admission when in hospital
          epsilon3 = qlogis(0.23)) # Probability of death in hospital bed
fit <- mle2(nll, start = as.list(free), fixed = as.list(fixed), method = "Nelder-Mead")
fit_ci <- confint(fit)
#saveRDS(fit, "fit.rds")
#saveRDS(fit_ci, "fit_ci.rds")
#fit <- readRDS("fit.rds")
#fit_ci <- readRDS("fit_ci.rds")

# Calculate bootstrap compatibility intervals
n_sim <- 1e3
m <- coef(fit, exclude.fixed = TRUE)
sigma <- vcov(fit)
sim_coef <- data.frame(rmvnorm(n_sim, mean = m, sigma = sigma))
sim_coef$seed <- exp(sim_coef$seed)
sim_coef$R_0 <- exp(sim_coef$R_0)
sim_coef$kappa <- plogis(sim_coef$kappa)
sim_coef$epsilon1 <- plogis(sim_coef$epsilon1)
sim_coef$epsilon2 <- plogis(sim_coef$epsilon2)
sim_coef$epsilon3 <- plogis(sim_coef$epsilon3)

est_seed <- c(trans(coef(fit))["seed"], exp(fit_ci[1, ]))
est_R_0 <- c(trans(coef(fit))["R_0"], exp(fit_ci[2, ]))
est_kappa <- c(trans(coef(fit))["kappa"], plogis(fit_ci[6, ]))
est_R_e <- c(trans(coef(fit))["R_0"]*trans(coef(fit))["kappa"], quantile(sim_coef$R_0*sim_coef$kappa, probs = c(0.025, 0.975)))
names(est_R_e)[1] <- "R_e" 
est_epsilon1 <- c(trans(coef(fit))["epsilon1"], plogis(fit_ci[3, ]))
est_epsilon2 <- c(trans(coef(fit))["epsilon2"], plogis(fit_ci[4, ]))
est_epsilon3 <- c(trans(coef(fit))["epsilon3"], plogis(fit_ci[5, ]))

est_seed
est_R_0
est_kappa
est_R_e
est_epsilon1
est_epsilon2
est_epsilon3

# Simulate bootstrap samples for prediction intervals
mod_length <- as.numeric(lastpoint - firstcase) + 1
mod_C1 <- matrix(NA, nrow = mod_length, ncol = n_sim)
mod_CI1 <- matrix(NA, nrow = mod_length - 1, ncol = n_sim)
mod_C2 <- matrix(NA, nrow = mod_length, ncol = n_sim)
mod_CI2 <- matrix(NA, nrow = mod_length - 1, ncol = n_sim)
mod_HU <- matrix(NA, nrow = mod_length, ncol = n_sim)
mod_U <- matrix(NA, nrow = mod_length, ncol = n_sim)
mod_D <- matrix(NA, nrow = mod_length, ncol = n_sim)
mod_DI <- matrix(NA, nrow = mod_length - 1, ncol = n_sim)
mod_R_e <- matrix(NA, nrow = mod_length, ncol = n_sim)

for(i in 1:n_sim) {
    parms <- c(unlist(sim_coef[i, ]), fixed)
    times <- c(0, parms["seed"] + 0:(mod_length - 1))
    sim <- as.data.frame(ode(inits, times, model, parms))
    sim <- sim[-1, ]
    mod_CI1[, i] <- rpois(length(diff(sim$C1)), lambda = diff(sim$C1))
    mod_C1[, i] <- c(0, cumsum(mod_CI1[, i]))
    mod_CI2[, i] <- rpois(length(diff(sim$C2)), lambda = diff(sim$C2))
    mod_C2[, i] <- c(0, cumsum(mod_CI2[, i]))
    mod_HU[, i] <- rnbinom(length(sim$H1), size = 0.5*(sim$H1 + sim$H2 + sim$H3 + sim$U), mu = (sim$H1 + sim$H2 + sim$H3 + sim$U))
    mod_U[, i] <- rnbinom(length(sim$U), size = 0.5*sim$U, mu = sim$U)
    mod_DI[, i] <- rpois(length(diff(sim$D)), lambda = diff(sim$D))
    mod_D[, i] <- c(0, cumsum(mod_DI[, i]))
    tt <- times - parms["seed"] - parms["control"] - parms["shift"]
    mod_R_e[, i] <- ((1 - (1 - parms["kappa"])/(1 + exp(- parms["slope"]*tt)))*parms["R_0"]/popsize)[-1]*sim$S
}

mod_C1 <- apply(mod_C1, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))
mod_CI1 <- apply(mod_CI1, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))
mod_C2 <- apply(mod_C2, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))
mod_CI2 <- apply(mod_CI2, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))
mod_HU <- apply(mod_HU, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))
mod_U <- apply(mod_U, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))
mod_D <- apply(mod_D, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))
mod_DI <- apply(mod_DI, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))
mod_R_e <- apply(mod_R_e, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))

# Calculate seroprevalence
tail(sim$C2, 1)/popsize
mod_C2[, length(sim$C2)]/popsize

# Plot best-fit model and predition intervals
parms <- trans(coef(fit))
times <- c(0, parms["seed"] + 0:(mod_length - 1))
sim <- as.data.frame(ode(inits, times, model, parms))
sim <- sim[-1, ]

timepoints1 <- firstcase + 0:(mod_length - 1)
timepoints2 <- firstcase + 1:(mod_length - 1)

# Set colors
cols <- brewer.pal(4, "Set1")
t.cols <- cols
for(i in 1:length(cols)) {
    x <- col2rgb(cols[i])
    t.cols[i] <- rgb(x[1, ], x[2, ], x[3, ], alpha = 125, maxColorValue = 255)
}

# Figures
par(mfrow = c(1, 1))
plot(timepoints1, sim$C2/popsize,
     ylim = c(0, 0.025),
     ty = "l",
     col = cols[4],
     xlab = NA, ylab = "Seroprevalence", frame = FALSE)
polygon(x = c(timepoints1, rev(timepoints1)), y = c(mod_C2[1,]/popsize, rev(mod_C2[2,]/popsize)), col = t.cols[4], border = NA)
plotCI(sero$date, sero$prev, ui = sero$upper, li = sero$lower, pch = 21, add = TRUE)

par(mfrow = c(1, 3))
plot(timepoints1, sim$H1 + sim$H2 + sim$H3 + sim$U,
     ylim = c(0, 3e2),
     ty = "l",
     col = cols[3],
     xlab = NA, ylab = "Number of hospitalized patients", main = "Hospitalized patients", frame = FALSE)
polygon(x = c(timepoints1, rev(timepoints1)), y = c(mod_HU[1,], rev(mod_HU[2,])), col = t.cols[3], border = NA)
points(reported$date, reported$hospitalized, pch = 21, bg = "white")
mtext("A", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)

plot(timepoints1, sim$U,
     ylim = c(0, 1e2),
     ty = "l",
     col = cols[2],
     xlab = NA, ylab = "Number of ventilated patients", main = "Ventilated patients", frame = FALSE)
polygon(x = c(timepoints1, rev(timepoints1)), y = c(mod_U[1,], rev(mod_U[2,])), col = t.cols[2], border = NA)
points(reported$date, reported$icu, pch = 21, bg = "white")
mtext("B", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)

plot(timepoints2, diff(sim$D),
     ylim = c(0, 1e1),
     ty = "l",
     col = cols[1],
     xlab = NA, ylab = "Daily number of deaths", main = "Daily deaths", frame = FALSE)
polygon(x = c(timepoints2, rev(timepoints2)), y = c(mod_DI[1,], rev(mod_DI[2,])), col = t.cols[1], border = NA)
points(reported$date, reported$deceased, pch = 21, bg = "white")
mtext("C", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)
