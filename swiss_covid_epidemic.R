# Time is of the essence: containment of the SARS-CoV-2 epidemic in Switzerland from February to May 2020
# Christian L. Althaus, 13 July 2020

# Load libraries
library(deSolve)
library(RColorBrewer)
library(lubridate)
library(bbmle)
library(mvtnorm)
library(skellam)
library(xtable)

# Set random number generator
set.seed(613467)

# Load data
# Source: https://github.com/daenuprobst/covid19-cases-switzerland
# Retrieved: 15 May 2020
reported <- read.csv("swiss_covid_epidemic.csv")
reported$date <- ymd(reported$date)
firstcase <- ymd(20200225)
lockdown <- ymd(20200317)
relaxation <- ymd(20200510)

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
    ll <- sum(dpois(reported$deceased, diff(simulation$D), log = TRUE)) + sum(dskellam(diff(c(0, reported$hospitalized)), diff(simulation$H_in), diff(simulation$H_out), log = TRUE)) + sum(dskellam(diff(c(0, reported$icu)), diff(simulation$U_in), diff(simulation$U_out), log = TRUE))
    return(-ll)
}

# Parameter transformations
trans <- function(pars) {
    pars["R_0"] <- exp(pars["R_0"])
    pars["seed"] <- exp(pars["seed"])
    pars["kappa"] <- plogis(pars["kappa"])
    pars["epsilon3"] <- plogis(pars["epsilon3"])
    return(pars)
}

# Fit the model to incidence of deaths
popsize <- 8.6e6
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
           epsilon1 = 0.0075, # Probability of hospital admission (re-estimated to obtain an infection fatality ratio of 0.75%)
           control = as.numeric(lockdown - firstcase),
           slope = 0.5, # Slope of logistic function
           shift = 0, # Shift of logistic function
           epsilon2 = 0.3, # Probability of ICU admission when in hospital
           epsilon4 = 0.5) # Probability of death in ICU
free <- c(seed = log(25), # Time at which model is initialized
          R_0 = log(2.5), # Basic reproduction number
          kappa = qlogis(0.25), # Relative transmission after implementation of NPIs
          epsilon3 = qlogis(0.25)) # Probability of death in hospital bed
fit <- mle2(nll, start = as.list(free), fixed = as.list(fixed), method = "Nelder-Mead")
fit_ci <- confint(fit)
#saveRDS(fit, "fit.rds")
#saveRDS(fit_ci, "fit_ci.rds")
#fit <- readRDS("fit.rds")
#fit_ci <- readRDS("fit_ci.rds")

# Calculate bootstrap compatibility intervals
n_sim <- 1e4
m <- coef(fit, exclude.fixed = TRUE)
sigma <- vcov(fit)
sim_coef <- data.frame(rmvnorm(n_sim, mean = m, sigma = sigma))
sim_coef$seed <- exp(sim_coef$seed)
sim_coef$R_0 <- exp(sim_coef$R_0)
sim_coef$kappa <- plogis(sim_coef$kappa)
sim_coef$epsilon3 <- plogis(sim_coef$epsilon3)

est_seed <- c(trans(coef(fit))["seed"], exp(fit_ci[1, ]))
est_R_0 <- c(trans(coef(fit))["R_0"], exp(fit_ci[2, ]))
est_epsilon3 <- c(trans(coef(fit))["epsilon3"], plogis(fit_ci[3, ]))
est_kappa <- c(trans(coef(fit))["kappa"], plogis(fit_ci[4, ]))
est_R_e <- c(trans(coef(fit))["R_0"]*trans(coef(fit))["kappa"], quantile(sim_coef$R_0*sim_coef$kappa, probs = c(0.025, 0.975)))
names(est_R_e)[1] <- "R_e" 
est_epsilon1 <- c(fixed[["epsilon1"]]/(fixed[["epsilon2"]]*fixed[["epsilon4"]] + (1 - fixed[["epsilon2"]])*trans(coef(fit))[["epsilon3"]]),
                  quantile(fixed[["epsilon1"]]/(fixed[["epsilon2"]]*fixed[["epsilon4"]] + (1 - fixed[["epsilon2"]])*sim_coef$epsilon3), probs = c(0.025, 0.975)))

est_seed
est_R_0
log(2)/(fixed["gamma"]*(sqrt(est_R_0) - 1))
est_kappa
est_R_e
log(2)/(fixed["gamma"]*(sqrt(est_R_e) - 1))
c(- (sqrt(est_R_0[1]) - 1)/(sqrt(est_R_e[1]) - 1), quantile(- (sqrt(sim_coef$R_0) - 1)/(sqrt(sim_coef$R_0*sim_coef$kappa) - 1), probs = c(0.025, 0.975)))
est_epsilon1
est_epsilon3

# Simulate bootstrap samples
mod_length <- as.numeric(relaxation - firstcase) + 1
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
    mod_HU[, i] <- c(0, cumsum(rskellam(length(diff(sim$H_in)), diff(sim$H_in), diff(sim$H_out))))
    mod_U[, i] <- c(0, cumsum(rskellam(length(diff(sim$U_in)), diff(sim$U_in), diff(sim$U_out))))
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

# Plot best-fit model and compatibility intervals
parms <- trans(coef(fit))
times <- c(0, parms["seed"] + 0:(mod_length - 1))
sim <- as.data.frame(ode(inits, times, model, parms))
sim <- sim[-1, ]

timepoints1 <- firstcase + 0:(mod_length - 1)
timepoints2 <- firstcase + 1:(mod_length - 1)

# Compute some numbers
timepoints2[which(diff(sim$C1) == max(diff(sim$C1)))]
max(diff(sim$C1))
max(sim$C1)
max(sim$H1 + sim$H2 + sim$H3 + sim$U)
timepoints1[which((sim$H1 + sim$H2 + sim$H3 + sim$U) == max(sim$H1 + sim$H2 + sim$H3 + sim$U))]
max(sim$U)
timepoints1[which(sim$U == max(sim$U))]
timepoints2[which(diff(sim$D) == max(diff(sim$D)))]
max(diff(sim$D))

tt <- times - parms["seed"] - parms["control"] - parms["shift"]
R_e <- (1 - (1 - parms["kappa"])/(1 + exp(- parms["slope"]*tt)))*parms["R_0"]/popsize
R_e <- R_e[-1]*sim$S

cols <- brewer.pal(4, "Set1")
t.cols <- cols
for(i in 1:length(cols)) {
    x <- col2rgb(cols[i])
    t.cols[i] <- rgb(x[1, ], x[2, ], x[3, ], alpha = 125, maxColorValue = 255)
}

# Infection attack rate
tail(sim$C1, 1)/popsize
mod_C1[, 76]/popsize

# Figure 2: Model fit
par(mfrow = c(2, 2))
plot(timepoints2, diff(sim$C1),
     ylim = c(1, 2e4),
     log = "y",
     ty = "l",
     col = cols[4],
     xlab = NA, ylab = "Daily number of infections", main = "Daily infections", frame = FALSE)
polygon(x = c(timepoints2, rev(timepoints2)), y = c(mod_CI1[1,], rev(mod_CI1[2,])), col = t.cols[4], border = NA)
abline(v = lockdown, lty = 2)
points(reported$date, reported$confirmed, pch = 22, bg = "gray")
mtext("A", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)
legend("bottomright", inset = 0.05, c("By date of infection", "By date of reporting"), col = c(cols[4], "black"), lty = c(1, NA), pch = c(NA, 22), bty = "n")

plot(timepoints1, sim$H1 + sim$H2 + sim$H3 + sim$U,
     ylim = c(0, 3e3),
     ty = "l",
     col = cols[3],
     xlab = NA, ylab = "Number of hospitalized patients", main = "Hospitalized patients", frame = FALSE)
polygon(x = c(timepoints1, rev(timepoints1)), y = c(mod_HU[1,], rev(mod_HU[2,])), col = t.cols[3], border = NA)
abline(v = lockdown, lty = 2)
points(reported$date, reported$hospitalized, pch = 21, bg = "white")
mtext("B", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)

plot(timepoints1, sim$U,
     ylim = c(0, 5e2),
     ty = "l",
     col = cols[2],
     xlab = NA, ylab = "Number of patients in ICU", main = "Intensive care unit (ICU)", frame = FALSE)
polygon(x = c(timepoints1, rev(timepoints1)), y = c(mod_U[1,], rev(mod_U[2,])), col = t.cols[2], border = NA)
abline(v = lockdown, lty = 2)
points(reported$date, reported$icu, pch = 21, bg = "white")
mtext("C", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)

plot(timepoints2, diff(sim$D),
     ylim = c(0, 1e2),
     ty = "l",
     col = cols[1],
     xlab = NA, ylab = "Daily number of deaths", main = "Daily deaths", frame = FALSE)
polygon(x = c(timepoints2, rev(timepoints2)), y = c(mod_DI[1,], rev(mod_DI[2,])), col = t.cols[1], border = NA)
abline(v = lockdown, lty = 2)
points(reported$date, reported$deceased, pch = 21, bg = "white")
mtext("D", side = 3, line = 1.5, adj = - 0.1, cex = 1, font = 2)

# Figure 3: Effective reproduction number
par(mfrow = c(1, 1))
plot(timepoints1, R_e,
     ylim = c(0, 3),
     ty = "l",
     col = cols[4],
     xlab = NA, ylab = "Effective reproduction number", frame = FALSE)
polygon(x = c(timepoints1, rev(timepoints1)), y = c(mod_R_e[1,], rev(mod_R_e[2,])), col = t.cols[4], border = NA)
abline(v = lockdown, lty = 2)
abline(h = 1, lty = 2)

# Figure 4: Underreporting
par(mfrow = c(1, 1))
ma <- function(x, n = 7){filter(x, rep(1/n, n), sides = 2)}
plot(timepoints2, ma(reported$confirmed[-1]/diff(sim$C2)),
     ylim = c(0, 0.25),
     ty = "l",
     col = cols[4],
     xlab = NA, ylab = "Proportion confirmed cases", frame = FALSE)
x <- timepoints2[4:(length(timepoints2) - 3)]
ui <- (ma(reported$confirmed[-1]/mod_CI2[1,]))[4:(length(timepoints2) - 3)]
li <- (ma(reported$confirmed[-1]/mod_CI2[2,]))[4:(length(timepoints2) - 3)]
polygon(x = c(x, rev(x)), y = c(ui, rev(li)), col = t.cols[4], border = NA)
abline(v = lockdown, lty = 2)

# Simulate counterfactual scenarios
shift_range <- seq(- 7, 7, 1)
n_sim <- 1e3

parms <- trans(coef(fit))
times <- 0:1e3
mod_length <- length(times)

scenarios <- array(NA, dim = c(length(shift_range), 6, 3))

for(i in 1:length(shift_range)) {
    parms["shift"] <- shift_range[i]
    sim <- as.data.frame(ode(inits, times, model, parms))
    
    mod_C1 <- array(NA, n_sim)
    mod_CI1 <- array(NA, n_sim)
    mod_HU <- array(NA, n_sim)
    mod_U <- array(NA, n_sim)
    mod_D <- array(NA, n_sim)
    mod_DI <- array(NA, n_sim)
    
    for(j in 1:n_sim) {
        parms <- c(unlist(sim_coef[j, ]), fixed)
        parms["shift"] <- shift_range[i]
        sim <- as.data.frame(ode(inits, times, model, parms))
        
        x <- rpois(length(diff(sim$C1)), lambda = diff(sim$C1))
        mod_CI1[j] <- max(x)
        mod_C1[j] <- sum(x)
        mod_HU[j] <- max(cumsum(rskellam(length(diff(sim$H_in)), diff(sim$H_in), diff(sim$H_out))))
        mod_U[j] <- max(cumsum(rskellam(length(diff(sim$U_in)), diff(sim$U_in), diff(sim$U_out))))
        x <- rpois(length(diff(sim$D)), lambda = diff(sim$D))
        mod_DI[j] <- max(x)
        mod_D[j] <- sum(x)
    }
    scenarios[i, 1, ] <- quantile(mod_CI1, probs = c(0.025, 0.5, 0.975))
    scenarios[i, 2, ] <- quantile(mod_C1, probs = c(0.025, 0.5, 0.975))
    scenarios[i, 3, ] <- quantile(mod_HU, probs = c(0.025, 0.5, 0.975))
    scenarios[i, 4, ] <- quantile(mod_U, probs = c(0.025, 0.5, 0.975))
    scenarios[i, 5, ] <- quantile(mod_DI, probs = c(0.025, 0.5, 0.975))
    scenarios[i, 6, ] <- quantile(mod_D, probs = c(0.025, 0.5, 0.975))
}

# Table 2: Counterfactual scenarios
m <- array(NA, c(length(shift_range), 6))
for(i in 1:length(shift_range)) {
    for(j in 1:6) {
        m[i, j] <- paste0(
            format(round(scenarios[i, j, 2], 0), big.mark = ","),
            " (",
            format(round(scenarios[i, j, 1], 0), big.mark = ","),
            "--",
            format(round(scenarios[i, j, 3], 0), big.mark = ","),
            ")")
    }
}
m <- cbind(paste(shift_range, "days"), m)
m <- rbind(c("of NPIs", "daily infections", "", "hospitalized patients", "in ICU", "daily deaths", ""), m)
m <- rbind(c("Implementation", "Maximal number of", "Total infections", "Maximal number of", "Maximal number", "Maximal number of", "Total deaths"), m)
m <- xtable(m,
            caption = "Counterfactual scenarios for the SARS-CoV-2 epidemic in Switzerland from February to May 2020.",
            label = "tab:scenarios", align = "cccccccc")
m <- print(m, hline.after = c(0, 2, length(shift_range) + 2), size = "footnotesize", include.rownames = FALSE, include.colnames = FALSE)
write(m, file = "scenarios.tex")
