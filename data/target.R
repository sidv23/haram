

log.lik <- function(param, data) {

  delta <- param[1] # bounded real value
  c <- param[2] # bounded real value
  mu <- param[3] # bounded real value
  sigma2 <- param[4] # sigma squared (>0)
  tau <- param[5] # timescale (>0)

  time <- data[, 1]
  leng.time <- length(time)

  lcA <- data[, 2]
  se.lcA <- data[, 3]
  lcB <- data[, 4]
  se.lcB <- data[, 5]

  # sorting time given delta
  time.d <- time - delta
  time.temp <- c(time, time.d)
  ord <- order(time.temp)
  time.comb <- time.temp[ord]
  leng.time.comb <- length(time.comb)
  
  # indicator taking on 1 for X(t - delta) and 0 for X(t)
  ind <- c(rep(0, leng.time), rep(1, leng.time))[ord]

  lc.temp <- c(lcA, lcB - c)
  lc.comb <- lc.temp[ord]
  se.lc.temp <- c(se.lcA, se.lcB)
  se.lc.comb <- se.lc.temp[ord]
   
  # x.star.i, i = 1, 2, ..., 2n
  x <- lc.comb - mu

  # omega.i, i = 1, 2, ..., 2n to be saved
  B <- rep(NA, leng.time.comb)

  # x.hat.i, i = 1, 2, ..., 2n to be saved
  mu.i <- rep(NA, leng.time.comb)
  mu.star.i <- rep(NA, leng.time.comb)

  # a.i, i = 2, ..., 2n
  a.i <- exp( -diff(time.comb) / tau)

  # omega.i, i = 1, 2, ..., 2n
  var0 <- tau * sigma2 / 2
  B[1] <- se.lc.comb[1]^2 / (se.lc.comb[1]^2 + var0)

  for (k in 2 : leng.time.comb) {
    B[k] <- se.lc.comb[k]^2 / ( se.lc.comb[k]^2 + 
                                a.i[k - 1]^2 * (1 - B[k - 1]) * se.lc.comb[k - 1]^2 +
                                var0 * (1 - a.i[k - 1]^2) ) 
  }  

  # x.hat.i, i = 1, 2, ..., 2n
  mu.i[1] <- (1 - B[1]) * x[1]
  for (k in 2 : leng.time.comb) {
    mu.i[k] <- (1 - B[k]) * x[k] + B[k] * a.i[k - 1] * mu.i[k - 1]
  }

  mu.star.i[1] <- 0
  mu.star.i[2 : leng.time.comb] <- a.i * mu.i[-leng.time.comb]

  var.star.i <- se.lc.comb^2 / B

  # log-likelihood
  sum(dnorm(x, mean = mu.star.i, sd = sqrt(var.star.i), log = TRUE))

}


log.prior <- function(param) {

  delta <- param[1] # Unif(-1178.939, 1178.939)
  c <- param[2] # Unif(-60, 60)
  mu <- param[3] # Unif(-30, 30)
  sigma2 <- param[4] # sigma^2 ~ inverse-Gamma(1, 2 * 10^(-7))
  tau <- param[5] # inverse-Gamma(1, 1)

  if (delta < -1178.939 | delta > 1178.939 | c < -60 | c > 60 | mu < -30 | mu > 30) {
    -Inf
  } else {
    -2 * log(sigma2) - 2 * 10^(-7) / sigma2 - 2 *log(tau) - 1 / tau
  }
}

setwd("/Users/hyungsuktak/Desktop/")
dat <- read.csv("q0957usno.csv", header = TRUE)

param <- c(420, 0.6, c(-18, 0.01^2, 200))
log.lik(param, dat)
log.prior(param)

