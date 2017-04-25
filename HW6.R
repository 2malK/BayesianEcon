mu <- .1
sigma <- .2
alpha <- .05
N <- 111
y <- rnorm(N, mu, sigma)

mu.hat <- cumsum(y)[-1 * 1:9] / 10:N
sigma.hat <- sqrt((cumsum(y^2)[-1 * 1:9] - 2 * mu.hat * cumsum(y)[-1 * 1:9] + 
                     (10:N) * mu.hat^2) / (10:N - 1))

sd.test <- c()
for (i in 10:111) {
  sd.test <- c(sd.test, sd(y[1:i]))
}
all.equal(sd.test, sigma.hat)



HW6b <- function() {
  mu <- .1
  sigma <- .2
  alpha <- .05
  N <- 111
  y <- rnorm(N, mu, sigma)
  
  mu.hat <- cumsum(y)[-1 * 1:9] / 10:N
  sigma.hat <- sqrt((cumsum(y^2)[-1 * 1:9] - 2 * mu.hat * cumsum(y)[-1 * 1:9] + 
                       (10:N) * mu.hat^2) / (10:N - 1))
  
  exceed <- sum(qnorm(1 - alpha, mu.hat, sigma.hat)[-102] < y[-1 * 1:10]) / 101
  exceed
}

system.time(alpha.est <- replicate(n = 10000, HW6b()))

hist(alpha.est, freq = FALSE, breaks = 30)
mean(alpha.est)

## reason for exceedances: sqrt not a linear transformation!

##### c) #####
rinvgamma <- function(n, shape, rate) {
  1 / rgamma(n, shape, rate)
}
dinvgamma <- function(x, shape, rate) {
  dgamma(1/x, shape, rate) / x^2
}

qt_ls <- function(p, df, loc, scale) {
  qt(p, df) * scale + loc
}
## reasonable values for a,b?
a <- 3
b <- .08
r <- sqrt(rinvgamma(10, a, b))
r < .4 | r > .1

b/(a-1)
sigma^2

curve(dinvgamma(x, a, b))

HW6c <- function() {
  m <- .1
  M <- 1
  alpha <- .05
  
  a <- 3
  b <- .08
  sigma <- sqrt(rinvgamma(1, a, b))
  mu <- rnorm(1, m, sigma)
  N <- 111
  y <- rnorm(N, mu, sigma)
  
  mu.hat <- cumsum(y)[-1 * 1:9] / 10:N
  sigma.hat <- sqrt((cumsum(y^2)[-1 * 1:9] - 2 * mu.hat * cumsum(y)[-1 * 1:9] + 
                       (10:N) * mu.hat^2) / (10:N - 1))
  
  m.prime <- (m + M * Tt*mu.hat) / (M*Tt + 1)
  M.prime <- M / (M*Tt + 1)
  a.prime <- a + .5 * Tt
  b.prime <- b + sigma.hat^2 * (10:N - 1) * .5 +
              (Tt * (mu.hat - m)^2) / (2 * (M*Tt + 1))
  
  exceed <- sum(qt_ls(1-alpha, 2*a.prime, m.prime, (b.prime * (M.prime + 1)) / 
                        a.prime)[-102] < y[-1 * 1:10]) / 101
  exceed
}

system.time(alpha.est <- replicate(n = 10000, HW6c()))
hist(alpha.est, freq = FALSE, breaks = 30)
mean(alpha.est)

