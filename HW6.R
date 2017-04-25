##### a) b) #####

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

dt_ls <- function(x, df, loc, scale) {
  1/scale * dt((x - loc)/scale, df)
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
  sigma.true <- .2
  
  N <- 111
  y <- rnorm(N, m, sigma.true)
  alpha <- .05
  
  Tt <- 10:111
  a <- 3
  b <- .08
  
  sigma <- sqrt(rinvgamma(1, a, b))
  mu <- rnorm(1, m, sigma)
  
  mu.hat <- cumsum(y)[-1 * 1:9] / 10:N
  sigma.hat <- sqrt((cumsum(y^2)[-1 * 1:9] - 2 * mu.hat * cumsum(y)[-1 * 1:9] + 
                       (10:N) * mu.hat^2) / (10:N - 1))
  
  m.prime <- (m + M * Tt*mu.hat) / (M*Tt + 1)
  M.prime <- M / (M*Tt + 1)
  a.prime <- a + .5 * Tt
  b.prime <- b + sigma.hat^2 * (10:N - 1) * .5 +
              (Tt * (mu.hat - m)^2) / (2 * (M*Tt + 1))
  
  exceed <- sum(qt_ls(1-alpha, 2*a.prime, m.prime, sqrt((b.prime * (M.prime + 1)) / 
                        a.prime))[-102] < y[-1 * 1:10]) / 101
  exceed
}

system.time(alpha.est <- replicate(n = 10000, HW6c()))
hist(alpha.est, freq = FALSE, breaks = 30)
mean(alpha.est)

##### d) #####
K <- 10000 # call M K

m <- .1
M <- 1
sigma.true <- .2

N <- 111
y <- rnorm(N, m, sigma.true)
alpha <- .05

Tt <- 10:111
a <- 3
b <- .08

mu.hat <- cumsum(y)[-1 * 1:9] / 10:N
sigma.hat <- sqrt((cumsum(y^2)[-1 * 1:9] - 2 * mu.hat * cumsum(y)[-1 * 1:9] + 
                     (10:N) * mu.hat^2) / (10:N - 1))

#### from slide 22, we know the cond. posteriors
## sigma does not depend on mu -> direct sampling possible

HW6d <- function() {
  shape.sigma <- a + .5 * Tt
  rate.sigma <- b + sigma.hat^2 * (10:N - 1) * .5 +
    (Tt * (mu.hat - m)^2) / (2 * (M*Tt + 1))
  
  sigmasq.m <- rinvgamma(length(shape.sigma), shape.sigma, rate.sigma)
  
  loc.mu <- (m + M*Tt*mu.hat) / (M*Tt + 1)
  scale.mu <- (M * sigmasq.m) / (M*Tt + 1)
  
  mu.m <- rnorm(length(sigma.m), loc.mu, sqrt(scale.mu))
  
  ## we return 1 simulation of y_t for each T in 10:111
  rnorm(length(mu.m), mu.m, sqrt(sigmasq.m))
}

system.time(dist.out <- replicate(K, HW6d()))
rownames(dist.out) <- paste0('T', 10:111)
colnames(dist.out) <- paste0('m', 1:K)

#### true density ####
m.prime <- (m + M * Tt*mu.hat) / (M*Tt + 1)
M.prime <- M / (M*Tt + 1)
a.prime <- a + .5 * Tt
b.prime <- b + sigma.hat^2 * (10:N - 1) * .5 +
  (Tt * (mu.hat - m)^2) / (2 * (M*Tt + 1))
######################

ind <- match(paste0('T', c(10, 50, 100)), rownames(dist.out))

for(i in ind) {
  op <- par(cex = .8) #, mfcol = c(3,1))
  plot(density(dist.out[i,]), xlab = '', main = rownames(dist.out)[i], lwd = 2)
  legend('topright', c('Simulation', 'True Density'), 
         lty = 1, cex = .8, col = c("black", "dodgerblue3"), bty = 'n')
  curve(dt_ls(x, 2*a.prime[i], m.prime[i], sqrt((b.prime[i] * (M.prime[i] + 1)) / 
                                            a.prime[i])), add = TRUE, col = 'dodgerblue3')
  par(op)
}
