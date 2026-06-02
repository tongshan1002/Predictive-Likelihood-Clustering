library(cluster)
library(fossil)
library(splines)
library(fda.usc)

simul <- 100
times <- 50
endtime <- 5
dt <- endtime / times

sigma <- 0.5
drift <- 1

myk <- 4
objects <- 20

# true clusters: 4 groups, each 5
clussizes <- c(5, 10, 15, objects)
trueclust <- c(
  rep(1, 5),
  rep(2, 5),
  rep(3, 5),
  rep(4, 5)
)

randtot <- numeric(simul)

for (iter in 1:simul) {
  
  timevec <- seq(0, endtime, length.out = times)
  
  # OU noise
  normalvec <- matrix(rnorm(objects * times, sd = sigma),
                      nrow = objects, ncol = times)
  
  ouvalue <- matrix(0, nrow = objects, ncol = times)
  
  for (timeint in 2:times) {
    ouvalue[, timeint] <-
      -drift * ouvalue[, timeint - 1] * dt +
      sigma * normalvec[, timeint] * sqrt(dt)
  }
  
  # signal curves
  curves <- matrix(0, nrow = objects, ncol = times)
  
  curves[1:5, ]   <- -sin(timevec - 1) * log(timevec + 0.5)
  curves[6:10, ]  <-  cos(timevec) * log(timevec + 0.5)
  curves[11:15, ] <- -0.25 - 0.1 * cos(0.5 * (timevec - 1)) *
    sqrt(5 * sqrt(timevec) + 0.5) * (timevec^1.5)
  curves[16:20, ] <- -1 + 0.3 * timevec
  
  obscurves <- curves + ouvalue
  
  # 单次 kmeans
  VDP_fdata <- fdata(obscurves, argvals = timevec)
  
  fit <- suppressWarnings(
    kmeans.fd(VDP_fdata,
              ncl = myk,
              cluster.size = 1,
              draw = FALSE)
  )
  
  randtot[iter] <- rand.index(trueclust, fit$cluster)
}

# summary
mean_rand <- mean(randtot)
sd_rand   <- sd(randtot)
mcse      <- sd_rand / sqrt(simul)

cat("Single-start kmeans.fd\n")
cat("Mean Rand =", mean_rand, "\n")
cat("SD        =", sd_rand, "\n")
cat("MCSE      =", mcse, "\n")