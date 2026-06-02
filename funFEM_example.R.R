library(fossil)   # rand.index
library(funFEM)
library(fda)

# 不设 seed（每次运行结果会变）

# -----------------------
# simulation config
simul <- 100
times <- 50
endtime <- 5
dt <- endtime / times

sigma <- 0.5
drift <- 1

Ktrue <- 2
objects <- 20

# true clusters: 2 groups, each 10
trueclust <- rep(1:2, each = 10)

# smoothing basis (Fourier)
basis <- create.fourier.basis(c(1, times), nbasis = 4)

randtot_fem <- numeric(simul)

# -----------------------
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
  
  # signal curves (ONLY type 1 and type 3)
  curves <- matrix(0, nrow = objects, ncol = times)
  
  # cluster 1: curves 1-10 (type 1)
  curves[1:10, ] <-
    -sin(timevec - 1) * log(timevec + 0.5)
  
  # cluster 2: curves 11-20 (type 3)
  curves[11:20, ] <-
    -0.25 - 0.1 * cos(0.5 * (timevec - 1)) *
    sqrt(5 * sqrt(timevec) + 0.5) * (timevec^1.5)
  
  # observed
  obscurves <- curves + ouvalue
  
  # -----------------------
  # funFEM
  fdobj <- smooth.basis(
    argvals = 1:times,
    y = t(obscurves),      # 50 x 20
    fdParobj = basis
  )$fd
  
  fem_fit <- suppressWarnings(
    funFEM(
      fdobj,
      K = Ktrue,
      model = "AkjBk",
      init = "kmeans",
      lambda = 0,
      disp = FALSE
    )
  )
  
  randtot_fem[iter] <- rand.index(trueclust, fem_fit$cls)
}

# -----------------------
# summary
mean_rand <- mean(randtot_fem)
sd_rand   <- sd(randtot_fem)
mcse      <- sd_rand / sqrt(simul)
ci95      <- mean_rand + c(-1,1) * 1.96 * mcse

cat("funFEM | 2 clusters (type1 vs type3) | 20 curves (10 each)\n")
cat("Mean Rand =", mean_rand, "\n")
cat("SD        =", sd_rand, "\n")
cat("MCSE      =", mcse, "\n")
cat("95% CI    = [", ci95[1], ", ", ci95[2], "]\n")