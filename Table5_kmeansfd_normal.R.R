library(fda.usc)
library(fossil)

# 不设 seed（每次运行结果会变）

simul <- 100

times <- 50
endtime <- 5

sigma_e <- 0.3   # NEW: iid normal noise sd

myk <- 4
objects <- 20

# true clusters: 4 groups, each 5
trueclust <- rep(1:4, each = 5)

randtot <- numeric(simul)

for (iter in 1:simul) {
  
  timevec <- seq(0, endtime, length.out = times)
  
  # -----------------------
  # Signal curves (4 groups, 5 each)
  curves <- matrix(0, nrow = objects, ncol = times)
  
  # cluster 1: 1-5
  curves[1:5, ] <-
    -sin(timevec - 1) * log(timevec + 0.5)
  
  # cluster 2: 6-10
  curves[6:10, ] <-
    cos(timevec) * log(timevec + 0.5)
  
  # cluster 3: 11-15
  curves[11:15, ] <-
    -0.25 - 0.1 * cos(0.5 * (timevec - 1)) *
    sqrt(5 * sqrt(timevec) + 0.5) * (timevec^1.5)
  
  # cluster 4: 16-20
  curves[16:20, ] <-
    -1 + 0.3 * timevec
  
  # -----------------------
  # NEW: iid Normal noise N(0, sigma_e^2)
  eps <- matrix(rnorm(objects * times, mean = 0, sd = sigma_e),
                nrow = objects, ncol = times)
  
  # observed curves
  obscurves <- curves + eps   # 20 x 50
  
  # -----------------------
  # fdata object: rows=curves, cols=time points
  VDP_fdata <- fdata(obscurves, argvals = timevec)
  
  # kmeans for functional data (single start)
  fit <- suppressWarnings(
    kmeans.fd(VDP_fdata, ncl = myk, cluster.size = 1)
  )
  
  # Rand index
  randtot[iter] <- rand.index(trueclust, fit$cluster)
}

# -----------------------
# Summary
mean_rand <- mean(randtot)
sd_rand   <- sd(randtot)
mcse      <- sd_rand / sqrt(simul)
ci95      <- mean_rand + c(-1, 1) * 1.96 * mcse

cat("kmeans.fd | K=4 | n=20 (5 each) | iid N(0,0.3^2) noise\n")
cat("Mean Rand =", mean_rand, "\n")
cat("SD        =", sd_rand, "\n")
cat("MCSE      =", mcse, "\n")
cat("95% CI    = [", ci95[1], ", ", ci95[2], "]\n")