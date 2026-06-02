library(fda.usc)
library(fossil)

# 不设 seed（每次运行结果会变）

simul <- 100

times <- 50
endtime <- 5

myk <- 4
objects <- 20

# true clusters: 4 groups, each 5
trueclust <- rep(1:4, each = 5)

# target error variance
tau2 <- 0.125

# t distribution settings
nu <- 10
var_t <- nu / (nu - 2)               # = 1.25
scale_t <- sqrt(tau2 / var_t)        # scale so variance becomes tau2

randtot <- numeric(simul)

for (iter in 1:simul) {
  
  timevec <- seq(0, endtime, length.out = times)
  
  # -----------------------
  # Signal curves (4 groups, 5 each)
  curves <- matrix(0, nrow = objects, ncol = times)
  
  curves[1:5, ] <-
    -sin(timevec - 1) * log(timevec + 0.5)
  
  curves[6:10, ] <-
    cos(timevec) * log(timevec + 0.5)
  
  curves[11:15, ] <-
    -0.25 - 0.1 * cos(0.5 * (timevec - 1)) *
    sqrt(5 * sqrt(timevec) + 0.5) * (timevec^1.5)
  
  curves[16:20, ] <-
    -1 + 0.3 * timevec
  
  # -----------------------
  # NEW: scaled t noise with df=10 and Var = tau2
  eps <- scale_t * matrix(rt(objects * times, df = nu),
                          nrow = objects, ncol = times)
  
  # observed curves
  obscurves <- curves + eps
  
  # -----------------------
  # functional kmeans
  VDP_fdata <- fdata(obscurves, argvals = timevec)
  
  fit <- suppressWarnings(
    kmeans.fd(VDP_fdata, ncl = myk, cluster.size = 1)
  )
  
  randtot[iter] <- rand.index(trueclust, fit$cluster)
}

# -----------------------
# Summary
mean_rand <- mean(randtot)
sd_rand   <- sd(randtot)
mcse      <- sd_rand / sqrt(simul)
ci95      <- mean_rand + c(-1, 1) * 1.96 * mcse

cat("kmeans.fd | K=4 | n=20 (5 each) | scaled t(df=10) noise | Var=tau^2\n")
cat("tau^2     =", tau2, "\n")
cat("scale_t   =", scale_t, "\n")
cat("Mean Rand =", mean_rand, "\n")
cat("SD        =", sd_rand, "\n")
cat("MCSE      =", mcse, "\n")
cat("95% CI    = [", ci95[1], ", ", ci95[2], "]\n")