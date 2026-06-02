library(fossil)   # rand.index
library(funFEM)
library(fda)

# -----------------------
# robust funFEM runner (retry until no empty clusters)
run_funFEM_until_success <- function(fdobj, K,
                                     model = "AkjBk",
                                     lambda = 0,
                                     init_methods = c("kmeans", "random"),
                                     ntry = 50,
                                     min_cluster_size = 2,
                                     disp = FALSE) {
  
  for (trial in 1:ntry) {
    init_now <- init_methods[(trial - 1) %% length(init_methods) + 1]
    set.seed(sample.int(1e9, 1))
    
    fit <- tryCatch(
      suppressWarnings(
        funFEM(fdobj, K = K, model = model, init = init_now, lambda = lambda, disp = disp)
      ),
      error = function(e) NULL
    )
    
    if (!is.null(fit) && !is.null(fit$cls)) {
      tab <- table(fit$cls)
      if (length(tab) == K && all(tab >= min_cluster_size)) {
        return(fit)
      }
    }
  }
  NULL
}

# -----------------------
# simulation config
simul <- 100
times <- 50
endtime <- 5

Ktrue <- 4
objects <- 20

# true labels: 4 clusters, 5 each
trueclust <- rep(1:4, each = 5)

# target error variance
tau2 <- 0.125
sigma_e <- sqrt(tau2)   # Normal sd so Var = tau2

# smoothing basis
basis <- create.fourier.basis(c(1, times), nbasis = 7)

randtot <- rep(NA_real_, simul)

# -----------------------
for (iter in 1:simul) {
  
  timevec <- seq(0, endtime, length.out = times)
  
  # signal curves
  curves <- matrix(0, nrow = objects, ncol = times)
  
  curves[1:5, ]  <- -sin(timevec - 1) * log(timevec + 0.5)
  curves[6:10, ] <-  cos(timevec)     * log(timevec + 0.5)
  curves[11:15, ] <- -0.25 - 0.1 * cos(0.5 * (timevec - 1)) *
    sqrt(5 * sqrt(timevec) + 0.5) * (timevec^1.5)
  curves[16:20, ] <- -1 + 0.3 * timevec
  
  # -----------------------
  # iid Normal errors: N(0, tau2)
  eps <- matrix(rnorm(objects * times, mean = 0, sd = sigma_e),
                nrow = objects, ncol = times)
  
  obscurves <- curves + eps
  
  # fd object (smooth onto Fourier basis)
  fdobj <- smooth.basis(
    argvals = 1:times,
    y = t(obscurves),      # 50 x 20
    fdParobj = basis
  )$fd
  
  fem_fit <- run_funFEM_until_success(
    fdobj = fdobj,
    K = Ktrue,
    model = "AkjBk",
    lambda = 0,
    init_methods = c("kmeans", "random"),
    ntry = 60,
    min_cluster_size = 2,
    disp = FALSE
  )
  
  if (!is.null(fem_fit)) {
    randtot[iter] <- rand.index(trueclust, fem_fit$cls)
  } else {
    randtot[iter] <- NA_real_
  }
}

# -----------------------
# summary (exclude failures)
ok <- !is.na(randtot)
cat("funFEM success rate =", sum(ok), "/", simul, "\n")
cat("Normal noise: tau^2 =", tau2, " (sd =", sigma_e, ")\n")

mean_rand <- mean(randtot[ok])
sd_rand   <- sd(randtot[ok])
mcse      <- sd_rand / sqrt(sum(ok))
ci95      <- mean_rand + c(-1, 1) * 1.96 * mcse

cat("funFEM | K=4 | n=20 (5 each) | iid Normal(0,tau^2)\n")
cat("Mean Rand =", mean_rand, "\n")
cat("SD        =", sd_rand, "\n")
cat("MCSE      =", mcse, "\n")
cat("95% CI    = [", ci95[1], ", ", ci95[2], "]\n")