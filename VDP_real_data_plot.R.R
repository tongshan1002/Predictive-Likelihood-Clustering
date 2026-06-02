
VDPdata <- read.table("C:/research/VDPdata.txt", header = TRUE)
X <- VDPdata[, 1]
Y <- VDPdata[, -1]
n <- ncol(Y)


set.seed(123)  
colors <- sample(rainbow(n), n)


pdf("C:/research/VDP_plot_colored.pdf", width = 7, height = 5)


plot(range(X), range(Y), type = "n",
     xlab = "Depth", ylab = "Board density",
     main = "Vertical Density Profiles",
     cex.lab = 1.4, cex.axis = 1.2, cex.main = 1.4, las = 1)


for (i in 1:n) {
  lines(X, Y[, i], col = colors[i], lwd = 0.8)
}

dev.off()


VDPdata <- read.table("C:/research/VDPdata.txt", header = TRUE)
X <- VDPdata[, 1]
Y <- VDPdata[, -1]


cluster_labels <- rep(NA, 24)
cluster_labels[c(1,2,3,4,5,8,14,15,20)] <- 1
cluster_labels[c(6,16)] <- 2
cluster_labels[c(9,11,17,18,24)] <- 3
cluster_labels[c(7,12,19)] <- 4
cluster_labels[10] <- 5
cluster_labels[13] <- 6
cluster_labels[21] <- 7


cluster_colors <- c("red", "blue", "forestgreen", "purple", "orange", "brown", "cyan")


pdf("C:/research/VDP_mean_curves_clean.pdf", width = 7, height = 5)


plot(range(X), range(Y), type = "n",
     xlab = "Depth", ylab = "Board density",
     main = "Mean Curves of 7 Clusters",
     cex.lab = 1.4, cex.axis = 1.2, cex.main = 1.4, las = 1)


for (k in 1:7) {
  idx <- which(cluster_labels == k)
  if (length(idx) > 0) {
    mean_curve <- rowMeans(Y[, idx, drop = FALSE])
    lines(X, mean_curve, col = cluster_colors[k], lwd = 0.8, lty = 1)
  }
}


legend("topright", legend = paste("Cluster", 1:7),
       col = cluster_colors, lty = 1, lwd = 1.5, bty = "n", cex = 1)

dev.off()