
genedata <- read.table("C:/research/GENE/genemult.txt", header=FALSE)

lgenemult <- log(genedata)


meanlgenemult <- apply(lgenemult, 1, mean)
meanlgenemat <- matrix(meanlgenemult, nrow=78, ncol=18, byrow=FALSE)
lgenemultc <- lgenemult - meanlgenemat


starttime <- 7/60
endtime <- 7*18/60
timevec <- seq(starttime, endtime, length.out=18)


tranobscurves <- t(lgenemultc)


pdf("C:/research/GENE/GeneExpressionCurves.pdf", width=7, height=5)


plot(range(timevec), range(tranobscurves), type="n",
     xlab="Time", ylab="Centered log(expression)",
     main="Gene Expression Curves",
     cex.lab=1.4, cex.axis=1.2, cex.main=1.4, las=1)


for (i in 1:ncol(tranobscurves)) {
  lines(timevec, tranobscurves[,i], col=rainbow(78)[i], lwd=0.6)
}

dev.off()


genedata <- read.table("C:/research/GENE/genemult.txt", header = FALSE)
lgenemult <- log(genedata)


meanlgenemult <- apply(lgenemult, 1, mean)
meanlgenemat <- matrix(meanlgenemult, nrow = 78, ncol = 18, byrow = FALSE)
lgenemultc <- lgenemult - meanlgenemat


starttime <- 7/60
endtime <- 7*18/60
timevec <- seq(starttime, endtime, length.out = 18)


tranobscurves <- t(lgenemultc)


cluster_labels <- read.csv("C:/research/GENE/for1208bestgroup.csv", header = FALSE)$V1[-1]
cluster_labels <- as.numeric(cluster_labels)


cluster_colors <- c("#FF9999", "#99CCFF", "#99FF99", "#CC99FF", "#FFD699")


pdf("C:/research/GENE/Gene_Mean_5Clusters_pastel.pdf", width = 7, height = 5)


plot(range(timevec), range(tranobscurves), type = "n",
     xlab = "Time", ylab = "Centered log(expression)",
     main = "Mean Curves of 5 Gene Clusters",
     cex.lab = 1.4, cex.axis = 1.2, cex.main = 1.4, las = 1)

for (k in 1:5) {
  cluster_idx <- which(cluster_labels == k)
  mean_curve <- rowMeans(tranobscurves[, cluster_idx, drop = FALSE])
  lines(timevec, mean_curve, col = cluster_colors[k], lwd = 1.2)
}


legend("topright", legend = paste("Cluster", 1:5),
       col = cluster_colors, lty = 1, lwd = 1.5, bty = "n", cex = 1)

dev.off()
