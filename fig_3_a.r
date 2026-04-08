#--- init ----------------------------------------------------------------------

library(readxl)
library(annmatrix)
library(kklibrary)
library(matrixTests)
library(EpiDISH)

X <- readRDS("data/koncevicius/3_annmatrix.rds")


#--- clean ---------------------------------------------------------------------

# select cytosines
cgs <- c("cg05575921", "cg21566642", "cg01940273", "cg21161138")
X   <- X[cgs, X$celltype == "neu"]


#--- plot ----------------------------------------------------------------------

plot_trends <- function(x, y) {
  x <- x[,order(x$age), drop = FALSE] * 100
  y <- y[,order(y$age), drop = FALSE] * 100
  fitx <- loess(x[1,] ~ x$age)
  fity <- loess(y[1,] ~ y$age)

  plot(c(x$age, y$age), c(x,y), pch = 1, col = rep(c("gray20", "#ff3232"), c(length(x), length(y))), xlim = c(18,70), ylim = c(0,100), xlab = "Age (years)", ylab = "")
  lines(x$age, fitx$fitted, col = "gray20", lwd = 2)
  lines(y$age, fity$fitted, col = "#ff3232", lwd = 2)
  title(rownames(x), line = 0.5)
  title(ylab  = "Modification (%)", line = 1.6)
}

pdf("fig_3_a.pdf", width = 14/2.54, height = 4/2.54, pointsize = 8)

par(mfrow = c(1,4), mar = c(3,3,3,0), oma = c(0,0,0,1), tck = -0.02, mgp = c(2,0.5,0), las = 1)

plot_trends(X[1,!X$smoking,drop=FALSE], X[1,X$smoking,drop=FALSE])
mtext("a", 3, adj = -0.2, font = 2, line = 1.7)
plot_trends(X[2,!X$smoking,drop=FALSE], X[2,X$smoking,drop=FALSE])
plot_trends(X[3,!X$smoking,drop=FALSE], X[3,X$smoking,drop=FALSE])
plot_trends(X[4,!X$smoking,drop=FALSE], X[4,X$smoking,drop=FALSE])

legend("bottomright", c("Smoker", "Non-smoker"), text.width = strwidth(c("Smoker", "Non-smoker")), pch = 1, col = c("#ff3232", "gray20"), bty = 'n', inset = c(0, 1.1), xpd = NA, horiz = TRUE)


invisible(dev.off())
