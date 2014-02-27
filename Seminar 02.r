rnorm(n=10)
rnorm(10)
rnorm(10,100,30)
set.seed(540)
rnorm(10)
n <- 10
B <- 4
x <- matrix(rnorm(n * B), nrow = n)
str(x)
head(x)
x <- matrix(rnorm(n * B), nrow = n)
rownames(x) <- sprintf("obs%02d", 1:n)
colnames(x) <- sprintf("samp%02d", 1:B)
x
dim(x)
mean(x[,2])
colmean(x)
colMean(x)
colMeans(x)
apply(x,2,mean)
mean(x)
var(x)
sd(x)
B <- 1000
x10 <- matrix(rnorm(10 * B), nrow = 10)
x100 <- matrix(rnorm(100 * B), nrow = 100)
x1000 <- matrix(rnorm(1000 * B), nrow = 1000)
x10000 <- matrix(rnorm(10000 * B), nrow = 10000)
xBar10 <- colMeans(x10)
xBar100 <- colMeans(x100)
xBar1000 <- colMeans(x1000)
xBar10000 <- colMeans(x10000)
xBarSd10 <- sd(colMeans(x10))
xBarSd100 <- sd(colMeans(x100))
xBarSd1000 <- sd(colMeans(x1000))
xBarSd10000 <- sd(colMeans(x10000))
cbind(sampSize = c(10, 100, 1000, 10000),
      trueSEM = 1 / sqrt(c(10, 100, 1000, 10000)),
      obsSEM = c(xBarSd10, xBarSd100, xBarSd1000, xBarSd10000))
colMeans(x)
pnorm(10)
library(lattice)
n <- 35
x <- rnorm(n)
densityplot(~x)
densityplot(~x, n = 200, ylim = dnorm(0) * c(-0.1, 1.15),
            panel = function(x, ...) {
              panel.densityplot(x, ...)
              panel.mathdensity(n = 200, col.line = "grey74")
            })
B <- 1000
n <- round(10^(seq(from = 1, to = 2.5, length = 4)), 0)
names(n) <- paste0("n", n)
getSampleMeans <- function(n, B) colMeans(matrix(rnorm(n * B), nrow = n))
x <- data.frame(sapply(n, getSampleMeans, B))
jFormula <- as.formula(paste("~", paste(names(n), sep = "", collapse = " + ")))
densityplot(jFormula, x, xlab = "sample means",
            auto.key = list(x = 0.9, y = 0.9, corner = c(1, 1),
                            reverse.rows = TRUE))
xTallSkinny <- stack(x)
names(xTallSkinny) <- c("x","n")
xTallSkinny$n <- factor(xTallSkinny$n, levels = colnames(x))
densityplot(~ x, xTallSkinny, xlab = "sample means", groups = n,
            auto.key = list(x = 0.9, y = 0.9, corner = c(1, 1),
                            reverse.rows = TRUE))
