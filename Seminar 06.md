Seminar 06
========================================================

```r
source("http://www.bioconductor.org/biocLite.R")
```

```
## Bioconductor version 2.13 (BiocInstaller 1.12.0), ?biocLite for help
```

```r
biocLite("limma")
```

```
## BioC_mirror: http://bioconductor.org
## Using Bioconductor version 2.13 (BiocInstaller 1.12.0), R version 3.0.2.
## Installing package(s) 'limma'
```

```
## package 'limma' successfully unpacked and MD5 sums checked
## 
## The downloaded binary packages are in
## 	C:\Users\Administrator\AppData\Local\Temp\Rtmp6Lp9UF\downloaded_packages
```

```
## Old packages: 'Hmisc', 'maptools', 'multcomp', 'plyr', 'testthat'
```

```r
biocLite("statmod")
```

```
## BioC_mirror: http://bioconductor.org
## Using Bioconductor version 2.13 (BiocInstaller 1.12.0), R version 3.0.2.
## Installing package(s) 'statmod'
```

```
## package 'statmod' successfully unpacked and MD5 sums checked
## 
## The downloaded binary packages are in
## 	C:\Users\Administrator\AppData\Local\Temp\Rtmp6Lp9UF\downloaded_packages
```

```
## Old packages: 'Hmisc', 'maptools', 'multcomp', 'plyr', 'testthat'
```


To load limma, lattice and ggplot2.

```r
library(limma)
library(lattice)
library(ggplot2)
```




```r
prDat <- read.table("GSE4051_data.tsv")
str(prDat, max.level = 0)
```

```
## 'data.frame':	29949 obs. of  39 variables:
```

```r
prDes <- readRDS("GSE4051_design.rds")
str(prDes)
```

```
## 'data.frame':	39 obs. of  4 variables:
##  $ sidChar : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
##  $ sidNum  : num  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
##  $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
```


DIY funtions

```r
prepareData <- function(myGenes) {
    miniDat <- t(wtDat[myGenes, ])
    miniDat <- data.frame(gExp = as.vector(miniDat), gene = factor(rep(colnames(miniDat), 
        each = nrow(miniDat)), levels = colnames(miniDat)))
    miniDat <- suppressWarnings(data.frame(wtDes, miniDat))
    miniDat
}
stripplotIt <- function(myData, ...) {
    stripplot(gExp ~ devStage | gene, myData, jitter.data = TRUE, auto.key = TRUE, 
        type = c("p", "a"), grid = TRUE, ...)
}
```


Simulation

```r
m <- 1000
n <- 3
x <- matrix(rnorm(m * n), nrow = m)
```



```r
obsVars <- apply(x, 1, var)
summary(obsVars)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   0.003   0.318   0.729   1.010   1.430   6.380
```

```r
mean(obsVars < 1/3)
```

```
## [1] 0.26
```

```r
densityplot(~obsVars, n = 200)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 


ANOVA

```r
wtDes <- subset(prDes, gType == "wt")
str(wtDes)
```

```
## 'data.frame':	20 obs. of  4 variables:
##  $ sidChar : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
##  $ sidNum  : num  20 21 22 23 24 25 26 27 28 29 ...
##  $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 2 2 2 2 3 3 ...
##  $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 1 1 1 1 1 1 ...
```

```r
wtDat <- subset(prDat, select = prDes$gType == "wt")
str(wtDat, max.level = 0)
```

```
## 'data.frame':	29949 obs. of  20 variables:
```




```r
wtDesMat <- model.matrix(~devStage, wtDes)
str(wtDesMat)
```

```
##  num [1:20, 1:5] 1 1 1 1 1 1 1 1 1 1 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:20] "12" "13" "14" "15" ...
##   ..$ : chr [1:5] "(Intercept)" "devStageP2" "devStageP6" "devStageP10" ...
##  - attr(*, "assign")= int [1:5] 0 1 1 1 1
##  - attr(*, "contrasts")=List of 1
##   ..$ devStage: chr "contr.treatment"
```




```r
wtFit <- lmFit(wtDat, wtDesMat)
wtEbFit <- eBayes(wtFit)
```




```r
topTable(wtEbFit)
```

```
##              X.Intercept. devStageP2 devStageP6 devStageP10
## 1423641_s_at        12.18    -0.0175     0.0750      0.0675
## 1438940_x_at        12.86     0.0850     0.1325      0.3425
## 1438657_x_at        12.78     0.1400     0.1250     -0.1850
## 1456736_x_at        12.32     0.1625     0.3050      0.2075
## 1436884_x_at        12.93     0.1775     0.3225      0.0300
## 1419700_a_at        12.32     0.1650     0.6475      0.8175
## 1435800_a_at        12.28     0.0450     0.6825      0.9000
## 1454613_at          12.47    -0.1075    -0.0500     -0.1025
## 1451240_a_at        13.00     0.3100     0.2800      0.2800
## 1450084_s_at        12.63     0.0825     0.0525      0.1725
##              devStage4_weeks AveExpr     F   P.Value adj.P.Val
## 1423641_s_at          0.1800   12.24 45350 3.574e-36 5.201e-32
## 1438940_x_at          0.3500   13.04 44957 3.865e-36 5.201e-32
## 1438657_x_at         -0.4500   12.71 43486 5.210e-36 5.201e-32
## 1456736_x_at          0.0725   12.47 39509 1.233e-35 6.725e-32
## 1436884_x_at          0.0250   13.04 39269 1.302e-35 6.725e-32
## 1419700_a_at          0.6825   12.78 39121 1.347e-35 6.725e-32
## 1435800_a_at          1.0200   12.81 36668 2.410e-35 1.031e-31
## 1454613_at           -0.3825   12.34 35835 2.962e-35 1.078e-31
## 1451240_a_at         -0.3700   13.10 35481 3.239e-35 1.078e-31
## 1450084_s_at          0.2600   12.75 34411 4.265e-35 1.234e-31
```

```r
topTable(wtEbFit, coef = 2:5)
```

```
##              devStageP2 devStageP6 devStageP10 devStage4_weeks AveExpr
## 1440645_at       0.3990     0.1953      0.9200           3.961   6.528
## 1416041_at       0.1580     0.4797      0.3327           5.114   9.383
## 1425222_x_at     0.8820     0.7995      1.5488           5.532   7.028
## 1451635_at       1.3025     1.1900      2.0160           6.188   8.319
## 1429028_at      -2.4433    -3.4073     -4.3105          -4.602   8.045
## 1422929_s_at    -2.9118    -3.6183     -3.5473          -3.661   7.278
## 1424852_at       0.4575     0.2298      0.5740           3.979   7.454
## 1425171_at       0.9980     3.0530      5.2788           6.079   9.620
## 1451617_at       0.7255     2.5128      4.9838           6.685   8.817
## 1451618_at       0.6028     2.8903      5.0508           6.288   9.431
##                  F   P.Value adj.P.Val
## 1440645_at   425.4 1.588e-17 4.755e-13
## 1416041_at   195.5 1.522e-14 2.280e-10
## 1425222_x_at 173.4 4.348e-14 4.341e-10
## 1451635_at   157.3 1.013e-13 7.585e-10
## 1429028_at   148.8 1.646e-13 9.203e-10
## 1422929_s_at 146.9 1.844e-13 9.203e-10
## 1424852_at   143.2 2.290e-13 9.799e-10
## 1425171_at   138.8 3.002e-13 1.124e-09
## 1451617_at   136.5 3.485e-13 1.160e-09
## 1451618_at   134.2 4.032e-13 1.207e-09
```

```r
colnames(coef(wtEbFit))
```

```
## [1] "(Intercept)"     "devStageP2"      "devStageP6"      "devStageP10"    
## [5] "devStage4_weeks"
```

```r
(dsHits <- topTable(wtEbFit, coef = grep("devStage", colnames(coef(wtEbFit)))))
```

```
##              devStageP2 devStageP6 devStageP10 devStage4_weeks AveExpr
## 1440645_at       0.3990     0.1953      0.9200           3.961   6.528
## 1416041_at       0.1580     0.4797      0.3327           5.114   9.383
## 1425222_x_at     0.8820     0.7995      1.5488           5.532   7.028
## 1451635_at       1.3025     1.1900      2.0160           6.188   8.319
## 1429028_at      -2.4433    -3.4073     -4.3105          -4.602   8.045
## 1422929_s_at    -2.9118    -3.6183     -3.5473          -3.661   7.278
## 1424852_at       0.4575     0.2298      0.5740           3.979   7.454
## 1425171_at       0.9980     3.0530      5.2788           6.079   9.620
## 1451617_at       0.7255     2.5128      4.9838           6.685   8.817
## 1451618_at       0.6028     2.8903      5.0508           6.288   9.431
##                  F   P.Value adj.P.Val
## 1440645_at   425.4 1.588e-17 4.755e-13
## 1416041_at   195.5 1.522e-14 2.280e-10
## 1425222_x_at 173.4 4.348e-14 4.341e-10
## 1451635_at   157.3 1.013e-13 7.585e-10
## 1429028_at   148.8 1.646e-13 9.203e-10
## 1422929_s_at 146.9 1.844e-13 9.203e-10
## 1424852_at   143.2 2.290e-13 9.799e-10
## 1425171_at   138.8 3.002e-13 1.124e-09
## 1451617_at   136.5 3.485e-13 1.160e-09
## 1451618_at   134.2 4.032e-13 1.207e-09
```




```r
stripplotIt(prepareData(rownames(dsHits)[c(3, 6, 9)]))
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11.png) 




```r
cutoff <- 1e-05
dsHits <- topTable(wtEbFit, coef = grep("devStage", colnames(coef(wtEbFit))), 
    p.value = cutoff, n = Inf)
(numBHhits <- nrow(dsHits))
```

```
## [1] 350
```




```r
dsHits[63, c("F", "adj.P.Val", "devStageP6")]
```

```
##                  F adj.P.Val devStageP6
## 1451633_a_at 64.01 1.049e-07      2.069
```



```r
P2Hits <- topTable(wtEbFit, coef = "devStageP2", n = Inf, sort = "none")
P10Hits <- topTable(wtEbFit, coef = "devStageP10", n = Inf, sort = "none")
xyplot(P10Hits$t ~ P2Hits$t, aspect = 1, xlab = "t-statistic for P2 effect", 
    ylab = "t-statistic for P10 effect", xlim = c(-20, 16), ylim = c(-20, 16), 
    panel = function(x, y, ...) {
        panel.smoothScatter(x, y, nbin = 100, ...)
        panel.abline(a = 0, b = 1, col = "orange")
    })
```

```
## KernSmooth 2.23 loaded
## Copyright M. P. Wand 1997-2009
## (loaded the KernSmooth namespace)
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14.png) 



```r
densityplot(~P10Hits$adj.P.Val + P2Hits$adj.P.Val, auto.key = TRUE, plot.points = FALSE, 
    n = 300)
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 



```r
cutoff <- 0.001
foo <- data.frame(P2 = P2Hits$adj.P.Val < cutoff, P10 = P10Hits$adj.P.Val < 
    cutoff)
addmargins(with(foo, table(P2, P10)))
```

```
##        P10
## P2      FALSE  TRUE   Sum
##   FALSE 29201   695 29896
##   TRUE      1    52    53
##   Sum   29202   747 29949
```



```r
P10pVals <- data.frame(raw = P10Hits$P.Value, BH = P10Hits$adj.P.Val, BY = p.adjust(P10Hits$P.Value, 
    method = "BY"))
splom(P10pVals, panel = function(x, y, ...) {
    panel.xyplot(x, y, pch = ".", ...)
    panel.abline(a = 0, b = 1, col = "orange")
})
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17.png) 



```r
colnames(wtDesMat)
```

```
## [1] "(Intercept)"     "devStageP2"      "devStageP6"      "devStageP10"    
## [5] "devStage4_weeks"
```

```r
(cont.matrix <- makeContrasts(P10VsP6 = devStageP10 - devStageP6, fourweeksVsP10 = devStage4_weeks - 
    devStageP10, levels = wtDesMat))
```

```
## Warning: Renaming (Intercept) to Intercept
```

```
##                  Contrasts
## Levels            P10VsP6 fourweeksVsP10
##   Intercept             0              0
##   devStageP2            0              0
##   devStageP6           -1              0
##   devStageP10           1             -1
##   devStage4_weeks       0              1
```



```r
wtFitCont <- contrasts.fit(wtFit, cont.matrix)
```

```
## Warning: row names of contrasts don't match col names of coefficients
```

```r
wtEbFitCont <- eBayes(wtFitCont)
```




```r
topTable(wtEbFitCont)
```

```
##              P10VsP6 fourweeksVsP10 AveExpr     F   P.Value adj.P.Val
## 1440645_at    0.7247          3.041   6.528 632.7 2.224e-17 6.662e-13
## 1416041_at   -0.1470          4.782   9.383 302.4 1.473e-14 2.206e-10
## 1425222_x_at  0.7493          3.983   7.028 235.4 1.300e-13 1.297e-09
## 1424852_at    0.3443          3.405   7.454 225.1 1.910e-13 1.430e-09
## 1420726_x_at  0.1733          3.551   7.190 203.5 4.555e-13 2.640e-09
## 1451635_at    0.8260          4.172   8.319 200.0 5.289e-13 2.640e-09
## 1429394_at   -0.0980          2.410   7.848 167.5 2.416e-12 1.034e-08
## 1455447_at   -0.9765         -1.800   9.973 153.5 5.063e-12 1.896e-08
## 1429791_at    0.2480          1.658   8.026 145.7 7.877e-12 2.621e-08
## 1422612_at    0.4838          3.426   8.833 142.2 9.676e-12 2.840e-08
```



```r
foo <- topTable(wtEbFitCont)
stripplotIt(prepareData(rownames(foo)[1:4]))
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21.png) 




```r
cutoff <- 1e-04
wtResCont <- decideTests(wtEbFitCont, p.value = cutoff, method = "global")
summary(wtResCont)
```

```
##    P10VsP6 fourweeksVsP10
## -1       4              8
## 0    29945          29895
## 1        0             46
```



```r
(hits1 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] < 0)])
```

```
## [1] "1416635_at" "1437781_at" "1454752_at" "1455260_at"
```

```r
stripplotIt(prepareData(hits1))
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23.png) 



```r
(hits2 <- rownames(prDat)[which(wtResCont[, "fourweeksVsP10"] < 0)])
```

```
## [1] "1416021_a_at" "1423851_a_at" "1434500_at"   "1435415_x_at"
## [5] "1437502_x_at" "1448182_a_at" "1452679_at"   "1455447_at"
```

```r
stripplotIt(prepareData(hits2[1:4]))
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-24.png) 





```r
intersect(hits1, hits2)
```

```
## character(0)
```





```r
(hits3 <- rownames(prDat)[which(wtResCont[, "fourweeksVsP10"] > 0)])
```

```
##  [1] "1416041_at"   "1417280_at"   "1418406_at"   "1418710_at"  
##  [5] "1418789_at"   "1419069_at"   "1420725_at"   "1420726_x_at"
##  [9] "1421061_at"   "1421818_at"   "1422612_at"   "1424852_at"  
## [13] "1424895_at"   "1425222_x_at" "1426059_at"   "1426223_at"  
## [17] "1427388_at"   "1428763_at"   "1429394_at"   "1429791_at"  
## [21] "1430580_at"   "1431174_at"   "1433699_at"   "1434297_at"  
## [25] "1434573_at"   "1435436_at"   "1435679_at"   "1435727_s_at"
## [29] "1436265_at"   "1436287_at"   "1440402_at"   "1440645_at"  
## [33] "1441518_at"   "1442243_at"   "1443252_at"   "1446484_at"  
## [37] "1449170_at"   "1449393_at"   "1451042_a_at" "1451635_at"  
## [41] "1452243_at"   "1453385_at"   "1455493_at"   "1457878_at"  
## [45] "1458418_at"   "1459904_at"
```

```r
stripplotIt(prepareData(hits3[1:4]))
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-26.png) 





```r
intersect(hits1, hits3)
```

```
## character(0)
```

```r
intersect(hits2, hits3)
```

```
## character(0)
```





```r
cutoff <- 0.01
nHits <- 8
wtResCont <- decideTests(wtEbFitCont, p.value = cutoff, method = "global")
summary(wtResCont)
```

```
##    P10VsP6 fourweeksVsP10
## -1      40             49
## 0    29897          29636
## 1       12            264
```

```r
hits1 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] < 0)]
stripplotIt(prepareData(hits1[1:nHits]))
```

![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-281.png) 

```r
hits2 <- rownames(prDat)[which(wtResCont[, "fourweeksVsP10"] < 0)]
stripplotIt(prepareData(hits2[1:nHits]))
```

![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-282.png) 

```r
hits3 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] > 0)]
stripplotIt(prepareData(hits3[1:nHits]))
```

![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-283.png) 

```r
hits4 <- rownames(prDat)[which(wtResCont[, "fourweeksVsP10"] > 0)]
stripplotIt(prepareData(hits4[1:nHits]))
```

![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-284.png) 

```r
vennDiagram(wtResCont)
```

![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-285.png) 

```r
hits5 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] != 0 & wtResCont[, "fourweeksVsP10"] != 
    0)]
stripplotIt(prepareData(hits5))
```

![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-286.png) 

```r
hits6 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] > 0 & wtResCont[, "fourweeksVsP10"] < 
    0)]
stripplotIt(prepareData(hits6))
```

![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-287.png) 





```r
lateStuff <- topTable(wtEbFitCont, n = Inf, sort = "none")
earlyStuff <- topTable(wtEbFit, coef = grep("devStageP[26]", colnames(coef(wtEbFit))), 
    n = Inf, sort = "none")
pVals <- data.frame(earlyStuff = earlyStuff$adj.P.Val, lateStuff = lateStuff$adj.P.Val)
xyplot(lateStuff ~ earlyStuff, pVals)
```

![plot of chunk unnamed-chunk-29](figure/unnamed-chunk-291.png) 

```r
discHits <- with(pVals, which(earlyStuff < quantile(earlyStuff, probs = 0.05) & 
    lateStuff > quantile(lateStuff, probs = 0.95)))
length(discHits)
```

```
## [1] 75
```

```r
set.seed(123)
stripplotIt(prepareData(miniDat <- sample(discHits, 6)), scales = list(y = list(relation = "free")))
```

![plot of chunk unnamed-chunk-29](figure/unnamed-chunk-292.png) 





