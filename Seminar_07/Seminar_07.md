Seminar 07

========================================================
Optional material: from BAM file to count data
BAM/SAM aligned sequence data file format


```r
# source ('http://bioconductor.org/biocLite.R') biocLite('ShortRead')
# biocLite('Rsamtools') biocLite('easyRNASeq')
# biocLite('BSgenome.Dmelanogaster.UCSC.dm3') biocLite('biomaRt')
```



```r
library(knitr)
library(ShortRead)
```

```
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist
## 
## Loading required package: IRanges
## Loading required package: GenomicRanges
## Loading required package: XVector
## Loading required package: Biostrings
## Loading required package: lattice
## Loading required package: Rsamtools
```

```r
library(Rsamtools)
bamDat <- readAligned("drosophilaMelanogasterSubset.bam", type = "BAM")
str(bamDat)
```

```
## Formal class 'AlignedRead' [package "ShortRead"] with 8 slots
##   ..@ chromosome  : Factor w/ 15 levels "chrYHet","chrM",..: 2 2 2 2 2 2 2 2 2 2 ...
##   ..@ position    : int [1:64206] 548 1497 1506 1528 1540 1552 1552 1555 1559 1566 ...
##   ..@ strand      : Factor w/ 3 levels "+","-","*": 2 1 1 1 1 1 1 1 2 2 ...
##   ..@ alignQuality:Formal class 'NumericQuality' [package "ShortRead"] with 1 slots
##   .. .. ..@ quality: int [1:64206] 132 132 127 130 130 122 132 132 132 132 ...
##   ..@ alignData   :Formal class 'AlignedDataFrame' [package "ShortRead"] with 4 slots
##   .. .. ..@ varMetadata      :'data.frame':	1 obs. of  1 variable:
##   .. .. .. ..$ labelDescription: chr "Type of read; see ?scanBam"
##   .. .. ..@ data             :'data.frame':	64206 obs. of  1 variable:
##   .. .. .. ..$ flag: int [1:64206] 16 0 0 0 0 0 0 0 16 16 ...
##   .. .. ..@ dimLabels        : chr [1:2] "readName" "alignColumn"
##   .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
##   .. .. .. .. ..@ .Data:List of 1
##   .. .. .. .. .. ..$ : int [1:3] 1 1 0
##   ..@ quality     :Formal class 'FastqQuality' [package "ShortRead"] with 1 slots
##   .. .. ..@ quality:Formal class 'BStringSet' [package "Biostrings"] with 5 slots
##   .. .. .. .. ..@ pool           :Formal class 'SharedRaw_Pool' [package "XVector"] with 2 slots
##   .. .. .. .. .. .. ..@ xp_list                    :List of 1
##   .. .. .. .. .. .. .. ..$ :<externalptr> 
##   .. .. .. .. .. .. ..@ .link_to_cached_object_list:List of 1
##   .. .. .. .. .. .. .. ..$ :<environment: 0x04dd3e38> 
##   .. .. .. .. ..@ ranges         :Formal class 'GroupedIRanges' [package "XVector"] with 7 slots
##   .. .. .. .. .. .. ..@ group          : int [1:64206] 1 1 1 1 1 1 1 1 1 1 ...
##   .. .. .. .. .. .. ..@ start          : int [1:64206] 1 37 73 109 145 181 217 253 289 325 ...
##   .. .. .. .. .. .. ..@ width          : int [1:64206] 36 36 36 36 36 36 36 36 36 36 ...
##   .. .. .. .. .. .. ..@ NAMES          : NULL
##   .. .. .. .. .. .. ..@ elementType    : chr "integer"
##   .. .. .. .. .. .. ..@ elementMetadata: NULL
##   .. .. .. .. .. .. ..@ metadata       : list()
##   .. .. .. .. ..@ elementType    : chr "BString"
##   .. .. .. .. ..@ elementMetadata: NULL
##   .. .. .. .. ..@ metadata       : list()
##   ..@ sread       :Formal class 'DNAStringSet' [package "Biostrings"] with 5 slots
##   .. .. ..@ pool           :Formal class 'SharedRaw_Pool' [package "XVector"] with 2 slots
##   .. .. .. .. ..@ xp_list                    :List of 1
##   .. .. .. .. .. ..$ :<externalptr> 
##   .. .. .. .. ..@ .link_to_cached_object_list:List of 1
##   .. .. .. .. .. ..$ :<environment: 0x04dd3e38> 
##   .. .. ..@ ranges         :Formal class 'GroupedIRanges' [package "XVector"] with 7 slots
##   .. .. .. .. ..@ group          : int [1:64206] 1 1 1 1 1 1 1 1 1 1 ...
##   .. .. .. .. ..@ start          : int [1:64206] 1 37 73 109 145 181 217 253 289 325 ...
##   .. .. .. .. ..@ width          : int [1:64206] 36 36 36 36 36 36 36 36 36 36 ...
##   .. .. .. .. ..@ NAMES          : NULL
##   .. .. .. .. ..@ elementType    : chr "integer"
##   .. .. .. .. ..@ elementMetadata: NULL
##   .. .. .. .. ..@ metadata       : list()
##   .. .. ..@ elementType    : chr "DNAString"
##   .. .. ..@ elementMetadata: NULL
##   .. .. ..@ metadata       : list()
##   ..@ id          :Formal class 'BStringSet' [package "Biostrings"] with 5 slots
##   .. .. ..@ pool           :Formal class 'SharedRaw_Pool' [package "XVector"] with 2 slots
##   .. .. .. .. ..@ xp_list                    :List of 1
##   .. .. .. .. .. ..$ :<externalptr> 
##   .. .. .. .. ..@ .link_to_cached_object_list:List of 1
##   .. .. .. .. .. ..$ :<environment: 0x04dd3e38> 
##   .. .. ..@ ranges         :Formal class 'GroupedIRanges' [package "XVector"] with 7 slots
##   .. .. .. .. ..@ group          : int [1:64206] 1 1 1 1 1 1 1 1 1 1 ...
##   .. .. .. .. ..@ start          : int [1:64206] 1 29 57 85 114 144 173 202 230 260 ...
##   .. .. .. .. ..@ width          : int [1:64206] 28 28 28 29 30 29 29 28 30 30 ...
##   .. .. .. .. ..@ NAMES          : NULL
##   .. .. .. .. ..@ elementType    : chr "integer"
##   .. .. .. .. ..@ elementMetadata: NULL
##   .. .. .. .. ..@ metadata       : list()
##   .. .. ..@ elementType    : chr "BString"
##   .. .. ..@ elementMetadata: NULL
##   .. .. ..@ metadata       : list()
```








































