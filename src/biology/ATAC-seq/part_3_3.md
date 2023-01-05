# ATAC-seq分析：Motif 鉴定（14）



## 1. 基序识别

由于 `ATACseq` 只告诉我们哪些区域是开放/可访问的，我们可以使用 `ATACseq` 峰区域内的基序来识别可能在那里起作用的潜在转录因子。

`MotifDb` 和 `JASPAR2020` 等数据库为我们提供了一组已知的 `motif` 模式，供我们使用我们的 `motif` 识别。



## 2. motifmatchr

为了识别已知的 `motifs`，我们将使用 `motimatchr` 包，它是 `MOODS c++` 库的包装器。

这意味着 `motifmathr` 为我们提供了一种快速识别 `ATACseq` 数据中基序的方法。

我们将使用具有默认 p-value cut-off 的 `motifmathr/MOODs`。

首先，我们可以检索一组合理的 `motifs` 以在我们的小鼠组织 ATACseq 数据中进行扫描。

在这里，我们检索脊椎动物，`JASPAR CORE motifs`。我们额外指定 `all_versions` 为 `FALSE` 以仅包含 `motifs` 的最新版本。

```R
opts <- list()
opts[["tax_group"]] <- "vertebrates"
opts[["collection"]] <- "CORE"
opts[["all_versions"]] <- FALSE
motifsToScan <- getMatrixSet(JASPAR2020, opts)
```

`motimatchr` 包的主要函数是 `matchMotifs()`。与许多 `Bioconductor` 函数一样，`matchMotifs` 使用其他 `Bioconductor` 对象，例如 `BSGenome`、`GRanges` 和 `summaryExperiment` 对象。



## 3. ENCODE ATAC 数据

我们需要一些数据来审查，以便我们可以使用 `ENCODE` 小鼠组织 `ATAC` 数据。我们可以从数据目录中将计数加载到 `summarizedExperiment` 对象中。

```R
load("data/myCounts.RData")
myCounts
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105204809976.png)



我们可以使用标准访问器 `rowRanges` 直接从 `SummarizedExperiment` 对象中检索峰的范围。

```R
peakRanges <- rowRanges(myCounts)
peakRanges[1, ]
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105204935911.png)



我们可以使用 `getSeq()` 函数从 `BSGenome` 对象中提取序列，使用 `GRanges` 对象定义感兴趣的区域。这是一种类似于 `ChIPseq` 从头发现基序的方法。

我们加载 `BSgenome.Mmusculus.UCSC.mm10` 库并使用调整大小功能将我们的峰重新居中到 `100bp`。调整大小后，我们可以使用 `getSeq()` 函数提取所需的序列。

```R
library(BSgenome.Mmusculus.UCSC.mm10)
peakRangesCentered <- resize(peakRanges, fix = "center", width = 100)
peakSeqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, peakRangesCentered)
names(peakSeqs) <- as.character(peakRangesCentered)
peakSeqs
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105205131737.png)



## 4. 寻找基序位置

正如我们在它的帮助中看到的，`matchMotifs` 函数可以提供 `motif` 匹配的输出作为匹配、分数或位置。

在这里，我们使用前 100 个 `ATACseq` 峰下的序列从 `JASPAR2020` 中扫描 4 个选定的基序，并将输出指定为位置。

```R
motif_positions <- matchMotifs(motifsToScan[1:4], peakSeqs[1:100], out = "positions")
class(motif_positions)
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105205231433.png)



```R
length(motif_positions)
```

![motif_positions](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105205623438.png)



结果包含一个列表，其长度与测试的 `motifs` 数量相同。每个元素都包含一个 `IRangeslist`，其中包含每个测试序列的条目和峰序列内基序位置的 `IRanges`。

```R
motif_positions$MA0029.1
```

![motif_positions](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105205707454.png)



我们可以将我们的 `IRangeslist` 取消列出到标准列表中，以便于工作。

```R
MA0029hits <- motif_positions$MA0029.1
names(MA0029hits) <- names(peakSeqs[1:100])
unlist(MA0029hits, use.names = TRUE)
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105205741572.png)



## 5. 寻找 motif hits

我们可能只是想将基序映射到它们的 `ATACseq` 峰。为此，我们可以将输出参数设置为匹配。这将返回一个 `SummarizedExperiment` 对象。

```R
motifHits <- matchMotifs(motifsToScan, peakSeqs, out = "matches")
class(motifHits)
```

![motifHits](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105210139437.png)



```R
motifHits
```

![motifHits](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105210151170.png)



我们可以使用 `motifMatches` 函数按 `motif` 和 `peak` 检索匹配矩阵。

```R
mmMatrix <- motifMatches(motifHits)
dim(mmMatrix)
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105210359702.png)



```R
mmMatrix[1:8, 1:8]
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105210434792.png)



虽然是稀疏矩阵，但我们仍然可以使用我们的矩阵运算从这个对象中提取有用的信息。我们可以使用 `colSums()` 来识别峰序列中基序的总出现次数。

```R
totalMotifOccurence <- colSums(mmMatrix)
totalMotifOccurence[1:4]
```

![totalMotifOccurence](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105210506157.png)



我们还可以识别包含选定基序命中的峰。

```R
peaksWithMA0912 <- peakRangesCentered[mmMatrix[, "MA0912.2"] == 1]
peaksWithMA0912
```

![peaksWithMA0912](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105210535630.png)