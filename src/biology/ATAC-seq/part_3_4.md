# ATAC-seq分析：信号汇总（15）



## 1. 将 ATAC 信号汇总到 Motifs

`chromVar` 包允许总结 `ATACseq` 信号变化到峰内的 `motifs`。通过这个总结，我们可以潜在地确定与其他样本相比，哪些基序在一组 ATAC 样本中可能具有重要作用。`chromVar` 包来自与 `motimatchr` 包相同的实验室和作者，因此可以很好地协同工作。

```R
library(chromVAR)
```



## 2. 配置

为了识别 `chomVar` 的 `motifs`，我们将使用 `motifMatchr` 和一组不同的输入。在这里，我们将直接向 `matchMotifs` 函数提供包含峰值计数的 `RangedSummarizedExperiment` 对象。

首先，我们将删除所有样本中读取次数少于 5 次的所有峰。

```R
myCounts <- myCounts[rowSums(assay(myCounts)) > 5, ]
```



## 3. GC 校正

接下来我们可以纠正任何可能在测序中出现的潜在 GC 偏差。之前我们已经看到在我们的 FastQ 质量检查中可以观察到 GC 偏差。

为了比较具有不同序列组成的不同组的峰，我们需要纠正这种偏差。我们可以使用 `addGCBias` 并将我们的基因组指定为 `BSgenome` 对象以在 `chromVar` 中进行校正。

```R
myCounts <- addGCBias(myCounts, genome = BSgenome.Mmusculus.UCSC.mm10)
```



## 4. 峰中的基序

纠正偏差后，我们可以再次使用 `matchMotifs` 函数来识别 `ATACseq` 峰下的基序。在这里，我们将峰中计数的 `RangedSummarizedExperiment` 和感兴趣的基因组提供给 `matchMotifs` 函数，并使用默认的匹配项。

```R
motif_ix <- matchMotifs(motifsToScan, myCounts, genome = BSgenome.Mmusculus.UCSC.mm10)
motif_ix
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105211039716.png)



## 5. 运行 chromVar

在我们的峰中识别出 `motifs` 之后，我们可以使用 `computeDeviations` 和 `computeVariability` 函数将 `ATACseq` 信号汇总到  `motifs`。

```R
deviations <- computeDeviations(object = myCounts, annotations = motif_ix)
variability_Known <- computeVariability(deviations)
```



## 6. chromVar 偏差

偏差结果包含一个带有 Z 分数的 `SummarizedExperiment` 对象，显示每个基序的每个样本中 `ATACseq` 信号的富集。

```R
devZscores <- deviationScores(deviations)
devZscores[1:2, ]
```

![devZscores](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105211252619.png)



## 7. chromVar 可变性

可变性结果包含按样本间可变性对基序的排名。高度可变的基序可能表示与特定样本组相关的基序或跨所有组的变量。

```R
variability_Known <- variability_Known[order(variability_Known$p_value), ]
variability_Known[1:10, ]
```

![variability_Known](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105211323986.png)



## 8. chromVar 结果

我们可以使用可变性的结果和我们的 Z 分数偏差来确定我们的 `motifs` 在哪个样本中得到了丰富。

```R
topVariable <- variability_Known[1:20, ]
devTop <- merge(topVariable[, 1, drop = FALSE], devZscores, by = 0)
devTop[1:2, ]
```

![devTop](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105211411547.png)



可视化这些结果的一种有用方法是使用热图。虽然我们将在稍后的会议中介绍这一点，但在这里我们可以使用具有默认设置的 pheatmap 库来说明我们最可变的 `motifs` 在哪里活跃。

```R
devToPlot <- as.matrix(devTop[, -c(1:2)])
rownames(devToPlot) <- devTop[, 2]
library(pheatmap)
pheatmap(devToPlot)
```

![devToPlot](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105211437845.png)