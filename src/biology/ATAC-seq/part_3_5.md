# ATAC-seq分析：从头识别 Motifs（16）



## 1. 从头识别

到目前为止，我们已经回顾了如何识别 ATACseq 峰中的已知基序。在我们的 ChIPseq 培训中，我们介绍了如何使用 Meme-ChIP 软件执行从头基序识别。在这里，我们将以判别模式使用 `Meme-ChIP` 软件来识别在一组样本中比另一组样本丰富的 De novo 基序。



## 2. 差异峰

首先，我们需要确定一组相对于另一组富含 `ATACseq` 信号的峰。为此，我们将使用 `DESeq2` 来识别 `ATACseq `峰在一组中比另一组更丰富。

我们可以将之前的肝脏、肾脏和大脑计数加载到 `DEseq` 对象中以进行差异分析。

```R
require(DESeq2)
load("data/myCounts.RData")
Group <- factor(c("HindBrain", "HindBrain", "Kidney", "Kidney", "Liver", "Liver"))
colData(myCounts) <- DataFrame(data.frame(Group, row.names = colnames(myCounts)))
dds <- DESeqDataSet(myCounts, design = ~Group)
dds <- DESeq(dds)
```

![dds](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105211837186.png)



然后我们可以使用 `DESeq2` 函数结果来识别不同的 `ATAC` 区域并将格式参数设置为 `GRanges` 以允许将结果作为 `GRanges` 对象返回。

```R
myRes <- results(dds, contrast = c("Group", "Liver", "Kidney"), format = "GRanges")
myRes <- myRes[order(myRes$padj), ]
upRegions <- myRes[myRes$log2FoldChange > 0][1:1000]
downRegions <- myRes[myRes$log2FoldChange < 0, ][1:1000]
upRegions
```

![upRegions](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105211911574.png)



## 3. 调整峰的大小

对于 `Meme-ChIP`，我们希望将区域大小调整为 100bp。默认情况下，`Meme-ChIP` 会为我们执行此修剪，但为了确保它位于中心，我们会提前执行此操作。

```R
upRegions <- resize(upRegions, fix = "center", width = 100)
downRegions <- resize(downRegions, fix = "center", width = 100)
```



## 4. 峰序列

我们现在可以像以前一样使用 `getSeq` 函数从 `GRange` 区域内提取信号，并使用 `writeXStringSet` 函数写入 `FASTA` 文件以用于 `Meme-ChIP`。

```R
library(BSgenome.Mmusculus.UCSC.mm10)
upStrings <- getSeq(BSgenome.Mmusculus.UCSC.mm10, upRegions)
downStrings <- getSeq(BSgenome.Mmusculus.UCSC.mm10, downRegions)
names(upStrings) <- as.character(upRegions)
names(downStrings) <- as.character(downRegions)
writeXStringSet(upStrings, file = "UpRegions.fa")
writeXStringSet(downStrings, file = "DownStrings.fa")
```



## 5. Meme-ChIP

然后我们可以将样本提交到[此处](http://meme-suite.org/tools/meme-chip "Meme-ChIP")的 Meme-ChIP 在线提交表格。与 `ChIPseq` 相比，我们将在差分模式下运行它，以比较我们在肝脏富集区域和肝脏中耗尽区域的序列。结果可以在 data/memeResult_LiverVsBrain 中找到



## 6. 使用 de novo motifs

我们可能想在其他 `Bioconductor` 软件中使用 `Meme-ChIP` 发现的基序。为此，我们可以使用 `universalmotif `包。

```R
library(universalmotif)
```



## 7. 导入 motifs

`universalmotif` 包提供了许多用于导入 `motif` 集的有用函数。一个有用且相关的函数是 `read_meme()` 函数。

```R
memeMotifs <- read_meme("data/memeResult_LiverVsBrain/combined.meme")
memeMotifs
```

![memeMotifs](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105212433595.png)



`universalmotif` 包还提供了在来自不同包的`motif`对象之间进行转换的函数。这里我们转换为 `TFBStools `  `motif`对象。

```R
memeMotifsTFBStools <- convert_motifs(memeMotifs, "TFBSTools-PWMatrix")
memeMotifsTFBStools
```

![memeMotifsTFBStools](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105212508290.png)