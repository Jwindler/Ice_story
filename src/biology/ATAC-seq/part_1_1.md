# 1_1 比对





### 指控

在比对之前，我们建议花一些时间查看 FASTQ 文件。一些基本的 QC 检查可以帮助我们了解您的测序是否存在任何偏差，例如读取质量的意外下降或非随机 GC 内容。



## Greenleaf

在本节中，我们将稍微处理一下 Greenleaf 数据集。

我们将处理从 FASTQ 到 BAM 的 Greenleaf 数据的一个样本，以允许我们审查 ATACseq 数据的一些特征，并创建一些处理过的文件以供审查和进一步分析。



## 参考基因组

首先，我们需要创建一个参考基因组来比对我们的 ATACseq 数据。我们可以创建一个 FASTA 文件用于从 Bioconductor BSGenome 对象进行比对。

这次我们正在处理人类数据，因此我们将使用 `BSgenome.Hsapiens.UCSC.hg19` 库构建 `hg19` 基因组。

```R
library(BSgenome.Hsapiens.UCSC.hg19)
mainChromosomes <- paste0("chr",c(1:21,"X","Y","M"))
mainChrSeq <- lapply(mainChromosomes,
                     function(x)BSgenome.Hsapiens.UCSC.hg19[[x]])
names(mainChrSeq) <- mainChromosomes
mainChrSeqSet <- DNAStringSet(mainChrSeq)
writeXStringSet(mainChrSeqSet,
                "BSgenome.Hsapiens.UCSC.hg19.mainChrs.fa")
```



## Rsubread

对于 Rsubread，我们必须在 Rsubread 的对齐步骤之前建立我们的索引。

这里我额外指定参数 indexSplit 为 TRUE 并结合 memory 参数设置为 1000 (1000MB) 以控制 Rsubread 对齐步骤中的内存使用。

```R
library(Rsubread)
buildindex("BSgenome.Hsapiens.UCSC.hg19.mainChrs",
           "BSgenome.Hsapiens.UCSC.hg19.mainChrs.fa",
           indexSplit = TRUE,
           memory = 1000)
```



## 比对

