# ATACseq分析：教程简介（1）



## 简介

[本课程](https://rockefelleruniversity.github.io/RU_ATACseq/ "Source")介绍 `Bioconductor` 中的 `ATACseq` 分析。

该课程由 2 个部分组成。这将引导您完成正常 `ATACseq` 分析工作流程的每个步骤。它涵盖比对、`QC`、`peak calling`、基因组富集测试、基序富集和差异可及性测试。



## 环境准备

### IGV

IGV 可以从 `BROAD` 网站安装。 》 https://www.broadinstitute.org/igv/



### MACS2

[MACS2](https://macs3-project.github.io/MACS/ "MACS") 没有 R 包，但 MACS2 可在适用于 Linux 或 MacOS 的 Anaconda 包存储库中找到。安装 MACS2 的最简单方法是使用 R 包 Herper。 Herper 允许您从 R 中管理和安装 Anaconda 包。

```R
BiocManager::install("Herper")
library(Herper)
```

安装 Herper 后，您可以使用 install_CondaTools 函数安装 MACS2。在幕后，Herper 将安装最小版本的 conda（称为 miniconda），然后创建一个新环境来安装 MACS2。当您运行该函数时，它会打印出 MACS2 的安装位置。

env 参数是您要为创建的环境指定的名称。 pathToMiniConda 指定您要安装 Miniconda 的位置，以及所有 conda 工具（如 MACS2）。

```R
install_CondaTools(tools="macs2", env="PeakCalling_analysis", pathToMiniConda="/path/to/install")
```



### R

见



### RStudio

见



### 包

- 课程包

```R
install.packages('BiocManager')
BiocManager::install('RockefellerUniversity/RU_ATACseq',subdir='atacseq')
```



- 来自 CRAN 和 Bioconductor

```R
install.packages('BiocManager')
BiocManager::install('methods')
BiocManager::install('ggplot2')
BiocManager::install('rmarkdown')
BiocManager::install('ShortRead')
BiocManager::install('ashr')
BiocManager::install('ChIPQC')
BiocManager::install('DiffBind')
BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')
BiocManager::install('Rsubread')
BiocManager::install('Rbowtie2')
BiocManager::install('R.utils')
BiocManager::install('Rsamtools')
BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
BiocManager::install('rtracklayer')
BiocManager::install('ChIPseeker')
BiocManager::install('soGGi')
BiocManager::install('GenomicAlignments')
BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')
BiocManager::install('DESeq2')
BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')
BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')
BiocManager::install('tracktables')
BiocManager::install('clusterProfiler')
BiocManager::install('TxDb.Mmusculus.UCSC.mm10.knownGene')
BiocManager::install('devtools')
BiocManager::install('tidyr')
BiocManager::install('DT')
BiocManager::install('dplyr')
BiocManager::install('rGREAT')
BiocManager::install('MotifDb')
BiocManager::install('Biostrings')
BiocManager::install('GenomicRanges')
BiocManager::install('pheatmap')
BiocManager::install('universalmotif')
BiocManager::install('seqLogo')
BiocManager::install('org.Mm.eg.db')
BiocManager::install('ATACseqQC')
BiocManager::install('JASPAR2020')
BiocManager::install('motifmatchr')
BiocManager::install('chromVAR')
BiocManager::install('ggseqlogo')
BiocManager::install('TFBSTools')
BiocManager::install('motifStack')
BiocManager::install('knitr')
BiocManager::install('testthat')
BiocManager::install('yaml')
```



## 内容

### Part_1

本节介绍Bioconductor Session部分对ATACseq数据的分析：

- 在 R 中预处理 ATACseq 数据
- 数据比对
- 为可视化创建 bigWig

![Session 1](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221226213009920.png)



### Part_2

本节演示如何使用 `ATACseq` 数据评估可访问性的全局变化。会话部分：

- 在 R 中注释 ATACseq 数据
- 绘制无核小体和单核小体信号
- 绘制 DNA 结合蛋白周围的切割位点

![Session 2](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221226213131724.png)



### Part_3

本节演示如何评估 ATAC-seq 数据中的基序。

- 从数据库中检索 `motifs`。
- 绘制 `motifs`。
- 识别已知的 `motifs`。

![Session 3](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221226213247919.png)



