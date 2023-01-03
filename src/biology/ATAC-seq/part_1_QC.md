# ATAC-seq分析：数据质控（6）



## 1. 质控

ATACseqQC 库允许我们在一个步骤中运行我们已经看到的许多 ATACseq QC 步骤。它可能会消耗更多内存，但会允许包含两个更有用的指标，称为 PCR 瓶颈系数（PBC1 和 PBC2）。

首先我们必须安装库。

```R
BiocManager::install("ATACseqQC")
```

与 ChIPQC 一样，ATACseqQC 函数包含一个工作流函数，它将通过 BAM 文件路径的单个参数获取大部分所需的 QC。由于这可能会占用大量内存，因此我只是在一个 BAM 文件中对其进行说明，该文件仅包含 ATACseq 数据的 17 号染色体读数。

```R
library(ATACseqQC)
ATACQC <- bamQC("~/Downloads/Sorted_ATAC_50K_2_ch17.bam")
```

生成的 ATACQC 对象具有许多 QC 信息，包括重复率、非冗余分数、跨染色体的信号分布、线粒体分数等。其中包括 PCRbottleneckCoefficient_1 和 PCRbottleneckCoefficient_2 值。

```R
names(ATACQC)
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221229172312605.png)



## 2. PCR 偏差

PCR bottleneck coefficients 确定 ATAC 样品制备过程中可能发生的 PCR 偏差/过度放大。PCRbottleneckCoefficient_1 计算为基因组中恰好有 1 个读取唯一映射的位置数与至少有 1 个读取的位置数相比。例如，如果我们有 20 个读数。 16 独特的位置。 4 没有唯一映射，而是有 2 个位置，两个位置都有 2 个读数。这将引导我们计算 16/18。因此，我们的 PBC1 为 0.889。小于 0.7 的值表示严重瓶颈，0.7 和 0.9 之间表示中度瓶颈。大于 0.9 表示没有瓶颈。

```R
ATACQC$PCRbottleneckCoefficient_1
```

![PCRbottleneckCoefficient_1](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221229172654736.png)



PCRbottleneckCoefficient_2 是我们衡量瓶颈的次要指标。它计算为基因组中恰好 1 个读数唯一映射的位置数与恰好 2 个读数唯一映射的位置数之比。我们可以重用我们的例子。如果我们有 20 个读取，其中 16 个映射唯一。 4 没有唯一映射，而是有 2 个位置，两个位置都有 2 个读数。这将引导我们计算 16/2。因此，我们的 PBC2 为 8。小于 1 的值表示严重瓶颈，1 到 3 之间表示中度瓶颈。大于 3 表示没有瓶颈。

```R
ATACQC$PCRbottleneckCoefficient_2
```

![PCRbottleneckCoefficient_2](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221229172734736.png)