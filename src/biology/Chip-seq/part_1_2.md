# ChIP-seq 分析：数据比对（3）



- 读取 = reads（二者含义相同，下文不做区分）

## 1. ChIPseq reads 比对

在评估读取质量和我们应用的任何读取过滤之后，我们将希望将我们的读取与基因组对齐，以便识别任何基因组位置显示比对读取高于背景的富集。

由于 ChIPseq 读数将与我们的参考基因组连续比对，我们可以使用我们在之前中看到的基因组比对器。生成的 BAM 文件将包含用于进一步分析的对齐序列读取。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230205124335394.png)



## 2. 参考基因组生成

首先，我们需要以 FASTA 格式检索感兴趣的基因组的序列信息。我们可以使用 BSgenome 库来检索完整的序列信息。对于小鼠 mm10 基因组，我们加载包 BSgenome.Mmusculus.UCSC.mm10。

```R
library(BSgenome.Mmusculus.UCSC.mm10)
BSgenome.Mmusculus.UCSC.mm10
```

![BSgenome.Mmusculus.UCSC.mm10](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230205124508908.png)



我们将仅使用主要染色体进行分析，因此我们可能会排除随机和未放置的重叠群。在这里，我们循环遍历主要染色体，并根据检索到的序列创建一个 DNAStringSet 对象。

```R
mainChromosomes <- paste0("chr", c(1:19, "X", "Y", "M"))
mainChrSeq <- lapply(mainChromosomes, function(x) BSgenome.Mmusculus.UCSC.mm10[[x]])
names(mainChrSeq) <- mainChromosomes
mainChrSeqSet <- DNAStringSet(mainChrSeq)
mainChrSeqSet
```

![mainChrSeqSet](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230205124534702.png)



现在我们有了一个 DNAStringSet 对象，我们可以使用 writeXStringSet 来创建我们的 FASTA 序列文件来比对。

```R
writeXStringSet(mainChrSeqSet, "BSgenome.Mmusculus.UCSC.mm10.mainChrs.fa")
```



## 3. 索引创建

我们将使用 subread 背后的 subjunc 算法进行对齐。因此，我们可以使用 Rsubread 包。在我们尝试比对我们的 FASTQ 文件之前，我们需要首先使用 buildindex() 函数从我们的参考基因组构建一个索引。

buildindex() 函数仅采用我们所需的索引名称和要从中构建索引的 FASTA 文件的参数。

```R
library(Rsubread)
buildindex("mm10_mainchrs", "BSgenome.Mmusculus.UCSC.mm10.mainChrs.fa", memory = 8000,
    indexSplit = TRUE)
```

- 请记住：建立索引会占用大量内存，默认情况下设置为 8GB。这对于您的笔记本电脑或台式机来说可能太大了。



## 4. 比对

### 4.1. Rsubread

我们可以使用 Rsubread 包将 FASTQ 格式的原始序列数据与 mm10 基因组序列的新 FASTA 文件进行比对。具体来说，我们将使用 align 函数，因为它利用了 subread 基因组比对算法。

```R
myMapped <- align("mm10_mainchrs", "filtered_ENCFF001NQP.fastq.gz", output_format = "BAM",
    output_file = "Myc_Mel_1.bam", type = "dna", phredOffset = 64, nthreads = 4)
```



### 4.2. Rbowtie2

Bowtie 家族是最著名的对齐算法之一。我们可以使用 Rbowtie2 包访问 Bowtie2。QuasR 包允许访问原始的 Bowtie 对准器，但它有点慢并且需要内存。

```R
library(Rbowtie2)
```

与 Rsubread 一样，Rbowtie2 包要求我们首先创建一个要对齐的索引。我们可以使用 bowtie2_build() 函数来完成此操作，指定我们的 FASTA 文件和所需的索引名称。

```R
bowtie2_build(references = "BSgenome.Mmusculus.UCSC.mm10.mainChrs.fa", bt2Index = file.path("BSgenome.Mmusculus.UCSC.mm10.mainChrs"))
```

然后我们可以使用 bowtie2() 函数对齐我们的 FASTQ 数据，指定我们新创建的索引、SAM 输出的所需名称和未压缩的 FASTQ。

我们需要先解压缩我们的 FASTQ。这里我们使用 remove is FALSE 设置来保持原始压缩的 FASTQ。

```R
library(R.utils)
gunzip("filtered_ENCFF001NQP.fastq.gz", remove = FALSE)

bowtie2(bt2Index = "BSgenome.Mmusculus.UCSC.mm10.mainChrs", samOutput = "ENCFF001NQP.sam",
    seq1 = "filtered_ENCFF001NQP.fastq")
```

由于 Rbowtie2 还输出 SAM 文件，我们需要将其转换为 BAM 文件。我们可以使用 RSamtools 的 asBam() 函数来做到这一点。

```R
bowtieBam <- asBam("ENCFF001NQP.sam")
```

使用 Rbowtie2 时的一个重要考虑因素是其未压缩文件的输入和输出。在命令行上，我们可以将输入流式传输到 Rbowtie2，但在 R 中这不是一个选项。我们需要确保删除任何创建的临时文件（SAM 和/或未压缩的 FASTQ）以避免填满我们的硬盘。我们可以使用 unlink() 函数删除 R 中的文件。

```R
unlink("ENCFF001NQP.sam")
```



### 4.3. 排序

和以前一样，我们分别使用 Rsamtools 包 sortBam() 和 indexBam() 函数对文件进行排序和索引。生成的排序和索引 BAM 文件现在可以用于外部程序，例如 IGV，也可以用于 R 中的进一步下游分析。

```R
library(Rsamtools)
sortBam("Myc_Mel_1.bam", "SR_Myc_Mel_rep1")
indexBam("SR_Myc_Mel_rep1.bam")
```

