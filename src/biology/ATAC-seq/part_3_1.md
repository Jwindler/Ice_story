# ATAC-seq分析：Motif 数据库（12）

与我们知道目标的 ChIPseq 不同，使用 ATACseq 我们知道一个区域是开放的/可访问的，但不知道可能存在哪些转录因子。

我们在上一节结束时简要了解了如何使用 `matchPWM` 扫描峰值下的特定 `motif`。

在接下来的会议中，我们将回顾如何使用 `motifmathr` 研究 ATACseq 中的 `motifs`。



## 1. Motif 来源

`Bioconductor` 提供两种主要的 `Motif` 来源作为数据库包。这些包括：

- [MotifDb](https://www.bioconductor.org/packages/release/bioc/html/MotifDb.html "MotifDb") 包
- [JASPAR](https://www.bioconductor.org/packages/release/data/annotation/html/JASPAR2020.html "JASPAR") 数据库



## 2. MotifDb

MotifDB 包从广泛的来源收集 `motif` 信息，并将它们存储在 DB 对象中，以便与其他 `Bioconductor` 包一起使用。

```R
library(MotifDb)
MotifDb
```

![MotifDb](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105201828035.png)



`MotifDb` 对象是称为 `MotifList` 的特殊对象类。

```R
class(MotifDb)
```

![MotifDb](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105201858055.png)



像标准的 `List` 对象一样，我们可以使用长度和名称来获取有关我们对象的一些信息。

```R
length(MotifDb)
```

![MotifDb](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105201925329.png)



```R
MotifNames <- names(MotifDb)
MotifNames[1:10]
```

![MotifNames](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105201939593.png)



## 3. 数据获取

我们还可以使用标准列表访问器直接从我们的列表中访问信息。这里 [ 将子集化为单个 `MotifList`。现在我们可以更清楚地看到 `MotifList` 中保存的信息。

```R
MotifDb[1]
```

![MotifDb](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105202132394.png)



与标准列表一样，将子集化为元素内的对象。这里我们提取位置概率矩阵。

```R
MotifDb[[1]]
```

![MotifDb](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105202207777.png)



```R
colSums(MotifDb[[1]])
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105202227279.png)



我们可以使用 `values()` 函数提取所有 `motif` 元数据信息的 `DataFrame`。

```R
values(MotifDb)[1:2, ]
```

![MotifDb](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105202303621.png)



我们可以使用查询函数通过元数据中的信息对我们的 `MotifList` 进行子集化。

```R
CTCFMotifs <- query(MotifDb, "CTCF")
CTCFMotifs
```

![CTCFMotifs](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105202331349.png)



对于更具体的查询，可以使用多个词进行过滤。

```R
CTCFMotifs <- query(MotifDb, c("CTCF", "hsapiens", "jaspar2018"))
CTCFMotifs
```

![CTCFMotifs](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105202348545.png)



## 4. JASPAR2020 包

`JASPAR` 包的更新频率更高，因此可能包含 `MotifDb` 包中未描述的 `motifs`。我们可以像往常一样加载 `JASPAR` 包并像以前一样查看对象。

```R
library(JASPAR2020)
JASPAR2020
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105202438605.png)



## 5. TFBSTools

为了与 `JASPAR` 包交互，我们将使用来自同一实验室的 `TFBSTools` 包。`JASPAR` 包包含有关 `Motif` 和位置概率矩阵 (PPM) 的信息，而 `TFBSTools` 具有操作这些工具并与之交互的功能。

`TFBStools` 提供的三个与 `JASPAR` 数据库交互的有用函数是 `getMatrixSet`、`getMatrixByID` 和 `getMatrixByID`。

```R
library(TFBSTools)
`?`(getMatrixByID)
```



## 6. 提取

`getMatrixByID` 和 `getMatrixByName` 分别采用 `JASPAR DB` 对象和 `JASPAR ID` 或转录因子名称。
这里我们使用转录因子 GATA2。结果是一个新的对象类 PFMatrix。

```R
GATA2mat <- getMatrixByName(JASPAR2020, "GATA2")
class(GATA2mat)
```

![GATA2mat](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105202708380.png)



JASPAR ID 是独一无二的，因此我们可以使用它们来选择我们想要的确切 `motifs`。

```R
GATA2mat <- getMatrixByID(JASPAR2020, "MA0036.2")
```

列表访问器在这里不起作用，但我们可以使用 ID 函数检索名称。

```R
ID(GATA2mat)
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105202751818.png)



## 7. PFM

要掌握位置频率矩阵 (PFM)，我们可以使用 `Matrix` 或 `as.matrix` 函数。

```R
myMatrix <- Matrix(GATA2mat)
myMatrixToo <- as.matrix(myMatrix)
myMatrix
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105203218620.png)



## 8. motif集

我们还可以使用 `getMatrixSet()` 函数来检索 `motifs` 集。我们可以为要检索的 `motifs` 指定一个选项列表。

要查看可用的过滤器，请使用 `getMatrixSet()` 函数的帮助，`?getMatrixSet`。

```R
opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "vertebrates"
motifList <- getMatrixSet(JASPAR2020, opts)
```

