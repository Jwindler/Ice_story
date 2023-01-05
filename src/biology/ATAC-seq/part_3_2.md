# ATAC-seq分析：Motif 可视化（13）



## 1. 可视化 motifs

`seqLogo` 包提供了一种简单直观的方式，可以使用 `seqLogo` 函数在 `out motifs` 中可视化我们的碱集频率。

一个简单的 `seqLogo`通过碱基相对于同一位置的其他碱基的相对大小显示每个基序位置的碱基的相对频率。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105203651968.png)



要创建 `seqLog`o，我们必须向 `seqLogo` 函数提供 `PPM`。在这里，我们直接从 `MotidDb` 为 `CTCF` 获取 `PPM`。
这里我们将 `ic.scale` 设置为 `FALSE` 以显示 Y 轴上的概率。

```R
library(seqLogo)
CTCFMotifs <- query(MotifDb, "CTCF")
seqLogo::seqLogo(CTCFMotifs[[1]], ic.scale = FALSE)
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105203736304.png)



将概率绘制为信息内容可能很有用。这里的信息内容会在 0 到 2 Bits 之间。每个碱基具有相同概率的位置将得分为 0，而只有 1 个可能碱基的位置将得分为 2。这使您可以快速识别重要的碱基。

```R
seqLogo::seqLogo(CTCFMotifs[[1]])
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105203803189.png)



`TFBSTools` 包和 `JASPAR` 为我们提供了我们不能直接在 `seqLogo` 包中使用的点频矩阵

```R
myMatrix
```

![myMatrix](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105203835520.png)



我们可以通过简单地将列除以它们的总和来将我们的点频率矩阵转换为点概率矩阵。

```R
ppm <- myMatrix/colSums(myMatrix)
ppm
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105203907447.png)



然后我们可以使用 seqLogo 绘制结果矩阵。

```R
seqLogo::seqLogo(ppm)
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105203921976.png)



幸运的是，`TFBSTools` 有自己的 `seqLogo` 函数版本，我们可以将其与它自己的 `ICMatrix` 类之一一起使用。我们只需要使用 `toICM()` 函数转换我们的对象。

```R
GATA2_IC <- toICM(GATA2mat)
TFBSTools::seqLogo(GATA2_IC)
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105203955020.png)



## 2. ggplot 可视化

我们还可以将 `ggplot2` 样式的图形语法与 `ggseqlogo` 包一起使用。

```R
library(ggseqlogo)
library(ggplot2)
ggseqlogo(myMatrix) + theme_minimal()
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105204029725.png)

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230105204036312.png)