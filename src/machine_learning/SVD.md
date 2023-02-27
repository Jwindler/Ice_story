# 降维算法: 奇异值分解SVD



## 1. 为什么降维

总所周知，在低维下，数据更容易处理，但是在通常情况下我们的数据并不是如此，往往会有很多的特征，进而就会出现很多问题：

1. 多余的特征会影响或误导学习器
2. 更多特征意味着更多参数需要调整，过拟合风险也越大
3. 数据的维度可能只是虚高，真实维度可能比较小
4. 维度越少意味着训练越快，更多东西可以尝试，能够得到更好的结果
5. 如果我们想要可视化数据，就必须限制在两个或三个维度上

因此，我们需要通过降维（dimensionality reduction）把无关或冗余的特征删掉。

- 现有降维方法：

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230227165736772.png)



## 2. SVD 概述

奇异值分解（Singular Value Decomposition）简称SVD，主要作用是简化数据，提取信息。

利用SVD实现，我们能够用小得多的数据集来表示原始数据集。这样做，实际上是去除了噪声和冗余信
息。当我们试图节省空间时，去除噪声和冗余信息就是很崇高的目标了，但是在这里我们则是从数据中
抽取信息。基于这个视角，我们就可以把SVD看成是**从有噪声数据中抽取相关特征**。

- SVD是如何从这些充满着大量噪声的数据中抽取相关特征呢？

SVD的公式：

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230227170107436.png)	



这个公式中， U 和 V 都是正交矩阵，即：

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230227170120149.png)



原始数据集A是一个m行n列的矩阵，它被分解成了三个矩阵，分别是：

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230227170141306.png)



这个公式用到的就是矩阵分解技术。在线性代数中还有很多矩阵分解技术。矩阵分解可以将原始矩阵
表示成新的易于处理的形式，这种新形式是两个或多个矩阵的乘积。

不同的矩阵分解技术具有不同的性质，其中有些更适合于某个应用，有些则更适合于其他应用。最常
见的一种矩阵分解技术就是SVD。

- Example

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230227170223318.png)



- Example

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230227170246222.png)



## 3. SVD 的应用

### 3.1. 信息检索

最早的SVD应用之一就是信息检索。利用SVD方法为隐形语义索引（Latent Semantic Indexing，LSI）或者隐形语义分析（Latent Semantic Analysis，LSA）。

在LSI中，一个矩阵是由文档和词语组成的。当我们在该矩阵上应用SVD时，就会构建出多个奇异值。这些奇异值代表了文档中的概念或主题，这一特点可以用于更高效的文档搜索。在词语拼写错误时，只基于词语存在与否的简单搜索方法会遇到问题。简单搜索的另一个问题就是同义词的使用。这就是说，当我们查找一个词时，其同义词所在的文档可能并不会匹配上。如果我们从上千篇相似的文档中抽取出概念，那么同义词就会映射为同一概念。 这样就可以大大提高文档搜索的效率。



### 3.2. 推荐系统

SVD的另外一个应用就是推荐系统。也是目前SVD最主要的一个应用简单版本的推荐系统能够计算项或者人之间的相似度。更先进的方法则先利用SVD从数据中构建一个主题空间，然后再在该空间下计算其相似度。