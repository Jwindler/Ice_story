# 如何高效利用 `t-SNE`



## 摘要

尽管`t-SNE`对于可视化高维数据非常有用，但有时其结果可能无法解读或具有误导性。通过探索它在简单情况下的表现，我们可以学会更有效地使用它。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221106215337998.png)



探索高维数据的一种流行方法是`t-SNE`，由 **[van der Maaten 和 Hinton](http://www.jmlr.org/papers/volume9/vandermaaten08a/vandermaaten08a.pdf "PDF")  在 2008 年提出**。该技术已在机器学习领域得到广泛应用，因为它具有几乎神奇的能力，可以从数百甚至数千维的数据中获取其二维的表示。尽管结果令人印象深刻，但这些结果很容易被误读。本文的目的就是指出一些常见的误解。

我们将通过一系列简单的示例来说明 `t-SNE` 图可以显示和不能显示的内容。`t-SNE` 技术确实很有用——但前提是你知道如何解释它。

深入研究之前：如果您以前没有遇到过 `t-SNE`，那么您需要了解它背后的数学知识。其目标是在高维空间中获取一组点，并在低维空间（通常是 2D 平面）中找到这些点的表示。该算法是非线性的，并适应底层数据，对不同区域执行不同的转换。这些差异可能是造成混乱的主要来源。

`t-SNE` 的第二个特征是可调整的参数，`perplexity`，它说明了如何在数据的局部和全局之间平衡注意力。从某种意义上说，该参数是对每个点的近邻数量的猜测。`perplexity`值对生成的图片有复杂的影响。

> 原论文说，“SNE 的性能对 `perplexity`的变化相当稳健，典型值在 5 到 50 之间。”

充分利用 `t-SNE` 可能意味着需要分析具有不同 `perplexity`的多个图。

例如，`t-SNE` 算法并不总是在连续运行中产生类似的输出，并且还有与优化过程相关的超参数。



## 1. 超参数

- 超参数的重要性

让我们从 `t-SNE` 的“hello world”开始：由两个相隔很远的 ``cluster`s` 组成的数据集。为了尽可能简单，我们将考虑二维平面中的`cluster`，如下左图所示。（为了对比，两个`cluster`采用不同的颜色表示。）右下图显示了五种不同 `perplexity` 的 `t-SNE` 图。

![`perplexity` ](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221106220601553.png)



van der Maaten 和 Hinton 建议的 `perplexity` 在 (5 - 50) 范围内，这些图确实显示了这些 ``cluster`s`，尽管形状非常不同。在这个范围之外的结果变得有点奇怪。对于 `perplexity` = 2，局部变化占主导地位。 `perplexity`=100 的图像表明：为了使算法正常运行，`perplexity`应该小于点的数量。否则，可能会产生意想不到的结果。

上面的每个图都是用 5,000 次迭计算作的，学习率（通常称为`epsilon`）为 10，并且在第 5,000 步时结果趋于稳定。

![超参数](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221106221140277.png)



Step 值有多大的影响呢？根据我们的经验，最重要的是继续迭代，直到趋于稳定。

上面的图像显示了在 `perplexity`=30 下的五次不同的运行。前四次在稳定之前停止。在 10、20、60 和 120 步之后，您可以看到具有看似一维甚至点状图像的``cluster`s`布局。如果您看到具有奇怪“挤压”形状的 `t-SNE` 图，则该过程可能过早停止。不幸的是，没有一个固定的`Step`值可以产生稳定的结果。不同的数据集可能需要不同数量的迭代才能收敛。

另一个问题是使用相同超参数的不同运行是否会产生相同的结果。在这个简单的两个簇示例以及我们讨论的大多数其他示例中，多次运行给出了相同的全局形状。然而，某些数据集在不同的运行中会产生明显不同的结果；稍后我们将给出其中之一的示例。

从现在开始，除非另有说明，否则我们将展示 5,000 次迭代的结果。这通常足以使本文中的（相对较小的）示例收敛。然而，我们将继续展示一系列的`perplexities`，因为这似乎在每种情况下都有很大的不同。



## 2. 簇

- `t-SNE` 图中的`cluster`(簇)大小没有任何意义

如果两个 `cluster` 有不同的标准差，大小也不同呢？（尺寸是指边界框测量值，而不是点数。）下面是平面上混合高斯的 `t-SNE` 图，其中一个的分散情况是另一个的 10 倍。

![`cluster`](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221106222413426.png)



令人惊讶的是，这两个 `cluster` 在 `t-SNE` 图中看起来大致相同。`t-SNE` 算法使其“距离”适应数据集中的区域密度变化。结果，它自然地扩展了密集的 `cluster` ，并收缩了稀疏的 `cluster` ，从而平衡了 `cluster`  的大小。需要明确的是，这与任何降维技术都会扭曲距离的情况不同。相反，密度均衡是通过设计产生的，并且是 `t-SNE` 的可预测特征。

然而，您无法在 `t-SNE` 图中看到 `cluster`  的相对大小。



## 3. 距离

- `cluster` 之间的距离可能没有任何意义

下图显示了三个高斯，每个 50 点，一对的距离是另一对的 5 倍。

![Distances ](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221106222854574.png)



`perplexity`为 50 时，该图很好地显示了全局情况。对于较低的 `perplexity`，`cluster` 看起来是等距的。当 `perplexity`为 100 时，我们可以看到全局形状很好，但其中一个 `cluster` 错误地出现，比其他 `cluster` 小得多。因为 `perplexity`=50 在这个例子中得到的结果很好，如果我们想看到全局情况，我们可以总是将 `perplexity`设置为 50 吗？

结果是不能。如果我们向每个`cluster` 添加更多点，则必须增加`perplexity`以进行补偿。这是三个高斯`cluster`的 `t-SNE` 图，每个`cluster` 有 200 个点，而不是 50 个。现在，没有一个`perplexity`给出了好的结果。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221106223452153.png)



想要看到全局情况需要微调`perplexity`，这是个坏消息。真实的数据可能会有多个具有不同数量元素的`cluster` 。可能没有一个 `perplexity`值可以捕获所有`cluster` 的距离——遗憾的是，`perplexity`是一个全局参数。解决这个问题可能是未来研究的一个方向。

因此， `t-SNE` 图中`cluster`之间的距离可能毫无意义。



## 4. 随机噪声

- 随机噪声并不总是看起来随机。

当你看到噪音时，识别它是一项关键技能，但需要时间来建立正确的直觉。`t-SNE` 的一个棘手之处在于它抛弃了很多现有的直觉。下图显示了真正的随机数据，从 100 维的单位高斯分布中抽取 500 个点。左图是前两个坐标的投影。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221106224203812.png)



`perplexity`=2 的情节似乎显示出戏剧性的`cluster`。如果您正在调整`perplexity`以显示数据的结构，您可能会认为自己中了大奖。

当然，因为我们知道点云是随机生成的，所以它没有统计上有趣的`cluster`：那些“`cluster`”没有意义。如果您回顾前面的示例，低`perplexity` 通常会导致这种分布。将这些团块识别为随机噪声是阅读 `t-SNE` 图的重要部分。

起初，`perplexity`=30 图看起来根本不像高斯分布：云的不同区域之间只有轻微的密度差异，而且这些点似乎是均匀分布的。事实上，这些特征说明了关于高维正态分布的有用信息，它们非常接近球体上的均匀分布，点之间的间距大致相等。从这个角度来看，`t-SNE` 图比任何线性投影都更准确。



## 5. 形状

- 有时你可以看到一些形状

很少有数据以完美对称的方式分布。我们来看一个 50 维的轴对齐高斯分布，其中坐标 i 的标准差为 1/`i`。也就是说，我们正在查看一个长椭圆形的点云。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221106224726573.png)



对于足够高的`perplexity`，细长的形状很容易阅读。另一方面，在低`perplexity`度下，局部效应和无意义的“聚集”占据中心位置。更极端的形状也出现了，但同样只是在正确的`perplexity`中。例如，这里有两个 75 个点的 2D 聚类，它们以平行线排列，带有一点噪音。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221106224857654.png)



对于一定范围的`perplexity`，长`cluster`看起来接近正确。

然而，即使在最好的情况下，也存在细微的失真：`t-SNE` 图中的线条略微向外弯曲。原因像往常一样，`t-SNE` 倾向于扩展更密集的数据区域。由于簇的中间比末端周围的空白空间少，因此算法会放大它们。



## 6. 拓扑

- 对于拓扑，您可能需要多次绘图

有时您可以从 `t-SNE` 图上读取拓扑信息，但这通常需要多个`perplexity`的视图。最简单的拓扑属性是包容性。下图显示了 50 维空间中的两组 75 个点。两者都是从以原点为中心的对称高斯分布中采样的，但其中一个的分散度是另一个的 50 倍。 “小”分布实际上包含在大分布中。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221106225114215.png)



`perplexity`=30 视图正确显示了基本的拓扑情况，但 `t-SNE` 再次夸大了较小点组的大小。在`perplexity`=50 处：外部组变成了一个圆圈，因为该图试图描绘它的所有点与内部组的距离大致相同。如果你单独看这张图片，很容易将这些外点误读为一维结构。

考虑一组点，这些点在三个维度上追踪一个链接或一个结。再一次，查看多个`perplexity`值给出了最完整的画面。低`perplexity`给出两个完全独立的循环；高的显示了一种全球连通性。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221106225300740.png)



` trefoil knot `是一个有趣的例子，说明了多次运行如何影响 `t-SNE` 的结果。下面是 `perplexity`为 2 时的五次运行结果。

该算法至少保留了原本的拓扑结构。但是其中的三个得到了不同的结果。使用点颜色作为对比，您可以看到第一次和第三次运行彼此相距很远。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221106225630173.png)



然而，在 `perplexity`=50 的五次运行结果（直到对称）在视觉上是相同的。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221106225657926.png)



## 总结

`t-SNE` 如此受欢迎是有原因的：它非常灵活，并且经常可以找到其他降维算法无法找到的结构。不幸的是，这种灵活性使其难以解释。在用户看不见的地方，该算法会进行各种调整以整理其可视化。好消息是，通过研究 `t-SNE` 在简单情况下的表现，可以对研究有一定帮助。
