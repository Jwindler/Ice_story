# 为什么以及如何在多重假设检验中调整 P 值

低于某个阈值的 P 值通常用作选择相关特征的方法。下面的建议建议如何正确使用它们。



当我们在多个特征上重复测试模型时，就会发生多重假设检验，因为获得一个或多个错误发现的概率会随着测试次数的增加而增加。例如，在基因组学领域，科学家们经常想测试数千个基因中的任何一个是否对感兴趣的结果具有显着不同的活性。

在[本文](https://towardsdatascience.com/why-and-how-to-adjust-p-values-in-multiple-hypothesis-testing-2ccf174cdbf8 "Source")中，我们将介绍几种通过调整模型 p 值来解释多重假设检验的流行方法：

1. False Positive Rate (FPR)
2. Family-Wise Error Rate (FWER)
3. False Discovery Rate (FDR)

并解释何时使用它们是有意义的。

本文档可以用下图概括：

![](https://s2.loli.net/2023/05/07/IlzFiOUcSYqB5E3.png)



## 测试数据

我们将创建一个模拟示例，以更好地理解对 p 值的各种操作如何导致不同的结论。要运行此代码，我们需要安装有 pandas、numpy、scipy 和 statsmodels 库的 Python。

出于本示例的目的，我们首先创建一个包含 1000 个特征的 Pandas DataFrame。其中 990 个 (99%) 的值将从均值 = 0 的正态分布生成，称为 Null 模型。 （在下面使用的函数 norm.rvs() 中，使用 loc 参数设置均值。）其余 1% 的特征将从均值 = 3 的正态分布生成，称为非空模型。我们将使用这些来表示我们想要发现的有趣特征。

```python
import pandas as pd
import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

np.random.seed(42)

n_null = 9900
n_nonnull = 100

df = pd.DataFrame({
    'hypothesis': np.concatenate((
        ['null'] * n_null,
        ['non-null'] * n_nonnull,
    )),
    'feature': range(n_null + n_nonnull),
    'x': np.concatenate((
        norm.rvs(loc=0, scale=1, size=n_null),
        norm.rvs(loc=3, scale=1, size=n_nonnull),
    ))
})
```

对于 1000 个特征中的每一个，如果我们假设它是从 Null 分布生成的，则 p 值是观察到该值至少一样大的概率。

P 值可以从累积分布（来自 scipy.stats 的 norm.cdf() ）计算，它表示获得等于或小于观察值的值的概率。然后计算 p 值，我们计算 1 - norm.cdf() 以找到大于观察到的概率：

```python
df['p_value'] = 1 - norm.cdf(df['x'], loc = 0, scale = 1)
df
```

![](https://s2.loli.net/2023/05/07/q4UI9iEgkrRoMFQ.png)



## False Positive Rate

第一个概念称为误报率，定义为我们标记为“显着”的原假设的一小部分（也称为 I 类错误）。我们之前计算的 p 值可以根据其定义解释为误报率：当我们对 Null 分布进行采样时，它们是获得至少与指定值一样大的值的概率。

出于说明目的，我们将应用一个常见的（神奇的）p 值阈值 0.05，但可以使用任何阈值：

```python
df['is_raw_p_value_significant'] = df['p_value'] <= 0.05
df.groupby(['hypothesis', 'is_raw_p_value_significant']).size()
```

![](https://s2.loli.net/2023/05/07/d7ReswfzIVHF4Tm.png)



请注意，在我们的 9900 个原假设中，有 493 个被标记为“重要”。因此，误报率为：FPR = 493 / (493 + 9940) = 0.053。

FPR 的主要问题是，在真实场景中，我们无法先验地知道哪些假设为零，哪些不为零。然后，原始 p 值本身（误报率）的用途有限。在我们的例子中，当非空特征的比例非常小时，大多数标记为重要的特征将是空的，因为它们的数量更多。具体来说，在 92 + 493 = 585 个标记为真（“正”）的特征中，只有 92 个来自我们的非空分布。这意味着大多数或大约 84% 的报告重要特征 (493 / 585) 是误报！

那么，我们能做些什么呢？有两种常见的方法可以解决这个问题：我们可以计算 Family-Wise Error Rate (FWER) 或 False Discovery Rate (FDR)，而不是误报率。这些方法中的每一种都将一组原始的、未经调整的 p 值作为输入，并生成一组新的“调整后的 p 值”作为输出。这些“调整后的 p 值”表示 FWER 和 FDR 上限的估计值。它们可以从 multipletests() 函数中获得，它是 statsmodels Python 库的一部分：

```python
def adjust_pvalues(p_values, method):
   return multipletests(p_values, method = method)[1]
```



## Family-Wise Error Rate

Family-Wise Error Rate 是错误地拒绝一个或多个原假设的概率，或者换句话说，将真 Null 标记为 Non-null，或者看到一个或多个误报的概率。

当只有一个假设被检验时，这等于原始 p 值（假阳性率）。然而，测试的假设越多，我们就越有可能得到一个或多个误报。有两种流行的估计 FWER 的方法：Bonferroni 和 Holm 程序。尽管 Bonferroni 和 Holm 程序都没有对运行测试对单个特征的依赖性做出任何假设，但它们会过于保守。例如，在所有特征都相同的极端情况下（相同模型重复 10,000 次），不需要校正。而在另一个极端，没有特征相关，则需要某种类型的校正。

### Bonferroni 程序

Bonferroni 程序是校正多重假设检验的最流行方法之一。这种方法之所以流行，是因为它非常容易计算，甚至可以手工计算。此过程将每个 p 值乘以执行的测试总数，或者如果此乘法会使它超过 1，则将其设置为 1。

```python
df['p_value_bonf'] = adjust_pvalues(df['p_value'], 'bonferroni')
df.sort_values('p_value_bonf')
```

![](https://s2.loli.net/2023/05/07/GclQqszXDur2CI8.png)



### Holm 程序

Holm 的程序提供了比 Bonferroni 的程序更强大的修正。唯一的区别是 p 值并未全部乘以测试总数（此处为 10000）。相反，每个排序的 p 值逐渐乘以递减序列 10000、9999、9998、9997、...、3、2、1。

```python
df['p_value_holm'] = adjust_pvalues(df['p_value'], 'holm')
df.sort_values('p_value_holm').head(10)
```

![](/home/jzj/.config/Typora/typora-user-images/image-20230507140852173.png)



我们可以自己验证这一点：此输出的最后 10 个 p 值乘以 9991：7.943832e-06 * 9991 = 0.079367。 Holm 修正也是 R 语言中 p.adjust() 函数默认调整 p 值的方法。

如果我们再次应用 0.05 的 p 值阈值，让我们看看这些调整后的 p 值如何影响我们的预测：

```python
df['is_p_value_holm_significant'] = df['p_value_holm'] <= 0.05
df.groupby(['hypothesis', 'is_p_value_holm_significant']).size()
```

![](https://s2.loli.net/2023/05/07/yIgUTH1vMSNmbBA.png)



这些结果与我们对原始 p 值应用相同阈值时的结果大不相同！现在，只有 8 个特征被标记为“重要”，并且所有 8 个都是正确的——它们是从我们的非空分布生成的。这是因为即使一个特征被错误标记的概率也只有 0.05 (5%)。

然而，这种方法有一个缺点：它无法将其他 92 个非空特征标记为重要。虽然确保没有空特征滑入是非常严格的，但它只能找到 8%（100 个中的 8 个）非空特征。这可以看作是采用与误报率方法不同的极端方法。

还有更多的中间立场吗？答案是“是”，中间立场是错误发现率。



## False Discovery Rate

如果我们可以接受一些误报，但捕获超过个位数百分比的真阳性怎么办？也许我们可以接受一些误报，只是没有那么多以至于它们淹没了我们标记为重要的所有特征——就像 FPR 示例中的情况一样。

这可以通过将错误发现率（而不是 FWER 或 FPR）控制在指定的阈值水平（例如 0.05）来实现。错误发现率定义为所有标记为阳性的特征中误报的分数：FDR = FP / (FP + TP)，其中 FP 是误报的数量，TP 是真阳性的数量。通过将 FDR 阈值设置为 0.05，我们表示我们可以接受在我们标记为正的所有特征中有 5%（平均）的误报。

有几种方法可以控制 FDR，这里我们将介绍如何使用两种流行的方法：Benjamini-Hochberg 和 Benjamini-Yekutieli 程序。这两个程序虽然比 FWER 程序更复杂，但很相似。他们仍然依赖于对 p 值进行排序，将它们与特定数字相乘，然后使用截止标准。



### Benjamini-Hochberg 程序

Benjamini-Hochberg (BH) 程序假定每个测试都是独立的。例如，如果被测试的特征相互关联，就会发生相关测试。让我们计算 BH 调整后的 p 值，并将其与我们之前使用 Holm 校正的 FWER 结果进行比较：

```python
df['p_value_bh'] = adjust_pvalues(df['p_value'], 'fdr_bh')
df[['hypothesis', 'feature', 'x', 'p_value', 'p_value_holm', 'p_value_bh']] \
    .sort_values('p_value_bh') \
    .head(10)
```

![](https://s2.loli.net/2023/05/07/bg9sdPpuwYNZnQh.png)



```python
df['is_p_value_holm_significant'] = df['p_value_holm'] <= 0.05
df.groupby(['hypothesis', 'is_p_value_holm_significant']).size()
```

![](https://s2.loli.net/2023/05/07/64rwXgDAFJHmpQ8.png)



```python
df['is_p_value_bh_significant'] = df['p_value_bh'] <= 0.05
df.groupby(['hypothesis', 'is_p_value_bh_significant']).size()
```

![](https://s2.loli.net/2023/05/07/MFBVbQ6R2P8tXpU.png)



BH 程序现在正确地将 100 个非空特征中的 33 个标记为重要特征——比 Holm 修正后的 8 个特征有所改进。但是，它也将 2 个无效特征标记为重要。因此，在标记为重要的 35 个特征中，不正确特征的比例为：2 / 33 = 0.06，即 6%。

请注意，在这种情况下，我们的 FDR 率为 6%，尽管我们的目标是将其控制在 5%。 FDR将平均控制在5%的利率：有时可能较低，有时可能较高。



### Benjamini-Yekutieli 程序

无论测试是否独立，Benjamini-Yekutieli (BY) 程序都会控制 FDR。同样，值得注意的是，所有这些程序都试图在 FDR（或 FWER）上建立上限，因此它们可能更保守或更保守。让我们将 BY 过程与上面的 BH 和 Holm 过程进行比较：

```python
df['p_value_by'] = adjust_pvalues(df['p_value'], 'fdr_by')
df[['hypothesis', 'feature', 'x', 'p_value', 'p_value_holm', 'p_value_bh', 'p_value_by']] \
    .sort_values('p_value_by') \
    .head(10)
```

![](https://s2.loli.net/2023/05/07/z2JQSiMyjoeqkZF.png)



```python
df['is_p_value_by_significant'] = df['p_value_by'] <= 0.05
df.groupby(['hypothesis', 'is_p_value_by_significant']).size()
```

![](https://s2.loli.net/2023/05/07/r7kQZKM5TPiaXSm.png)



BY程序对FDR的控制更为严格；在这种情况下，甚至比 Holm 控制 FWER 的过程更重要，它只将 7 个非空特征标记为重要！使用它的主要优点是当我们知道数据可能包含大量相关特征时。然而，在那种情况下，我们可能还想考虑过滤掉相关的特征，这样我们就不需要测试所有的特征。