# Pandas 高效数据选择指南



##  简介

存在从 pandas 对象中选择数据子集的不同方法。根据具体的操作，结果要么是指向原始数据的视图，要么是原始数据的副本。这直接关系到操作的效率。复制和查看规则部分源自 NumPy 高级索引规则。我们将研究不同的操作以及如何尽可能地提高性能和效率。

我们还将研究写入时复制将如何改变某些操作的行为以提高性能并尽可能避免复制。



## 数据集

我们将使用包含 FIFA 2021 所有球员的数据集。您可以在[此处](https://www.kaggle.com/datasets/stefanoleone992/fifa-21-complete-player-dataset "FIFA dataset")下载数据集。

```python
import pandas as pd

df = pd.read_csv("players_21.csv", index_col="team_position").sort_index()
```

我们将每个玩家的位置设置为索引，并根据它对 DataFrame 进行排序。这将允许更快、更容易地按位置访问球员，并将帮助我们举例说明。



## 选择行的子集

我们首先从我们的数据集中按位置选择球员。有几种方法可以实现这一点。最常见的可能是通过布尔掩码进行选择。我们可以通过以下方式计算布尔掩码以选择位置为“LS”的所有玩家：

```python
mask = df.index == "LS"
```

之后，我们可以通过以下方式从 DataFrame 中提取行：

```python
result1 = df[mask]
result2 = df.loc[mask]
```

在这种情况下，两种操作都会获得相同的结果。在查看修改我们的 DataFrame 时，我们将调查差异。

通过布尔掩码选择行总是创建数据的副本。根据数据集的大小，这可能会导致速度显着下降。或者，我们可以通过切片对象来选择数据：

```python
result = df.loc["LS"]
```

切片对象会创建底层数据视图，从而显着加快操作速度。您还可以通过以下方式选择每第二行/第 n 行：

```python
result = df.iloc[slice(1, len(df), 2)]
```

这也将创建一个指向原始对象的视图。获取视图通常更可取，因为它可以提高性能并减少内存使用。另一方面，您还可以创建一个与我们的切片相对应的整数列表：

```python
result = df.iloc[list(range(1, len(df), 2))]
```

通过整数列表选择行将创建一个副本，即使操作看起来相似并返回完全相同的数据。这再次源自 NumPy 的索引规则。

切片有很多应用，例如按整数位置、使用 DatetimeIndex 或使用字符串切片索引。如果可能的话，按切片选择数据比使用整数列表或布尔掩码要快得多。

总而言之，根据您的用例，您可以在选择行时显着提高性能。设置适当的索引可能会使您的操作更易于阅读和更有效。



## 选择列的子集

从 DataFrame 中选择列时，通常需要考虑两种情况：

- 选择单个列
- 选择多列

选择单个列相对简单，您可以为此使用常规的 getitem 或 loc。选择数据时，单列没有实质性区别，只有当我们要更新所述数据时。

```python
result = df["long_name"]
result = df.loc[:, "long_name"]
```

一旦将可迭代对象传递给两个调用之一，或者如果所选列重复，我们将返回一个 DataFrame，但会制作基础数据的副本，例如：

```python
result = df.loc[:, ["short_name", "long_name"]]
```

选择多个列通常会立即复制。启用 Copy-on-Write 后，所有这些操作都将返回视图。这将显着提高较大对象的性能。



## 将数据分配给 DataFrame 的子集

让我们看看如何有效地更新 DataFame 的一个子集。一般有两种可能性：常规 setitem 或使用 loc / iloc。

向 DataFrame 添加新列时，我建议使用常规的 setitem 操作。它更短，更容易阅读。两种操作没有本质区别，例如：

```python
df["new_column"] = 100
```

但是，更新 DataFrame 时有很大的不同。假设我们要为对象中位置为“LS”的所有玩家设置名称。常规的 setitem 操作永远不会写入底层数组。此列的数据在更新发生之前被复制。此外，无法在一次操作中更新特定行的子集。您必须使用链式赋值，它有其自身的缺陷。我们稍后会调查它们。

```python
long_name = df[["long_name"]]
long_name[long_name.index == "LS"] = "Testname"
```

在更新所有具有索引“LS”的行之前，我们正在复制整个列。这比使用 loc / iloc 慢得多。如果可能，这两种方法都会更新底层数组。此外，我们不必使用布尔掩码来实现这一点。

```python
df.loc["LS", "long_name"] = "Testname"
```

一般来说，iloc 比 loc 更高效。缺点是，您已经必须知道要插入新值的位置。但是如果你想更新一组特定的行，使用 iloc 比 loc 更有效。

如果要设置的值/值的 dtype 与基础数组的 dtype 兼容，则仅在不制作副本的情况下就地设置值才有效。例如，将整数值设置为 float 或 object dtype 列通常就地操作。将浮点值设置为整数 dtype 列也必须复制数据。整数列不能保存浮点值，因此必须将数据转换为可以保存这两个值的数据类型。作为旁注：如果将不兼容的值设置到列中，则正在进行关于弃用此行为并引发错误的讨论。在设置值之前，需要将列显式转换为浮动。

有一个特殊的例外：覆盖一整列时，使用常规 setitem 通常比使用 loc 更快。

```python
df["long_name"] = "Testname"
```

原因很简单：loc 写入底层数组，这意味着您必须更新此列的每一行。上面的操作只是换出旧列并将新列添加到对象中，而不复制任何内容。



## 链式赋值

链式赋值描述了用一条语句做两个索引操作，然后将数据分配给选定的子集，例如：

```python
df["long_name"][df.index == "LS"] = "Testname"
```

此操作会相应地更新 DataFrame。通常，不应使用链式赋值，因为它是 SettingWithCopyWarning 背后的常见罪魁祸首。此外，链式赋值将在全局启用写时复制或写时复制成为默认设置时引发错误。



## 性能比较

让我们看看这在性能方面意味着什么。这只是一个简单的示例，用于展示如何通过避免副本来提高数据选择的效率。您必须根据您的应用程序对其进行定制。 loc 和 iloc 非常灵活，因此用例会有很大差异。

我们需要更大的数据帧来避免操作中的噪音。我们用随机数实例化一个 DataFrame：

```python
import numpy as np

df = pd.DataFrame(
    np.random.randint(1, 100, (1_000_000, 30)), 
    columns=[f"col_{i}" for i in range(30)],
)
```

让我们看看切片与选择整数列表在性能方面的意义：

```python
%timeit df.loc[slice(10_000, 900_000)]
9.61 µs ± 493 ns per loop (mean ± std. dev. of 7 runs, 100,000 loops each)

%timeit df.loc[list(range(10_000, 900_000))]
68.2 ms ± 465 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
```

对于您的代码的小改动，这是一个非常显着的差异。使用 iloc 显示相同的差异。



## 总结

您可以通过选择最适合您的操作的方法来加快数据选择和数据修改方法。通常，使用切片从 DataFrame 中选择行比使用布尔掩码或整数列表要快得多。设置值时，必须小心使用兼容值。此外，如果修改底层数组没有问题，我们可以通过使用 loc 或 iloc 来提高性能。