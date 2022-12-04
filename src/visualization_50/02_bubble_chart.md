# Matplotlib可视化50图：气泡图（2）



## 导读

[本文](https://towardsdatascience.com/bubble-plots-in-matplotlib-3f0b3927d8f9 "Source")将学习如何使用 `Python` 的 `Matplotlib` 库通过示例绘制气泡图。



## 简介

气泡图是散点图的改进版本。在散点图中，有两个维度 x 和 y。在气泡图中，存在三个维度 x、y 和 z。其中第三维 z 表示权重。这样，气泡图比二维散点图在视觉上提供了更多信息。



## 数据准备

对于本教程，我将使用包含加拿大移民信息的数据集。它拥有从 1980 年到 2013 年的数据，其中包括来自 195 个国家/地区的移民人数。导入必要的包和数据集：

```python
import numpy as np  
import pandas as pd 
df = pd.read_excel('https://s3-api.us-geo.objectstorage.softlayer.net/cf-courses-data/CognitiveClass/DV0101EN/labs/Data_Files/Canada.xlsx',
                       sheet_name='Canada by Citizenship',
                       skiprows=range(20),
                       skipfooter=2)
```

数据集太大。所以，我不能在这里显示完整截图。让我们看看列的名称。

![dataset](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221204194259959.png)



我们不会使用很多列。我只是删除了这些列并将国家名称（“OdName”）设置为索引。

```python
df = df.drop(columns = ['Type', 'Coverage', 'AREA', 'AreaName',      'REG', 'RegName', 'DEV', 'DevName',]).set_index('OdName')
df.head()
```

![example](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221204194416785.png)



我为这个练习选择了爱尔兰和巴西的数据。没有特殊原因。我随机选择了它们。

```python
Ireland = df.loc['Ireland']
Brazil = df.loc['Brazil']
```



## 归一化

有几种不同的方法可以归一化数据。我们将数据归一化以使数据处于相似的范围内。爱尔兰和巴西的移民数据有不同的范围。我需要将它们调整到 0 到 1 的范围内。我只是将爱尔兰数据除以爱尔兰数据系列的最大值。我对巴西的数据系列做了同样的事情。

```python
i_normal = Ireland / Ireland.max()
b_normal = Brazil / Brazil.max()
```

我们将根据年份绘制爱尔兰和巴西的数据。将年份列在清单上会很有用。

```python
years = list(range(1980, 2014))
```



## 可视化

为了看看区别，让我们先绘制散点图。

```python
import matplotlib.pyplot as plt
plt.figure(figsize=(14, 8))
plt.scatter(years, Ireland, color='blue')
plt.scatter(years, Brazil, color='orange')
plt.xlabel("Years", size=14)
plt.ylabel("Number of immigrants", size=14)
plt.show()
```

![scatter](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221204194701072.png)



现在，绘制气泡图。我们必须输入我们之前定义的尺寸。

```python
plt.figure(figsize=(12, 8))
plt.scatter(years, Brazil, 
                  color='darkblue', 
                 alpha=0.5,
                 s = b_normal * 2000)
plt.scatter(years, Ireland, 
                  color='purple', 
                 alpha=0.5,
                 s = i_normal * 2000,
                 )
plt.xlabel("Years", size=14)
plt.ylabel("Number of immigrants", size=14)
```

![bubble](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221204194759670.png)



我们可以通过气泡的大小来了解移民的数量。气泡越小，移民人数越少。

我们也可以让结果更多彩多姿。为了让它有点意义，我们需要对数据系列进行排序。您很快就会看到原因。

```python
c_br = sorted(Brazil)
c_fr = sorted(France)
```

现在我们将传递这些值来改变颜色。

```python
plt.figure(figsize=(12, 8))
plt.scatter(years, Brazil, 
                  c=c_br,
                 alpha=0.5,
                 s = b_normal * 2000)
plt.scatter(years, Ireland, 
                  c=c_fr,
                 alpha=0.5,
                 s = i_normal * 2000,
                 )
plt.xlabel("Years", size=14)
plt.ylabel("Number of immigrants", size=14)
```

![result](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221204194906393.png)



现在我们添加了另一个维度，颜色。颜色随移民数量变化。但是当我们绘制两个变量时，它并没有那么好。因为在这个过程中我们没有明确定义各个变量的颜色。但是当我们在 y 轴上绘制一个变量时，它做得很好。让我们绘制每年来自巴西的移民人数，以了解多年来的趋势。

```python
plt.figure(figsize=(12, 8))
plt.scatter(years, Brazil, 
                  c=c_br,
                 alpha=0.5,
                 s = b_normal * 2000)
plt.xlabel("Years", size=14)
plt.ylabel("Number of immigrants of Brazil", size=14)
```

![color](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221204195056211.png)