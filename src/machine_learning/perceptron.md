# 使用 Python 探索 **感知**机 算法

从理论到实践，我们将从简要的理论介绍开始研究感知机(器)学习方法，然后实现。

在这篇[博文](https://towardsdatascience.com/exploring-the-perceptron-algorithm-using-python-c1d3af53a7c7 "Source")的最后，您将能够了解何时以及如何使用这种机器学习算法，清楚地了解它的所有优缺点。



## 1. 理论

### 1.1. 引言

感知器有其存在的生物学原因。我们的神经元不断从其他神经元接收能量，但只有在它们接收到的能量大于或等于一定量后，它们才会决定“激活”并发出自己的信号。

让我们从最后开始。给定一个 4 维输入，这个输入用 4 个不同的权重进行处理，总和进入激活函数，你就得到了结果。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230301144657198.png)



讲的更清除一点，假设您有这张功能表（列）X1、X2、X3 和 X4。这些特征是 4 个不同的值，用于表征数据集的单个实例（行）。

这个实例需要进行二进制分类，这样你就会有一个额外的值 t，它是目标，可以是 -1 或 1。

感知机算法将 X1、X2、X3 和 X4 乘以一组 4 个权重。出于这个原因，我们认为感知器是一种线性算法。

然后，激活函数将应用于此乘法的结果。

这是整个过程中的方程式：

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230301144933437.png)



其中 a 是所谓的激活函数。

当然，输入可以是 N 维的（N 不一定是四维），这样您也可以使用 N 权重 + 1 偏差。尽管如此，纯感知器算法旨在用于二进制分类。

当然，y=a(w_1x_1+…+w_4x_4)的结果需要在-1到1之间。换句话说，归根结底，所谓的激活函数需要能够给你一个分类。

N 维输入与 N 维权重的乘积将为您提供一个数字。那么如果这个数字大于 0，你的算法会说“1”，否则会说“-1”。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230301145313524.png)



这就是它的运行方式，也是它做出决定的方式。



### 1.2. 损失函数

我们都知道机器学习算法带有损失函数。在这种情况下，损失函数是错误分类点的加权和。

假设您有一个分类不正确的点。这意味着，例如，将您的参数与您的输入相乘，您将得到 -0.87 的最终结果。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230301145541015.png)



好的，重点来了，错误分类，记得吗？因此，这意味着该点 (t=1) 的目标确实为“1”。所以这意味着如果你做这个乘法：

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230301145620062.png)



你实际上得到了一个数，告诉你你错了多少，你应该改变你的权重和 bias 来做更好的分类工作。

一般来说，损失函数是所有错误分类点的负和：

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230301145704427.png)



其中S是错误分类点的集合。

我们将开始优化这个损失函数，当然我们想要最小化。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230301145752952.png)



您在上面看到的等式称为梯度下降。这意味着我们遵循损失达到最小值的方向，并按照这个方向更新参数。

由于损失函数取决于错误分类点的数量，这意味着我们将慢慢开始纠正实例，直到如果数据集是线性可分的，将不再有目标“正确”，我们的分类任务将是完美的。



## 2. 实现

当然，SkLearn Perceptron 是众所周知的现成实现。尽管如此，为了更好地理解它，让我们从头开始创建这个感知器。

让我们从库开始：

```python
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('ggplot')
plt.rcParams['font.family'] = 'sans-serif' 
plt.rcParams['font.serif'] = 'Ubuntu' 
plt.rcParams['font.monospace'] = 'Ubuntu Mono' 
plt.rcParams['font.size'] = 14 
plt.rcParams['axes.labelsize'] = 12 
plt.rcParams['axes.labelweight'] = 'bold' 
plt.rcParams['axes.titlesize'] = 12 
plt.rcParams['xtick.labelsize'] = 12 
plt.rcParams['ytick.labelsize'] = 12 
plt.rcParams['legend.fontsize'] = 12 
plt.rcParams['figure.titlesize'] = 12 
plt.rcParams['image.cmap'] = 'jet' 
plt.rcParams['image.interpolation'] = 'none' 
plt.rcParams['figure.figsize'] = (10, 10
                                 ) 
plt.rcParams['axes.grid']=True
plt.rcParams['lines.linewidth'] = 2 
plt.rcParams['lines.markersize'] = 8
colors = ['xkcd:pale range', 'xkcd:sea blue', 'xkcd:pale red', 'xkcd:sage green', 'xkcd:terra cotta', 'xkcd:dull purple', 'xkcd:teal', 'xkcd: goldenrod', 'xkcd:cadet blue',
'xkcd:scarlet']
bbox_props = dict(boxstyle="round,pad=0.3", fc=colors[0], alpha=.5)
```

让我们定义决策函数：

```python
def step_func(z):
        return 1.0 if (z > 0) else 0.0
```



### 2.1. 线性数据

让我们使用 SkLearn 创建一个线性可分的数据集。

```python
from sklearn import datasets
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler


X, y = datasets.make_blobs(n_samples=150,n_features=2,
                           centers=2,cluster_std=3.20)
y[y==0]=-1

#Plotting

min_max_scaler = MinMaxScaler()
X = min_max_scaler.fit_transform(X)

fig = plt.figure(figsize=(10,8))
plt.plot(X[:, 0][y == -1], X[:, 1][y == -1], 'r^')
plt.plot(X[:, 0][y == 1], X[:, 1][y == 1], 'bs')
plt.xlabel("feature 1")
plt.ylabel("feature 2")
plt.title('Random Classification Data with 2 classes')
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230301145944057.png)



### 2.2. 感知器函数

使用这个函数，实际上实现了之前讲解过的所有思路：

```python
def perceptron(X, y, lr, epochs):
    
    # X --> Inputs.
    # y --> labels/target.
    # lr --> learning rate.
    # epochs --> Number of iterations.
    
    # m-> number of training examples
    # n-> number of features 
    m, n = X.shape
    
    # Initializing parapeters(theta) to zeros.
    # +1 in n+1 for the bias term.
    theta = np.zeros((n+1,1))
    
    # Empty list to store how many examples were 
    # misclassified at every iteration.
    n_miss_list = []
    loss_list = []
    # Training.
    for epoch in range(epochs):
        
        # variable to store #misclassified.
        n_miss = 0
        
        # looping for every example.
        for idx, x_i in enumerate(X):
            
            # Insering 1 for bias, X0 = 1.
            x_i = np.insert(x_i, 0, 1).reshape(-1,1)
            
            # Calculating prediction/hypothesis.
            y_hat = step_func(np.dot(x_i.T, theta))
            if y_hat==0:
              y_hat = -1
            # Updating if the example is misclassified.
            if (np.squeeze(y_hat) - y[idx]) != 0:
                theta += lr*((y[idx] - y_hat)*x_i)

                # Incrementing by 1.
                n_miss += 1
        #Defining the loss function
        x1 = X[:,0]
        x2 = X[:,1]
        theta_array = theta
        loss_value = (theta_array[1]*x1+theta_array[2]*x2+theta_array[0])*y
        loss_value = loss_value.sum()/len(x1)
        loss_list.append(loss_value)
        # Appending number of misclassified examples
        # at every iteration.
        n_miss_list.append(n_miss)

    return theta, n_miss_list,loss_list
```

然后我们可以使用以下代码绘制决策边界：

```python
def plot_decision_boundary(X, theta):
    
    # X --> Inputs
    # theta --> parameters
    
    # The Line is y=mx+c
    # So, Equate mx+c = theta0.X0 + theta1.X1 + theta2.X2
    # Solving we find m and c
    x1 = [min(X[:,0]), max(X[:,0])]
    m = -theta[1]/theta[2]
    c = -theta[0]/theta[2]
    x2 = m*x1 + c
    
    # Plotting
    fig = plt.figure(figsize=(10,8))
    plt.plot(X[:, 0][y==-1], X[:, 1][y==-1], "r^")
    plt.plot(X[:, 0][y==1], X[:, 1][y==1], "bs")
    plt.xlabel("Feature 1")
    plt.ylabel("Feature 2")
    plt.title('Perceptron Algorithm')
    plt.plot(x1, x2, 'y-')
```

那么让我们看看玩具数据集中发生了什么：

```python
learning_rate , epoch = 0.005,200
theta, miss_l,loss_list= perceptron(X, y, learning_rate, epoch)
plot_decision_boundary(X, theta)
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230301150106276.png)



可以看出，所有的点都被很好地分类了（即使是小的红色三角形）。

让我们看看损失函数图：

```python
def plot_training(miss_l):
  plt.figure(figsize=(12,12))
  list_array = np.arange(0,len(miss_l),1)
  plt.xlabel('Number of Epochs')
  plt.ylabel('Number of Wrong Classified Points')
  plt.plot(list_array,miss_l)
plot_training(miss_l)
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230301150138601.png)



这意味着数据集现在已经完全分类了。



### 2.3. 非线性数据

让我们考虑一个更难的非线性可分数据集。

```python
from sklearn import datasets
X, y = datasets.make_blobs(n_samples=150,n_features=2,
                           centers=2,cluster_std=3.20)
y[y==0]=-1

#Plotting

min_max_scaler = MinMaxScaler()
X = min_max_scaler.fit_transform(X)

fig = plt.figure(figsize=(10,8))
plt.plot(X[:, 0][y == -1], X[:, 1][y == -1], 'r^')
plt.plot(X[:, 0][y == 1], X[:, 1][y == 1], 'bs')
plt.xlabel("feature 1")
plt.ylabel("feature 2")
plt.title('Random Classification Data with 2 classes')
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230301150317560.png)



让我们运行算法：

```python
def plot_decision_boundary(X, theta):
    
    # X --> Inputs
    # theta --> parameters
    
    # The Line is y=mx+c
    # So, Equate mx+c = theta0.X0 + theta1.X1 + theta2.X2
    # Solving we find m and c
    x1 = [min(X[:,0]), max(X[:,0])]
    m = -theta[1]/theta[2]
    c = -theta[0]/theta[2]
    x2 = m*x1 + c
    
    # Plotting
    fig = plt.figure(figsize=(10,8))
    plt.plot(X[:, 0][y==-1], X[:, 1][y==-1], "r^")
    plt.plot(X[:, 0][y==1], X[:, 1][y==1], "bs")
    plt.xlabel("Feature 1")
    plt.ylabel("Feature 2")
    plt.title('Perceptron Algorithm')
    plt.plot(x1, x2, 'y-')
    plt.xlim(-0.1,1.1)
    plt.ylim(-0.1,1.1)
    
theta, miss_l,loss = perceptron(X, y, 0.2, 10)

plot_decision_boundary(X, theta)

theta, miss_l,loss = perceptron(X, y, 1, 20)

plot_decision_boundary(X, theta)    
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230301150348220.png)



![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230301150402613.png)



好的，现在我们可能需要做一些工作才能获得最佳分类。

让我们运行不同数量的 epoch 和不同的学习率（所谓的超参数调整）以获得感知器的最佳版本：

```python
from sklearn.linear_model import Perceptron
num_of_epochs = [10,100,500,1000]
etas = np.linspace(1e-5,1,100)
scores = []
for e in etas:
  for num in num_of_epochs:
    clf = Perceptron(eta0=e,max_iter=num)
    clf.fit(X, y)
    scores.append({'Num':num,'Eta':e.round(5),'Score':clf.score(X, y)})
    
    
import pandas as pd
import seaborn as sns
scores=pd.DataFrame(scores)
pivot = scores.pivot('Num','Eta','Score')


sns.heatmap(data=pivot)
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230301150452138.png)



所以这是最佳的epoch和学习率：

```python
scores[scores.Score==scores.Score.max()]
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230301150524324.png)



## 总结

- 感知器算法很快。其实就是一个线性乘法+阶跃函数的应用。它非常简单易用。
- 当数据集不可线性分离时，算法不会根据损失函数收敛。这意味着该感知器旨在（完美地）仅在线性可分数据集上工作。尽管如此，我们可以对数据集应用转换，并将感知器算法应用于转换后的数据集
- 超参数调整部分可以大大提高算法的性能。



