# 机器学习算法简介: 支持向量机



## 1. 简介

我猜你现在已经习惯了线性回归和逻辑回归算法。如果没有，我建议你在继续学习支持向量机之前先看看它们。支持向量机是每个机器学习专家都应该拥有的另一种简单算法。支持向量机是许多人的首选，因为它以较少的计算能力产生显着的准确性。支持向量机，缩写为SVM，可用于回归和分类任务。但是，它广泛用于分类目标。

[本文](https://towardsdatascience.com/support-vector-machine-introduction-to-machine-learning-algorithms-934a444fca47 "Source")将介绍其主要内容和代码实现！



## 2. 什么是 SVM

支持向量机算法的目标是在 N 维空间（N — 特征数）中找到一个超平面，该超平面可以明确地对数据点进行分类。

![Possible hyperplanes](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230215211826129.png)

为了分离两类数据点，可以选择许多可能的超平面。我们的目标是找到一个具有最大边距的平面，即两个类的数据点之间的最大距离。最大化边缘距离提供了一些强化，以便可以更有信心地对未来的数据点进行分类。



## 3. 超平面和支持向量

![Hyperplanes in 2D and 3D feature space](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230215211914688.png)



超平面是帮助对数据点进行分类的决策边界。落在超平面两侧的数据点可以归为不同的类别。此外，超平面的维数取决于特征的数量。如果输入特征的数量是 2，那么超平面就是一条直线。如果输入特征的个数是3，那么超平面就变成了一个二维平面。当特征数超过 3 时，就变得难以想象了。



![Support Vectors](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230215211956528.png)



支持向量是距离超平面较近并影响超平面位置和方向的数据点。使用这些支持向量，我们最大化分类器的边缘。删除支持向量将改变超平面的位置。这些是帮助我们构建 SVM 的要点。



## 4. 大边距

在逻辑回归中，我们获取线性函数的输出并使用 sigmoid 函数将值压缩在 [0,1] 范围内。如果压缩值大于阈值 (0.5)，我们为其分配标签 1，否则我们为其分配标签 0。在 SVM 中，我们采用线性函数的输出，如果该输出大于 1，我们确定它是一个类，如果输出是 -1，我们确定是另一个类。由于在 SVM 中阈值更改为 1 和 -1，我们获得了作为边距的值的强化范围（[-1,1]）。



## 5. 成本函数和梯度更新

在 SVM 算法中，我们希望最大化数据点和超平面之间的边距。有助于最大化边距的损失函数是铰链损失。

![Hinge loss function](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230215212129462.png)



如果预测值和实际值具有相同的符号，则成本为 0。如果不是，则我们计算损失值。我们还在成本函数中添加了一个正则化参数。正则化参数的目标是平衡边缘最大化和损失。添加正则化参数后，成本函数如下所示。

![Loss function for SVM](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230215212149282.png)



现在我们有了损失函数，我们对权重求偏导数来找到梯度。使用梯度，我们可以更新我们的权重。

![Gradients](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230215212205892.png)



当没有错误分类时，即我们的模型正确预测了我们的数据点的类别，我们只需要更新正则化参数的梯度。

![Gradient Update](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230215212221951.png)



当出现错误分类时，即我们的模型在预测我们的数据点的类别时出错，我们将损失与正则化参数一起包括在内以执行梯度更新。

![Gradient Update](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230215212239992.png)



## 6. 代码实现

我们将用来实现 SVM 算法的数据集是 Iris 数据集。您可以从此[链接](https://www.kaggle.com/jchen2186/machine-learning-with-iris-dataset/data "dataset")下载它。

```python
import pandas as pd

df = pd.read_csv('/Users/rohith/Documents/Datasets/Iris_dataset/iris.csv')
df = df.drop(['Id'],axis=1)
target = df['Species']
s = set()
for val in target:
    s.add(val)
s = list(s)
rows = list(range(100,150))
df = df.drop(df.index[rows])
```

由于 Iris 数据集具有三个类，我们将删除其中一个类。这给我们留下了一个二元类分类问题。

```python
import matplotlib.pyplot as plt

x = df['SepalLengthCm']
y = df['PetalLengthCm']

setosa_x = x[:50]
setosa_y = y[:50]

versicolor_x = x[50:]
versicolor_y = y[50:]

plt.figure(figsize=(8,6))
plt.scatter(setosa_x,setosa_y,marker='+',color='green')
plt.scatter(versicolor_x,versicolor_y,marker='_',color='red')
plt.show()
```

![Visualizing data points](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230215212348184.png)



此外，还有四个功能可供我们使用。我们将仅使用两个特征，即萼片长度和花瓣长度。我们采用这两个特征并将它们绘制成可视化。从上图中，您可以推断出可以使用一条直线来分隔数据点。

```python
from sklearn.utils import shuffle
from sklearn.cross_validation import train_test_split
import numpy as np
## Drop rest of the features and extract the target values
df = df.drop(['SepalWidthCm','PetalWidthCm'],axis=1)
Y = []
target = df['Species']
for val in target:
    if(val == 'Iris-setosa'):
        Y.append(-1)
    else:
        Y.append(1)
df = df.drop(['Species'],axis=1)
X = df.values.tolist()
## Shuffle and split the data into training and test set
X, Y = shuffle(X,Y)
x_train = []
y_train = []
x_test = []
y_test = []

x_train, x_test, y_train, y_test = train_test_split(X, Y, train_size=0.9)

x_train = np.array(x_train)
y_train = np.array(y_train)
x_test = np.array(x_test)
y_test = np.array(y_test)

y_train = y_train.reshape(90,1)
y_test = y_test.reshape(10,1)
```

我们提取所需的特征并将其拆分为训练和测试数据。 90% 的数据用于训练，其余 10% 用于测试。现在让我们使用 numpy 库构建我们的 SVM 模型。

```python
## Support Vector Machine 
import numpy as np

train_f1 = x_train[:,0]
train_f2 = x_train[:,1]

train_f1 = train_f1.reshape(90,1)
train_f2 = train_f2.reshape(90,1)

w1 = np.zeros((90,1))
w2 = np.zeros((90,1))

epochs = 1
alpha = 0.0001

while(epochs < 10000):
    y = w1 * train_f1 + w2 * train_f2
    prod = y * y_train
    print(epochs)
    count = 0
    for val in prod:
        if(val >= 1):
            cost = 0
            w1 = w1 - alpha * (2 * 1/epochs * w1)
            w2 = w2 - alpha * (2 * 1/epochs * w2)
            
        else:
            cost = 1 - val 
            w1 = w1 + alpha * (train_f1[count] * y_train[count] - 2 * 1/epochs * w1)
            w2 = w2 + alpha * (train_f2[count] * y_train[count] - 2 * 1/epochs * w2)
        count += 1
    epochs += 1
```

α(0.0001) 是学习率，正则化参数 λ 设置为 1/epochs。因此，正则化值减少了 epochs 的数量增加。

```python
from sklearn.metrics import accuracy_score

## Clip the weights 
index = list(range(10,90))
w1 = np.delete(w1,index)
w2 = np.delete(w2,index)

w1 = w1.reshape(10,1)
w2 = w2.reshape(10,1)
## Extract the test data features 
test_f1 = x_test[:,0]
test_f2 = x_test[:,1]

test_f1 = test_f1.reshape(10,1)
test_f2 = test_f2.reshape(10,1)
## Predict
y_pred = w1 * test_f1 + w2 * test_f2
predictions = []
for val in y_pred:
    if(val > 1):
        predictions.append(1)
    else:
        predictions.append(-1)

print(accuracy_score(y_test,predictions))
```

我们现在裁剪权重，因为测试数据仅包含 10 个数据点。我们从测试数据中提取特征并预测值。我们获得预测并将其与实际值进行比较并打印我们模型的准确性。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230215212448998.png)



还有另一种简单的方法来实现 SVM 算法。我们可以使用Scikit learn库，调用相关函数即可实现SVM模型。代码行数明显减少太少的行数。

```python
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score

clf = SVC(kernel='linear')
clf.fit(x_train,y_train)
y_pred = clf.predict(x_test)
print(accuracy_score(y_test,y_pred))
```



## 总结

支持向量机是一种优雅而强大的算法。明智地使用它:)