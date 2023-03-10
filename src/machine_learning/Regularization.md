# 用 Dropout 正则化对抗 过拟合



过拟合是我们大多数人在训练和使用机器学习模型时已经或最终会遇到的常见挑战。自机器学习诞生以来，研究人员一直在努力对抗过拟合。他们提出的一种技术是 dropout 正则化，其中模型中的神经元被随机移除。在[本文](https://towardsdatascience.com/combating-overfitting-with-dropout-regularization-f721e8712fbe "Source")中，我们将探讨 dropout 正则化的工作原理、如何在您自己的模型中实施它，以及与其他方法相比它的优缺点。



## 1. 简介

### 1.1. 什么是过拟合

过拟合是指模型在其训练数据上过度训练，导致它在新数据上表现不佳。从本质上讲，在模型力求尽可能准确的过程中，它过分关注训练数据集中的细节和噪声。这些属性通常不存在于真实世界的数据中，因此模型往往表现不佳。当模型的参数相对于数据量而言太多时，就会发生过拟合。这可能导致模型过度关注与模型必须开发的一般模式无关的较小细节。例如，假设训练了一个复杂模型（许多参数）来识别图片中是否有马。在这种情况下，它可能会开始关注天空或环境的细节，而不是马本身。这可能发生在：

- 该模型太复杂（具有太多参数）而不利于其自身。
- 模型训练时间过长。
- 训练模型的数据集太小。
- 该模型在相同的数据上进行训练和测试。
- 训练模型的数据集具有重复的特征，使其容易过拟合。



### 1.2. 重要性

过拟合不仅仅是一个简单的烦恼——它会破坏整个模型。它给人一种模型表现良好的错觉，即使它无法对所提供的数据进行适当的概括。

过拟合会产生极其严重的后果，尤其是在人工智能越来越普及的医疗保健等领域。由于过拟合而未经过适当训练或测试的 AI 可能导致错误诊断。



## 2. 什么是Dropout

- Dropout 是一种正则化技术

理想情况下，对抗过拟合的最佳方法是在同一数据集上训练大量不同架构的模型，然后对它们的输出进行平均。这种方法的问题在于它非常耗费资源和时间。虽然相对较小的模型可能负担得起，但可能需要大量时间来训练的大型模型很容易压垮任何人的资源。

Dropout 的工作原理是从输入层或隐藏层中“丢弃”一个神经元。多个神经元从网络中移除，这意味着它们实际上不存在——它们的传入和传出连接也被破坏。这人为地创建了许多更小、更不复杂的网络。这迫使模型不再完全依赖于一个神经元，这意味着它必须使其方法多样化并开发多种方法来实现相同的结果。例如，回到马的例子，如果一个神经元主要负责马的树部分，它的被丢弃将迫使模型更多地关注图像的其他特征。 Dropout 也可以直接应用于输入神经元，这意味着整个特征都从模型中消失了。



- 将 Dropout 应用于神经网络

通过在每一层（包括输入层）中随机丢弃神经元，将 Dropout 应用于神经网络。预定义的丢弃率决定了每个神经元被丢弃的机会。例如，dropout rate 为 0.25 意味着神经元有 25% 的几率被丢弃。在模型训练期间的每个时期都会应用 Dropout。



## 3. 应用

### 3.1. 数据集

让我们从一个可能容易过拟合的数据集开始：

```python
# Columns: has tail, has face, has green grass, tree in background, has blue sky, 3 columns of noise | is a horse image (1) or not (0)
survey = np.array([
 [1, 1, 1, 1, 1, 1], # tail, face, green grass, tree, blue sky | is a horse image
 [1, 1, 1, 1, 1, 1], # tail, face, green grass, tree blue sky | is a horse image
 [0, 0, 0, 0, 0, 0], # no tail, no face, no green grass, no tree, no blue sky | is not a horse image
 [0, 0, 0, 0, 0, 0], # no tail, no face, no green grass, no tree, no blue sky | is not a horse image
])
```

此数据与我们的马及其环境示例相关。我们将图像的特性抽象为一种易于理解的简单格式。可以清楚地看到，数据并不理想，因为其中有马的图像也恰好包含树木、绿草或蓝天——它们可能在同一张照片中，但一个不影响另一个。



### 3.2. 模型

让我们使用 Keras 快速创建一个简单的 MLP：

```python
# Imports
from keras.models import Sequential
from keras.layers import Dense, Dropout
import numpy as np

# Columns: has tail, has face, has green grass, tree in background, has blue sky, 3 columns of noise | is a horse image (1) or not (0)
survey = np.array([
 [1, 1, 1, 1, 1, 1], # tail, face, green grass, tree, blue sky | is a horse image
 [1, 1, 1, 1, 1, 1], # tail, face, green grass, tree blue sky | is a horse image
 [0, 0, 0, 0, 0, 0], # no tail, no face, no green grass, no tree, no blue sky | is not a horse image
 [0, 0, 0, 0, 0, 0], # no tail, no face, no green grass, no tree, no blue sky | is not a horse image
])

# Define the model
model = Sequential([
    Dense(16, input_dim=5, activation='relu'),
    Dense(8, activation='relu'),
    Dense(1, activation='sigmoid')
])

# Compile the model
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])

# Train the model
X = survey[:, :-1]
y = survey[:, -1]
model.fit(X, y, epochs=1000, batch_size=1)

# Test the model on a new example
test_example = np.array([[1, 1, 0, 0, 0]])
prediction = model.predict(test_example)
print(prediction)
```

我强烈建议使用 Python notebook（例如 Jupyter Notebook）来组织代码，这样您就可以快速重新运行单元而无需重新训练模型。沿每个注释拆分代码。

让我们进一步分析我们正在测试模型的数据：

```python
test_example = np.array([[1, 1, 0, 0, 0]])
```

本质上，我们有一张包含马的所有属性的图像，但没有包含在数据中的任何环境因素（绿草、蓝天、树木等）。模型输出：

```sh
0.02694458
```

即使模型有脸和尾巴——我们用它来识别马——也只有 2.7% 的概率确定图像是马。



### 3.3. Dropout

Keras 使实施 dropout 以及其他防止过拟合的方法变得非常简单。我们只需要返回到包含模型层的列表：

```python
# Define the model
model = Sequential([
    Dense(16, input_dim=5, activation='relu'),
    Dense(8, activation='relu'),
    Dense(1, activation='sigmoid')
])
```

并添加一些 dropout 层！

```python
# Define the model
model = Sequential([
    Dense(16, input_dim=5, activation='relu'),
    Dropout(0.5),
    Dense(8, activation='relu'),
    Dropout(0.5),
    Dense(1, activation='sigmoid')
])
```

现在模型输出：

```sh
0.98883545
```

马图像即使不包含环境变量，也有 99% 的把握是马！

Dropout(0.5) 表示上层中的任何神经元都有 50% 的机会被“丢弃”或从存在中移除。通过实施 dropout，我们基本上以资源高效的方式在数百个模型上训练了 MLP。



### 3.4. Dropout Rate

为你的模型找到理想的 Dropout 率的最好方法是通过反复试验——没有万能的方法。从 0.1 或 0.2 左右的低丢失率开始，然后慢慢增加，直到达到所需的精度。使用我们的马 MLP，0.05 的 dropout 导致模型有 16.5% 的置信度图像是马的图像。另一方面，0.95 的 dropout 只是丢弃了太多神经元以使模型无法运行——但仍然达到了 54.1% 的置信度。这些值不适用于此模型，但这确实意味着它们可能适合其他模型。



## 4. 总结

dropout 是机器学习中用于防止过拟合和整体提高模型性能的一种强大技术。它通过从输入层和隐藏层的模型中随机“丢弃”神经元来实现这一点。这允许分类器在一次训练中训练成百上千个独特的模型，防止它过度关注某些特征。