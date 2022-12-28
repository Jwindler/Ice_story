# 利用mAP评估目标检测模型

在[本文](https://www.kdnuggets.com/2021/03/evaluating-object-detection-models-using-mean-average-precision.html "Source")中，我们将了解如何使用 `precision` 和召回率来计算平均精度 (mAP)。`mAP` 将真实边界框与检测到的框进行比较并返回分数。分数越高，模型的检测越准确。

之前我们详细研究了混淆矩阵、模型准确性、精确度和召回率。我们也使用 `Scikit-learn` 库来计算这些指标。现在我们将扩展讨论以了解如何使用精度和召回率来计算 `mAP`。



## 1. 从预测分数到类别标签

在本节中，我们将快速回顾一下如何从预测分数中派生出类标签。鉴于有两个类别，正类和负类，这里是 10 个样本的真实标签。

```python
y_true = ["positive", "negative", "negative", "positive", "positive", "positive", "negative", "positive", "negative", "positive"]
```

当这些样本被输入模型时，它会返回以下预测分数。基于这些分数，我们如何对样本进行分类（即为每个样本分配一个类标签）？

```python
pred_scores = [0.7, 0.3, 0.5, 0.6, 0.55, 0.9, 0.4, 0.2, 0.4, 0.3]
```

为了将分数转换为类别标签，使用了一个阈值。当分数等于或高于阈值时，样本被归为一类。否则，它被归类为其他类别。如果样本的分数高于或等于阈值，则该样本为阳性。否则，它是负面的。下一个代码块将分数转换为阈值为 0.5 的类别标签。

```python
import numpy

pred_scores = [0.7, 0.3, 0.5, 0.6, 0.55, 0.9, 0.4, 0.2, 0.4, 0.3]
y_true = ["positive", "negative", "negative", "positive", "positive", "positive", "negative", "positive", "negative", "positive"]

threshold = 0.5
y_pred = ["positive" if score >= threshold else "negative" for score in pred_scores]
print(y_pred)
```

- `y_pred`

```python
['positive', 'negative', 'positive', 'positive', 'positive', 'positive', 'negative', 'negative', 'negative', 'negative']
```

现在 `y_true` 和 `y_pred` 变量中都提供了真实标签和预测标签。基于这些标签，可以计算混淆矩阵、精度和召回率。

```python
r = numpy.flip(sklearn.metrics.confusion_matrix(y_true, y_pred))
print(r)

precision = sklearn.metrics.precision_score(y_true=y_true, y_pred=y_pred, pos_label="positive")
print(precision)

recall = sklea
```

- 结果

```python
# Confusion Matrix (From Left to Right & Top to Bottom: True Positive, False Negative, False Positive, True Negative)
[[4 2]
 [1 3]]

# Precision = 4/(4+1)
0.8

# Recall = 4/(4+2)
0.6666666666666666
```

在快速回顾了计算准确率和召回率之后，在下一节中我们将讨论创建准确率-召回率曲线。



## 2. PR 曲线

根据第 1 部分给出的精度和召回率的定义，请记住**精度越高，模型将样本分类为正时的置信度就越高。召回率越高，模型正确分类为正的正样本就越多**。

> 当一个模型的召回率高但精度低时，该模型会正确分类大部分正样本，但它有很多误报（即将许多负样本分类为正样本）。当模型具有高精度但召回率低时，模型将样本分类为正样本时是准确的，但它可能仅对部分正样本进行分类。



由于精度和召回率的重要性，精度-召回率曲线显示了不同阈值的精度和召回率值之间的权衡。该曲线有助于选择最佳阈值以最大化两个指标。

创建精确-召回曲线需要一些输入：

1. 真实标签。
2. 样本的预测分数。
3. 将预测分数转换为类别标签的一些阈值。

下面的代码块创建 y_true 列表来保存真实标签，pred_scores 列表用于预测分数，最后是不同阈值的阈值列表。

```python
import numpy

y_true = ["positive", "negative", "negative", "positive", "positive", "positive", "negative", "positive", "negative", "positive", "positive", "positive", "positive", "negative", "negative", "negative"]

pred_scores = [0.7, 0.3, 0.5, 0.6, 0.55, 0.9, 0.4, 0.2, 0.4, 0.3, 0.7, 0.5, 0.8, 0.2, 0.3, 0.35]
```

以下是保存在阈值列表中的阈值。因为有 10 个阈值，所以将创建 10 个精度和召回值。

```py
[0.2, 
 0.25, 
 0.3, 
 0.35, 
 0.4, 
 0.45, 
 0.5, 
 0.55, 
 0.6, 
 0.65]
```

下面是名为 `precision_recall_curve()` 的函数，其接受真实标签、预测分数和阈值。它返回两个代表精度和召回值的等长列表。

```python
import sklearn.metrics

def precision_recall_curve(y_true, pred_scores, thresholds):
    precisions = []
    recalls = []
    
    for threshold in thresholds:
        y_pred = ["positive" if score >= threshold else "negative" for score in pred_scores]

        precision = sklearn.metrics.precision_score(y_true=y_true, y_pred=y_pred, pos_label="positive")
        recall = sklearn.metrics.recall_score(y_true=y_true, y_pred=y_pred, pos_label="positive")
        
        precisions.append(precision)
        recalls.append(recall)

    return precisions, recalls
```

下一段代码在三个先前准备好的列表后调用 `precision_recall_curve()` 函数。它返回精度和召回列表，分别包含精度和召回的所有值。

```python
precisions, recalls = precision_recall_curve(y_true=y_true, 
                                             pred_scores=pred_scores,
                                             thresholds=thresholds)
```

- 以下是精度列表中的返回值。

```python
[0.5625,
 0.5714285714285714,
 0.5714285714285714,
 0.6363636363636364,
 0.7,
 0.875,
 0.875,
 1.0,
 1.0,
 1.0]
```

- 这是召回列表中的值列表。

```python
[1.0,
 0.8888888888888888,
 0.8888888888888888,
 0.7777777777777778,
 0.7777777777777778,
 0.7777777777777778,
 0.7777777777777778,
 0.6666666666666666,
 0.5555555555555556,
 0.4444444444444444]
```

给定两个长度相等的列表，可以在二维图中绘制它们的值，如下所示。

```python
matplotlib.pyplot.plot(recalls, precisions, linewidth=4, color="red")
matplotlib.pyplot.xlabel("Recall", fontsize=12, fontweight='bold')
matplotlib.pyplot.ylabel("Precision", fontsize=12, fontweight='bold')
matplotlib.pyplot.title("Precision-Recall Curve", fontsize=15, fontweight="bold")
matplotlib.pyplot.show()
```

准确率-召回率曲线如下图所示。请注意，随着召回率的增加，精度会降低。原因是当正样本数量增加（高召回率）时，正确分类每个样本的准确率降低（低精度）。这是预料之中的，因为当有很多样本时，模型更有可能预测出错。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221228204751903.png)



准确率-召回率曲线可以很容易地确定准确率和召回率都高的点。根据上图，最好的点是(recall, precision)=(0.778, 0.875)。

使用上图以图形方式确定精度和召回率的最佳值可能有效，因为曲线并不复杂。更好的方法是使用称为 f1 分数的指标，它是根据下一个等式计算的。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221228204905166.png)



f1 指标衡量准确率和召回率之间的平衡。当 f1 的值很高时，这意味着精度和召回率都很高。较低的 f1 分数意味着精确度和召回率之间的失衡更大。

根据前面的例子，f1 是根据下面的代码计算的。根据 f1 列表中的值，最高分是 0.82352941。它是列表中的第 6 个元素（即索引 5）。召回率和精度列表中的第 6 个元素分别为 0.778 和 0.875。相应的阈值为 0.45。

```python
f1 = 2 * ((numpy.array(precisions) * numpy.array(recalls)) / (numpy.array(precisions) + numpy.array(recalls)))
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221228205207969.png)



下图以蓝色显示了与召回率和准确率之间的最佳平衡相对应的点的位置。总之，平衡精度和召回率的最佳阈值是 0.45，此时精度为 0.875，召回率为 0.778。

```python
matplotlib.pyplot.plot(recalls, precisions, linewidth=4, color="red", zorder=0)
matplotlib.pyplot.scatter(recalls[5], precisions[5], zorder=1, linewidth=6)

matplotlib.pyplot.xlabel("Recall", fontsize=12, fontweight='bold')
matplotlib.pyplot.ylabel("Precision", fontsize=12, fontweight='bold')
matplotlib.pyplot.title("Precision-Recall Curve", fontsize=15, fontweight="bold")
matplotlib.pyplot.show()
```

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221228205420667.png)



## 3. AP

平均精度 (AP) 是一种将精度召回曲线汇总为表示所有精度平均值的单个值的方法。根据面等式计算 `AP`。使用遍历所有精度/召回率的循环，计算当前召回率和下一次召回率之间的差异，然后乘以当前精度。换句话说，`AP` 是每个阈值的精度加权和，其中权重是召回率的增加。

![AP](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221228205639175.png)



分别在召回率和准确率列表上附加 0 和 1 很重要。例如，如果召回列表为 0.8、0.60.8、0.6，则应在 0.8、0.6、0.00.8、0.6、0.0 后附加 0。精度列表也是如此，但附加了 1 而不是 0（例如 0.8、0.2、1.00.8、0.2、1.0）。

鉴于召回率和精度都是 `NumPy` 数组，前面的等式根据下面 `Python` 代码建模。

```python
AP = numpy.sum((recalls[:-1] - recalls[1:]) * precisions[:-1])
```

- 这是计算 `AP` 的完整代码。

```python
import numpy
import sklearn.metrics

def precision_recall_curve(y_true, pred_scores, thresholds):
    precisions = []
    recalls = []
    
    for threshold in thresholds:
        y_pred = ["positive" if score >= threshold else "negative" for score in pred_scores]

        precision = sklearn.metrics.precision_score(y_true=y_true, y_pred=y_pred, pos_label="positive")
        recall = sklearn.metrics.recall_score(y_true=y_true, y_pred=y_pred, pos_label="positive")
        
        precisions.append(precision)
        recalls.append(recall)

    return precisions, recalls

y_true = ["positive", "negative", "negative", "positive", "positive", "positive", "negative", "positive", "negative", "positive", "positive", "positive", "positive", "negative", "negative", "negative"]
pred_scores = [0.7, 0.3, 0.5, 0.6, 0.55, 0.9, 0.4, 0.2, 0.4, 0.3, 0.7, 0.5, 0.8, 0.2, 0.3, 0.35]
thresholds=numpy.arange(start=0.2, stop=0.7, step=0.05)

precisions, recalls = precision_recall_curve(y_true=y_true, 
                                             pred_scores=pred_scores, 
                                             thresholds=thresholds)

precisions.append(1)
recalls.append(0)

precisions = numpy.array(precisions)
recalls = numpy.array(recalls)

AP = numpy.sum((recalls[:-1] - recalls[1:]) * precisions[:-1])
print(AP)
```

- 这都是关于平均精度的。以下是计算 `AP` 的步骤摘要：

1. 使用模型生成预测分数。
2. 将预测分数转换为类别标签。
3. 计算混淆矩阵。
4. 计算精度和召回率指标。
5. 创建精确召回曲线。
6. 测量平均精度。



## 4. IoU

要训练目标检测模型，通常有 2 个输入：

1. 图片
2. 图像检测结果的真实框

该模型预测检测到的对象的边界框。预计预测框不会与真实框完全匹配。下图显示了猫的图像。对象的真实框为红色，而预测框为黄色。基于 2 个框的可视化，模型是否做出了高匹配分数的良好预测？

很难主观地评估模型预测。例如，有人可能会得出匹配率为 50% 的结论，而其他人则注意到匹配率为 60%。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221228210606174.png)



更好的替代方法是使用定量测量来对真实框和预测框的匹配程度进行评分。此度量是交并集 (IoU)。 `IoU` 有助于了解一个区域是否有对象。

`IoU` 是根据下面等式计算的，通过将 2 个框之间的交叉区域除以它们的联合区域。 `IoU` 越高，预测越好。

![IoU](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221228210804053.png)



下图显示了具有不同 `IoU` 的 3 个案例。请注意，每个案例顶部的 `IoU` 都是客观测量的，可能与现实略有不同，但它是有道理的。
对于案例 A，预测的黄色框远未与红色真值框对齐，因此 IoU 得分为 0.2（即两个框之间只有 20% 的重叠）。
对于情况 B，2 个框之间的交叉区域更大，但 2 个框仍然没有很好地对齐，因此 IoU 分数为 0.5。
对于案例 C，两个框的坐标非常接近，因此它们的 IoU 为 0.9（即两个框之间有 90% 的重叠）。
请注意，当预测框和真实框之间的重叠率为 0% 时，IoU 为 0.0。当 2 个框彼此 100% 匹配时，IoU 为 1.0。

![ABC](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221228210854341.png)



要计算图像的 IoU，这里有一个名为 `intersection_over_union()` 的函数。它接受以下 2 个参数：

1. `gt_box`：真实边界框。
2. `pred_box`：预测边界框。

它分别计算交集和并集变量中两个框之间的交集和并集。此外，IoU 是在 iou 变量中计算的。它返回所有这 3 个变量。

```python
def intersection_over_union(gt_box, pred_box):
    inter_box_top_left = [max(gt_box[0], pred_box[0]), max(gt_box[1], pred_box[1])]
    inter_box_bottom_right = [min(gt_box[0]+gt_box[2], pred_box[0]+pred_box[2]), min(gt_box[1]+gt_box[3], pred_box[1]+pred_box[3])]

    inter_box_w = inter_box_bottom_right[0] - inter_box_top_left[0]
    inter_box_h = inter_box_bottom_right[1] - inter_box_top_left[1]

    intersection = inter_box_w * inter_box_h
    union = gt_box[2] * gt_box[3] + pred_box[2] * pred_box[3] - intersection
    
    iou = intersection / union

    return iou, intersection, union
```

传递给函数的边界框是一个包含 4 个元素的列表，它们是：

1. 左上角的 x 轴。
2. 左上角的 y 轴。
3. 宽
4. 高

这是汽车图像的真实边界框和预测边界框。

```python
gt_box = [320, 220, 680, 900]
pred_box = [500, 320, 550, 700]
```

图像名为 cat.jpg，下面是在图像上绘制边界框的完整示例。

```python
import imageio
import matplotlib.pyplot
import matplotlib.patches

def intersection_over_union(gt_box, pred_box):
    inter_box_top_left = [max(gt_box[0], pred_box[0]), max(gt_box[1], pred_box[1])]
    inter_box_bottom_right = [min(gt_box[0]+gt_box[2], pred_box[0]+pred_box[2]), min(gt_box[1]+gt_box[3], pred_box[1]+pred_box[3])]

    inter_box_w = inter_box_bottom_right[0] - inter_box_top_left[0]
    inter_box_h = inter_box_bottom_right[1] - inter_box_top_left[1]

    intersection = inter_box_w * inter_box_h
    union = gt_box[2] * gt_box[3] + pred_box[2] * pred_box[3] - intersection
    
    iou = intersection / union

    return iou, intersection, union

im = imageio.imread("cat.jpg")

gt_box = [320, 220, 680, 900]
pred_box = [500, 320, 550, 700]

fig, ax = matplotlib.pyplot.subplots(1)
ax.imshow(im)

gt_rect = matplotlib.patches.Rectangle((gt_box[0], gt_box[1]),
                                       gt_box[2],
                                       gt_box[3],
                                       linewidth=5,
                                       edgecolor='r',
                                       facecolor='none')

pred_rect = matplotlib.patches.Rectangle((pred_box[0], pred_box[1]),
                                         pred_box[2],
                                         pred_box[3],
                                         linewidth=5,
                                         edgecolor=(1, 1, 0),
                                         facecolor='none')
ax.add_patch(gt_rect)
ax.add_patch(pred_rect)

ax.axes.get_xaxis().set_ticks([])
ax.axes.get_yaxis().set_ticks([])
```

下图显示了带有边界框的图像。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221228211421308.png)



要计算 IoU，只需调用 `intersection_over_union()` 函数。基于边界框，IoU 得分为 0.54。

```python
iou, intersect, union = intersection_over_union(gt_box, pred_box)
print(iou, intersect, union)
```

- 结果

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221228211504924.png)



IoU 分数 0.54 意味着真实边界框和预测边界框之间有 54% 的重叠。看着方框，有人可能在视觉上觉得它足以得出模型检测到猫对象的结论。其他人可能会觉得模型还不准确，因为预测框与真实框不太吻合。

为了客观地判断模型是否正确预测了框的位置，使用了一个阈值。如果模型预测 IoU 分数大于或等于阈值的框，则预测框与其中一个真实框之间存在高度重叠。这意味着该模型能够成功检测到一个对象。检测到的区域被归类为阳性（即包含一个对象）。

另一方面，当 IoU 分数小于阈值时，模型做出了错误的预测，因为预测框与真实框不重叠。这意味着检测到的区域被归类为负面（即不包含对象）。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221228211759807.png)



让我们举个例子来阐明 IoU 分数如何帮助将区域分类为对象。假设对象检测模型由下一张图像提供，其中有 2 个目标对象，其真实框为红色，预测框为黄色。
下一个代码读取图像（假设它被命名为 pets.jpg），绘制框，并计算每个对象的 IoU。左侧对象的 IoU 为 0.76，而另一个对象的 IoU 分数为 0.26。

```python
import matplotlib.pyplot
import matplotlib.patches
import imageio

def intersection_over_union(gt_box, pred_box):
    inter_box_top_left = [max(gt_box[0], pred_box[0]), max(gt_box[1], pred_box[1])]
    inter_box_bottom_right = [min(gt_box[0]+gt_box[2], pred_box[0]+pred_box[2]), min(gt_box[1]+gt_box[3], pred_box[1]+pred_box[3])]

    inter_box_w = inter_box_bottom_right[0] - inter_box_top_left[0]
    inter_box_h = inter_box_bottom_right[1] - inter_box_top_left[1]

    intersection = inter_box_w * inter_box_h
    union = gt_box[2] * gt_box[3] + pred_box[2] * pred_box[3] - intersection
    
    iou = intersection / union

    return iou, intersection, union, 

im = imageio.imread("pets.jpg")

gt_box = [10, 130, 370, 350]
pred_box = [30, 100, 370, 350]

iou, intersect, union = intersection_over_union(gt_box, pred_box)
print(iou, intersect, union)

fig, ax = matplotlib.pyplot.subplots(1)
ax.imshow(im)

gt_rect = matplotlib.patches.Rectangle((gt_box[0], gt_box[1]),
                                       gt_box[2],
                                       gt_box[3],
                                       linewidth=5,
                                       edgecolor='r',
                                       facecolor='none')

pred_rect = matplotlib.patches.Rectangle((pred_box[0], pred_box[1]),
                                         pred_box[2],
                                         pred_box[3],
                                         linewidth=5,
                                         edgecolor=(1, 1, 0),
                                         facecolor='none')
ax.add_patch(gt_rect)
ax.add_patch(pred_rect)

gt_box = [645, 130, 310, 320]
pred_box = [500, 60, 310, 320]

iou, intersect, union = intersection_over_union(gt_box, pred_box)
print(iou, intersect, union)

gt_rect = matplotlib.patches.Rectangle((gt_box[0], gt_box[1]),
                                       gt_box[2],
                                       gt_box[3],
                                       linewidth=5,
                                       edgecolor='r',
                                       facecolor='none')

pred_rect = matplotlib.patches.Rectangle((pred_box[0], pred_box[1]),
                                         pred_box[2],
                                         pred_box[3],
                                         linewidth=5,
                                         edgecolor=(1, 1, 0),
                                         facecolor='none')
ax.add_patch(gt_rect)
ax.add_patch(pred_rect)

ax.axes.get_xaxis().set_ticks([])
ax.axes.get_yaxis().set_ticks([])
```

鉴于 IoU 阈值为 0.6，则只有 IoU 分数大于或等于 0.6 的区域被归类为正（即有物体）。因此，IoU 得分为 0.76 的框为正，而另一个 IoU 为 0.26 的框为负。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221228212250295.png)



如果阈值更改为 0.2 而不是 0.6，则两个预测都是`Positive`的。如果阈值为 0.8，则两个预测都是负面的。

作为总结，IoU 分数衡量预测框与真实框的接近程度。它的范围从 0.0 到 1.0，其中 1.0 是最佳结果。当 IoU 大于阈值时，该框被分类为正，因为它围绕着一个对象。否则，它被归类为负面。



## 5. mAP

通常，目标检测模型使用不同的 IoU 阈值进行评估，其中每个阈值可能给出与其他阈值不同的预测。假设模型由一个图像提供，该图像具有分布在 2 个类中的 10 个对象。如何计算mAP？

要计算 mAP，首先要计算每个类的 AP。所有类别的 AP 的平均值是 mAP。

假设使用的数据集只有 2 个类。对于第一类，这里分别是 `y_true` 和 `pred_scores` 变量中的真实标签和预测分数。

```python
y_true = ["positive", "negative", "positive", "negative", "positive", "positive", "positive", "negative", "positive", "negative"]

pred_scores = [0.7, 0.3, 0.5, 0.6, 0.55, 0.9, 0.75, 0.2, 0.8, 0.3]
```

这是第二类的 `y_true` 和 `pred_scores` 变量。

```python
y_true = ["negative", "positive", "positive", "negative", "negative", "positive", "positive", "positive", "negative", "positive"]

pred_scores = [0.32, 0.9, 0.5, 0.1, 0.25, 0.9, 0.55, 0.3, 0.35, 0.85]
```

IoU 阈值列表从 0.2 到 0.9，步长为 0.25。

```python
thresholds = numpy.arange(start=0.2, stop=0.9, step=0.05)
```

要计算一个类的 AP，只需将其 `y_true` 和 `pred_scores` 变量提供给下一个代码。

```python
precisions, recalls = precision_recall_curve(y_true=y_true, 
                                             pred_scores=pred_scores, 
                                             thresholds=thresholds)

matplotlib.pyplot.plot(recalls, precisions, linewidth=4, color="red", zorder=0)

matplotlib.pyplot.xlabel("Recall", fontsize=12, fontweight='bold')
matplotlib.pyplot.ylabel("Precision", fontsize=12, fontweight='bold')
matplotlib.pyplot.title("Precision-Recall Curve", fontsize=15, fontweight="bold")
matplotlib.pyplot.show()

precisions.append(1)
recalls.append(0)

precisions = numpy.array(precisions)
recalls = numpy.array(recalls)

AP = numpy.sum((recalls[:-1] - recalls[1:]) * precisions[:-1])
print(AP)
```

对于第一类，这是它的精确召回曲线。基于此曲线，AP 为 0.949。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221228212932526.png)



第二类的`precision-recall`曲线如下图所示。它的 AP 是 0.958。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221228213007212.png)



基于 2 个类（0.949 和 0.958）的 AP，根据下面等式计算目标检测模型的 `mAP`。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221228213100305.png)



基于此等式，`mAP` 为 0.9535。

```python
mAP = (0.949 + 0.958)/2 = 0.9535
```



## 总结

本教程讨论了如何计算目标检测模型的平均精度 (mAP)。我们首先讨论如何将预测分数转换为类别标签。使用不同的阈值，创建精确召回曲线。从该曲线可以测量平均精度 (AP)。

对于目标检测模型，阈值是对检测到的对象进行评分的 IoU。一旦为数据集中的每个类测量了 AP，就会计算出 mAP。