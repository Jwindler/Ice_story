# 模型性能分析：ROC 分析和 AUC

[本文](https://towardsdatascience.com/roc-analysis-and-the-auc-area-under-the-curve-404803b694b9 "Source")将介绍模型性能分析的两个方法：ROC & AUC。

`ROC` 分析和曲线下面积 (AUC) 是数据科学中广泛使用的工具，借鉴了信号处理，用于评估不同参数化下模型的质量，或比较两个或多个模型的性能。

传统的性能指标，如准确率和召回率，在很大程度上依赖于正样本的观察。因此，ROC 和 AUC 使用真阳性率和假阳性率来评估质量，同时考虑到正面和负面观察结果。

从分解问题到使用机器学习解决问题的过程有多个步骤。它涉及数据收集、清理和特征工程、构建模型，最后是，评估模型性能。

当您评估模型的质量时，通常会使用精度和召回率等指标，也分别称为数据挖掘领域的置信度和灵敏度。

这些指标将预测值与通常来自保留集的实际观察值进行比较，使用混淆矩阵进行可视化。

![confusion matrix](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221229203021377.png)



让我们首先关注精度，也称为阳性预测值。使用混淆矩阵，您可以将 Precision 构建为所有真实阳性与所有预测阳性的比率。

![Precision](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221229203141883.png)



召回率，也称为真阳性率，表示真阳性与观察到的和预测的所有阳性的比率。

![Recall](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221229203234339.png)



使用混淆矩阵中的不同观察集来描述 `Precision` 和 `Recall`，您可以开始了解这些指标如何提供模型性能的视图。

值得注意的是 Precision 和 Recall 只关注正例和预测，而不考虑任何负例。此外，他们不会将模型的性能与中值场景进行比较，中值场景只是随机猜测。



## 1. ROC 曲线

ROC 作为汇总工具，用于可视化 Precision 和 Recall 之间的权衡。ROC 分析使用 ROC 曲线来确定二进制信号的值有多少被噪声污染，即随机性。它为连续预测器提供了一系列操作点的灵敏度和特异性摘要。ROC 曲线是通过绘制 x 轴上的假阳性率与 y 轴上的真阳性率来获得的。

由于真阳性率是检测信号的概率，而假阳性率是误报的概率，因此 ROC 分析也广泛用于医学研究，以确定可靠地检测疾病或其他行为的阈值。

![ROC](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221229204500119.png)



一个完美的模型将具有等于 1 的误报率和真阳性率，因此它将是 ROC 图左上角的单个操作点。而最差的可能模型将在 ROC 图的左下角有一个单一的操作点，其中误报率等于 1，真阳性率等于 0。

随机猜测模型有 50% 的机会正确预测结果，因此假阳性率将始终等于真阳性率。这就是为什么图中有一条对角线，代表检测信号与噪声的概率为 50/50。



## 2. AUC 面积

要全面分析 ROC 曲线并将模型的性能与其他几个模型进行比较，您实际上需要计算曲线下面积 (AUC)，在文献中也称为 c 统计量。曲线下面积 (AUC) 的值介于 0 和 1 之间，因为曲线绘制在 1x1 网格上，并且与信号理论平行，它是信号可检测性的度量。

这是一个非常有用的统计数据，因为它可以让我们了解模型对真实观察结果和错误观察结果的排名有多好。它实际上是 Wilcoxon-Mann-Whitney 秩和检验的归一化版本，它检验零假设，其中两个有序测量样本是从单个分布 中抽取的。

要绘制 ROC 曲线并计算曲线下面积 (AUC)，您决定使用 SckitLearn 的 RocCurveDisplay 方法并将多层感知器与随机森林模型进行比较，以尝试解决相同的分类任务。

```python
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, RocCurveDisplay

def plot_roc(model, test_features, test_targets):
    """
    Plotting the ROC curve for a given Model and the ROC curve for a Random Forests Models
    """

    # comparing the given model with a Random Forests model
    random_forests_model = RandomForestClassifier(random_state=42)
    random_forests_model.fit(train_features, train_targets)

    rfc_disp = RocCurveDisplay.from_estimator(random_forests_model, test_features, test_targets)
    model_disp = RocCurveDisplay.from_estimator(model, test_features, test_targets, ax=rfc_disp.ax_)
    model_disp.figure_.suptitle("ROC curve: Multilayer Perceptron vs Random Forests")

    plt.show()

# using perceptron model as input
plot_roc(ml_percetron_model, test_features, test_targets)
```

![ROC&AUC](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221229205406041.png)



