# 提高CV模型训练性能的 9 个技巧



## 导读

[本文](https://www.kaggle.com/competitions/happy-whale-and-dolphin/discussion/310105 "Source") 主要想分享一些可能有助于提高计算机视觉任务模型训练速度和准确性的一般技巧或建议，这些建议是通过课程、阅读顶级文章或论文学习所得来的。



## 1. 分辨率

- 从较小的分辨率开始：

前两个技巧侧重于实现更快的模型——您尝试的想法越多，获得更好结果的机会就越大。为了更快地迭代，我们需要从“小”着手以减少我们的训练时间：

Ayush 创建了一个指向所有共享数据集的[数据集线程](https://www.kaggle.com/c/happy-whale-and-dolphin/discussion/309691 "Datasets thread")。从较小的数据集尺寸或分辨率开始可以让您更快地迭代。

如果您使用较小的 GPU 内存，那么可以通过增加 `batch_size` 加快迭代速度。一旦您对自己的想法充满信心并看到了效果提升，您就可以扩展到更大的图像尺寸或分辨率。



## 2. 数据集

- 从数据子集开始：

您应该从少量的数据集或示例开始，并在那里验证您的训练模型。

例如：训练 10 个 `classes`，检查它是否提高了 CV ->  提交

扩展到 20 个`classes`，检查 CV，然后再次提交

如果一切顺利，就在完整数据集上进行训练。



## 3. 精度

- 使用 `FP16` 或半精度训练：

`NVIDIA GPU` 具有 `Tensor-Cores`，在使用“半精度”张量时可提供巨大的加速。这里写了一篇更详细的[博客](https://hackernoon.com/rtx-2080ti-vs-gtx-1080ti-fastai-mixed-precision-training-comparisons-on-cifar-100-761d8f615d7f "半精度训练")，尝试使用 `fp_16` 训练来观察任何 `GPU`（和 `TPU`！）上的加速情况。



## 4. TPU

- 使用 `TPU`：

`Kaggle` 每周提供 20 小时的 `TPU`。 `TPU` 有 8 个核心，这允许您的 `batch_sizes` 是 8 的倍数。这允许更快的训练和更快的迭代。

注意：最近发现了 `Hugging Face Accelerate`，它声称可以在 `TPU` 上使用 `PyTorch` 为您提供简单的工作流程。



## 5. 渐进式

- 渐进式调整大小：

`IIRC` 在 `Efficientnet` 论文中被介绍，也在 `fastai` 课程中教授。

`Chris Deotte` 发表了一篇关于 `CNN` 输入图像大小的[文章](https://www.kaggle.com/c/siim-isic-melanoma-classification/discussion/160147 "post")。这个[博客](https://medium.com/analytics-vidhya/novel-techniques-to-win-an-image-classification-hackathon-part-2-e33bf0ad5fe6 "blog")教你渐进调整大小在 `fastai` 中是如何工作的。长话短说：

- 训练模型尺寸：小
- 保存权重并在更大的图像尺寸上重新训练模型
- 再次保存权重并重新训练最终图像大小

这个过程将会获得更快的收敛和更好的性能。



## 6. `Depthwise Convs`

- 使用 `Depthwise Convs` 而不是 `Regular Convs`：

这个[概念](https://paperswithcode.com/method/depthwise-convolution "concept")首先是在 `MobileNet` 论文中引入的，最近与 `ConvNext` 架构相关的讨论中它再次出现。 `Depthwise Convolutions` 具有更少的 `filters`，因此训练速度更快。

请参阅[此处](https://discuss.pytorch.org/t/how-to-modify-a-conv2d-to-depthwise-separable-convolution/15843/4 "Depthwise Convolutions")以获取有关使其在 `PyTorch` 中运行的一些提示



## 7. 学习率

- 在模型训练期间更改 `learning_rate`：

慢的 lr 需要太长的时间，而快的 lr 可能无法帮助你的模型收敛，使用这个逻辑，我们应该使用动态学习率。

我建议使用 `fastai` 及其 `fine_tune()` 或 `fit_one_cycle()` 函数。有关更多详细信息，请参见[此处](https://forums.fast.ai/t/fine-tune-vs-fit-one-cycle/66029 "学习率")。



## 8. 热身

- 从论文 [Bag of Tricks](https://arxiv.org/pdf/1812.01187.pdf "Bag of Tricks")中，使用 LR 预热是亮点之一：

当你开始训练一个模型时，它具有更多的“随机性”，因为它刚刚开始学习特征，因此首先从较小的 `learning_rate` 开始允许它选择细节，然后你可以在“预热”后将其增加到预期的`schedule`。



## 9. 图像增强

`NNs` 受益于更多数据。图像中的微小变化确实可以帮助模型提高对图像内部特征的理解。使用正确的图像增强可以真正帮助您的模型。此外，在训练模型时可视化结果，以确保它们了解的是特征而不是背景！