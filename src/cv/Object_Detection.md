# 两阶段目标检测指南：R-CNN、FPN、Mask R-CNN

[Source](https://medium.com/codex/a-guide-to-two-stage-object-detection-r-cnn-fpn-mask-r-cnn-and-more-54c2e168438c "Source")



## 多阶段（Two-stage）物体检测

计算机视觉中最基本和最广泛研究的挑战之一是目标检测。该任务旨在在给定图像中绘制多个对象边界框，这在包括自动驾驶在内的许多领域非常重要。通常，这些目标检测算法可以分为两类：单阶段模型和多阶段模型。在这篇文章中，我们将通过回顾该领域一些最重要的论文，深入探讨用于对象检测的多阶段管道的关键见解。

![](https://s2.loli.net/2023/05/20/ycVE8TFJz6aZ75U.png)



物体检测器的一个分支是基于多阶段模型。源自 R-CNN 的工作，一个模型用于提取对象区域，第二个模型用于分类并进一步细化对象的定位。众所周知，此类方法相对较慢，但非常强大，但最近的进展（例如共享特征）改进了 2 级检测器，使其具有与单级检测器相似的计算成本。这些工作高度依赖于以前的工作，并且大多建立在以前的管道作为基线的基础上。因此，了解两级检测器中的所有主要算法非常重要。

这篇文章的论文选择主要基于调查。



## R-CNN

2014 年的论文提出了基于 CNN 的两阶段检测算法的朴素版本，该算法在以下论文中得到了改进和加速。如上图所述，整个流水线由三个阶段组成：

1. 生成区域建议：模型必须在图像中绘制候选对象，独立于类别。
2. 第二阶段是一个全卷积神经网络，计算每个候选区域的特征。
3. 最后阶段是全连接层，在论文中表示为 SVM。

![](https://s2.loli.net/2023/05/20/MiTVwsEA6cKu91f.png)



Region proposals可以使用多种方法生成，本文选择使用selective search来与之前的工作进行比较。尽管如此，该管道仍与大多数区域提案方法兼容。此处和本演示文稿中提供了选择性搜索的详细说明。

为了总结选择性搜索，将分割算法应用于图像，并根据分割图绘制区域建议（边界框）。分割图被迭代合并，更大的区域建议从细化的地图中提取，如下图所示。此处详细说明了合并和框绘制的工作原理。

![](https://s2.loli.net/2023/05/20/uSJK6IQebaX3OL5.png)



第二阶段和第三阶段一起可以看作是处理裁剪区域提案的传统 CNN。该论文使用 AlexNet 的卷积部分作为第二阶段，而可以使用任何其他 CNN 架构。由于区域建议的大小不同，本文采用最朴素的方式将所有边界框变形并调整为所需大小。

作者还使用经过训练的边界框分类器来进一步细化通过分割进行的边界框估计。另一个完全连接的网络被训练为输入特征图和回归边界框偏移量在 4 个元组 (r, c, h, w) 中表示相对平移和对数尺度宽度/高度缩放因子。该技术在 R-CNN BB 的消融研究中显示出性能提升。

![](https://s2.loli.net/2023/05/22/3svGI5z1hMtgwjk.png)



为了拒绝推理中的重叠区域提议，其中两个或多个边界框指向同一个对象，作者提出了一种贪婪算法，如果该区域与另一个具有更有信心的预测。

由于图像的域更改为扭曲窗口的图像，因此分类器模型在扭曲图像和新标签上进一步训练。在训练分类器时，与地面实况 (GT) 框具有 >0.5 IoU 的区域被认为是该类别，并被训练为输出 GT 框的类别。当该框与任何 GT 框没有显着重叠时，或者当该区域与每个框的 IoU <0.5 时，分类器必须将该区域分类为背景类。为了解决潜在的类不平衡问题，选择了 32 个正区域和 96 个背景区域来形成大小为 128 的小批量。

当 IoU >0.5 的区域被认为完全重叠时，本文认为 0.3<IoU<0.5 的区域部分重叠。这些案例通过提供背景和 GT 框类的混合标签进行特殊处理。

与其他方法相比，R-CNN 的性能优势来自执行自下而上样式选择性搜索的想法，也使用 CNN 来定位对象，以及用于在对象检测数据上微调网络的技术。这项工作结合了经典 CV 和深度学习的工作，以改进目标检测。但是 R-CNN 非常耗时，因为它将 CNN 应用于大约 2,000 个扭曲的选择性搜索区域。

- 总结

1. 为 2 阶段目标检测提出基线管道：生成区域建议并对其进行分类。
2. 使用选择性搜索生成区域提议
3. 分类网络调整区域提案的大小并预测类别概率（包括背景）和边界框细化。



## SPP-Net

该论文建议使用空间金字塔池 (SPP) 层，该层旨在适用于任何图像大小，而无需将它们调整为固定大小，这可能会导致信息丢失和图像失真。卷积，在 CNN 中被描述为特征提取器，并不是限制固定输入大小的卷积，而是输入大小限制是因为完全连接的分类层。

![](https://s2.loli.net/2023/05/20/GuEfhc5nyDtCLOW.png)



因此，作者提出了一个特殊的池化层，将不同大小的特征进行变换，并将其馈送到全连接层，以消除网络的固定大小约束，如上图所述。

基本上，SPP 层应用最大池化各种比例的输出，与图像大小成比例。 SPP 层使用大小与图像大小成比例的空间箱，允许将任何形状的图像映射到单一大小。每个空间 bin 最大池化其区域中的值，并且可以通过此过程保留空间信息。下图对此进行了描述。每个过滤器都使用覆盖一定比例图像的不同大小的池化进行处理，并将结果连接起来。 256 是特征图中过滤器的数量。

![](https://s2.loli.net/2023/05/20/FNcyCYuszHwDKQ5.png)



虽然 SPP 层不是作者提出的，但他们首先考虑在 CNN 中使用 SPP 层。 SPP 具有以下属性：

1. 无论输入大小如何，都生成固定长度的输出
2. 已知对物体变形（正则化）具有鲁棒性
3. 可以从各种尺度（分辨率）中提取信息

![](https://s2.loli.net/2023/05/20/stBRpxCnDLMA9eo.png)



该论文侧重于图像分类，并展示了对象检测的结果作为泛化性能的证明，但在应用于对象检测时具有一些不同于 R-CNN 算法的有趣特性。

SPP-Net的目标检测流水线如上图所示。 CNN 在完整图像上执行一次，并根据选择性搜索检测到的区域裁剪 CNN 的输出特征。 SPP 应用于每个作物，并根据 SPP 层的输出预测类别。这样，卷积层仅应用于图像一次，并且仅应用与检测到的区域数量相对应的较轻的 FC 层。

卷积特征检测器在图像分类任务上进行了预训练，而不是在对象检测上进行进一步训练。分类器 FC 层是根据地面实况窗口单独训练的。尺度不变性是通过两种方法对图像进行预处理来实现的，如本文所述。在对 FC 网络进行微调时，也应用了 R-CNN 的许多技术。

这篇论文的贡献真的很惊人，因为它减少了几个数量级的训练和推理时间，同时由于不必调整图像大小和扭曲图像甚至提高了性能。然而，我怀疑在图像分类上训练的特征图是否真的包含裁剪图像的空间信息。这在使用深度神经网络时可能是一个大问题，因为接受大小会很大，因此可能会限制 SPP-Net 管道使用更深的特征提取器。其他一些损失可用于在对象检测数据集上一起微调特征提取器。

- 概括

1. 建议应用空间金字塔池来为任意输入大小输出固定长度的特征。
2. 改进训练/推理过程以减少 # forward passes 从每个区域的一次传递（每张图像约 2,000 个区域）到整个图像的一次前向传递。



## Fast R-CNN

以前的目标检测算法，即 R-CNN 通常分别学习定位和分类阶段，这使得训练成本更高。此外，这些算法在测试时非常慢，阻碍了实时应用程序。 Fast R-CNN 联合学习检测对象的空间位置并对它们进行分类。

R-CNN 很慢，因为对每个对象提议都进行了前向传递。虽然 SPP-Nets 确实解决了这个问题并在测试时将 R-CNN 加速了 100 倍，但训练是一个多阶段过程，需要许多密集计算步骤，与 R-CNN 相比仅加速了 3 倍。此外，固定的卷积层对网络的准确性造成了限制。

![](https://s2.loli.net/2023/05/22/kNc67EpAus8XiDj.png)



上图说明了 Fast R-CNN 管道。 CNN 处理图像并根据对象提议裁剪特征图。然后，感兴趣区域 (RoI) 池化层提取固定长度的向量，然后通过全连接网络对其进行处理，以预测类别概率并细化边界框。

RoI 池化层是 SPP 层的一个特例，具有一个金字塔层级。 h × w RoI 窗口被划分为大小为 h/H × w/W 的 H × W 网格，最大池化应用于每个网格单元。输出始终是 H × W 形状的向量。 Fast R-CNN 流程与 SPP-Net 管道非常相似，只是稍作修改。

以前在 SPP-Nets 中，通过卷积层反向传播效率低下，因为感受野可能跨越整个图像，这非常大。 Fast R-CNN 通过同时从一张图像中训练多个 RoI 样本作为小批量来解决这个问题。这些特征可以在训练期间共享，从而加快训练速度并免除缓存特征。这个技巧被称为分层抽样。此外，Fast R-CNN 通过多任务损失联合优化分类器和边界框回归器，而不是单独训练。

![](https://s2.loli.net/2023/05/22/I5YNH6CmvXiy1o3.png)



还对 R-CNN 算法进行了一些额外的改进。例如，Fast R-CNN 使用稳健的 L1 损失而不是 L2 损失进行回归。超参数也有修改。该论文还结合了 R-CNN 和 SPP-Net 的技术。论文中提供了详细的解释。 Fast R-CNN 能够达到 S.O.T.A 精度，同时在训练和测试中都快了几个数量级。

- 概括

1. 将 SPP 修改为 RoI 池化
2. 通过从一张图像中采样多个补丁来进行高效训练 -> 仅在卷积层上进行一次前向/反向传递。
3. -> 通过反向传播启用卷积特征提取器的训练



## Faster R-CNN

论文指出，目标提议阶段是实时目标检测的计算瓶颈。作为一种解决方案，Faster R-CNN 实现了与特征提取器网络共享卷积层的区域提议网络 (RPN)，从而引入了计算对象提议的边际成本。管道与 Fast R-CNN 一致，只是对象提议是通过内部训练的 RPN 进行的，如下图所示。

![](https://s2.loli.net/2023/05/22/3wPBILVjqkA9F2x.png)



RPN 模型接收由特征提取器计算的特征图，并通过在特征图上滑动一个小的 CNN 来输出一个对象建议列表。在每个滑动窗口位置，网络预测 k 个参考框（锚点）的对象建议，每个对象建议由 4 个坐标和一个估计对象概率的分数组成。 RPN模型如下图所示。

![](https://s2.loli.net/2023/05/22/ecNdwTmBgko9RhV.png)



RPN 模型与 Fast R-CNN 分类管道分开训练。 Fast R-CNN 模型的训练与原始过程类似，包括以图像为中心的采样策略。一个区别是可以确定 RoI 大小而不是任意大小。因此，k 个边界框回归器，每个负责改进相应锚点类型的回归器都被训练，受益于锚点设计。

在训练 RPN 模型时，基于与地面实况边界框的 IoU，为每个锚点分配一个二进制标签。根据与真实框的 IoU，标签可以是正的、负的或中性的。 RPN 模型在分数和坐标估计上进行训练。本文讨论了通过梯度下降联合训练两个模型的三种方式。该论文使用替代训练来训练网络，其中首先训练 RPN，然后使用在此过程中计算的建议来训练 Fast R-CNN。

- 概括

1. 代替缓慢的选择性搜索，提出 RPN 来训练边界框提议过程。
2. RPN 模型预测对象在锚点上的概率、位置。
3. 比较各种训练方法，以便将 RPN 模型与原始的基于区域的检测网络一起有效训练。



## 特征金字塔网络 (FPN)

特征化图像金字塔（图 a）提供多尺度特征表示，通过支持尺度不变性可以方便地进行对象检测。该模型必须能够检测图像中物体的所有尺度，改变金字塔的层数可以很容易地抵消物体的尺度方差。但它显然需要相当长的时间来计算多层次的特征，并且没有用于 Fast/Faster R-CNN 等管道（图 b）。

![](https://s2.loli.net/2023/05/22/Hb6EWZjthJlYz5w.png)



卷积神经网络本质上计算多尺度特征表示，因为每一层分层计算不同分辨率的特征图。然而，以前利用 CNN 的层次属性制作具有较小计算量的特征化图像金字塔的工作（图 c）是不完整的。映射 CNN 中间特征的问题是特征根据网络的深度自然地传达不同的语义。为了充分利用 CNN 进行多尺度特征表示，重要的是这些层在所有尺度上都具有强语义。

FPN 旨在为高分辨率特征提供丰富的语义，是一种类似 U-Net 的架构。自下而上的路径（红色）是前馈 CNN。每个分辨率表示为一个阶段，并且为每个阶段定义一个金字塔级别。自上而下的路径（蓝色）通过对来自更高金字塔级别的语义更强的特征图进行上采样来产生更高分辨率的特征。直观上，更多的操作可以增强任意尺度特征图的语义，提供丰富的多尺度特征。使用通过横向连接投射的自下而上路径的特征进一步增强了这些特征。

![](https://s2.loli.net/2023/05/20/SUK6R3CxYNMedqP.png)



FPN 管道为生成具有丰富语义内容的多尺度特征图提供了通用解决方案。当应用于 Faster R-CNN 对象检测流水线时，FPN 架构既适用于生成边界框建议的 RPN 网络，也适用于 Fast R-CNN 基于区域的分类器主干。通过替换主干网络并提供 FPN 输出而不是单个特征图，FPN 被采用到 RPN。在应用锚点时，我们在金字塔输入的不同层次上应用锚点的每个尺度。例如{32² , 64² , 128² , 256² , 512²} 大小的锚点，每个用于特征图 {P2, P3, P4, P5, P6}。 Faster R-CNN 检测网络应用于特征图列表之一，根据边界框的大小确定。

- 概括

1. 提出新的 FPN 网络架构来计算语义丰富的多尺度特征表示。
2. 使用 CNN 的中间层作为多尺度特征和图像金字塔，并使用这些特征训练 RPN 和骨干网络。



## Mask R-CNN

![](https://s2.loli.net/2023/05/20/qF2IRbZmaDOSsyM.png)



Mask R-CNN 的提出是为了解决一个稍微不同的实例分割问题。简而言之，这个问题是对象检测和语义分割的结合。如上所示，该任务旨在生成划分对象的像素级边界。

Mask R-CNN 基于 Faster R-CNN 流水线，但每个对象提议有三个输出，而不是两个。附加分支预测 K(# classes) 个二进制对象掩码，用于分割图像中每个类的对象。使用分类分支的结果选择最终要绘制的实例分割图。这称为解耦掩码和类别预测。

全卷积网络 (FCN) 用于从每个 RoI 绘制 m×m 掩码。与绘制边界框不同，生成像素级掩码需要像素级空间信息。所以函数在生成mask分割时在折叠特征之前分支出来，如下图所示。

RoI 是小特征图，由 RoI 池化操作计算，该操作严格地将特征图切割成 bin。这是因为这会在 RoI 和提取的特征之间引入错位，这些错位在分类中被忽略，但会损害像素级掩模，而像素级掩模会在很大程度上受到小平移的影响。提出了一个 RoIAlign 层，平滑了 RoIPool 的硬切片。 RoIAlgin 层基本上是大地图到小地图的双线性插值。结果显示出巨大的性能提升，作者提出了更多证据表明问题出在对齐不一致上。

![](https://s2.loli.net/2023/05/20/GHqREwKfUAhjBn5.png)



为了训练掩码分支，在原始分类和边界框回归损失函数中添加了一个损失项 L_mask。 mask 损失项被计算为具有 k 类的地面真值分割图和第 k 个掩码之间的交叉熵损失。

![](https://s2.loli.net/2023/05/22/TnyzqBWJDaLvfZK.png)



这篇论文不仅实现了高性能的实例分割，而且在常规边界框对象检测和姿态估计等其他任务中也取得了令人惊讶的结果。上表显示了边界框对象检测的结果，其中 Mask R-CNN 优于更快的 R-CNN。 Faster R-CNN，RoIAlgin 显示了在训练期间未使用掩码损失时的结果。结果表明，在使用掩码预测目标进行训练时，对象检测管道可以学习到更通用、更丰富的特征。

- 概括

1. 通过引入掩码分支，提出了基于 Faster R-CNN 的实例分割的通用框架。
2. 通过解决切片中的未对齐问题来修复 RoIPooling 层。
3. 简单但令人惊叹的论文:)



## Cascade R-CNN

![](https://s2.loli.net/2023/05/22/GChc8l2sNfb5u7P.png)



如果 IoU 高于阈值 u，则该补丁被认为是一个类的示例，或者被认为是背景类。当使用松散的 IoU 阈值（如 u=0.5）对数据集进行训练时，边界框预测会变得嘈杂。但是增加 IoU 阈值并不能解决问题，因为用于训练/推理的最佳 IoU 不匹配。它还将显着减少正样本的数量，引入不平衡数据的问题，这在右图中红色图表的低性能中得到了说明。区分“接近但不正确”的边界框很重要，但在以前的工作中没有研究过。

这些图说明了在 u = 0.5、0.6、0.7 的 IoU 阈值上训练的三个检测器。如左图所示，每个模型在不同的 IoU 范围内表现最佳。该论文提供了更多关于为什么单个分类器难以在整体 IoU 水平上表现一致的原因。基于单个检测器对于单个质量水平是最佳的假设，Cascade R-CNN 训练了一系列用增加的 IoU 阈值训练的检测器。

![](https://s2.loli.net/2023/05/20/4tfOUAx3vZk9o5r.png)



在 Faster R-CNN（图 a）中，RPN 网络提供了用于细化框和分类的 RoI。在 Cascade R-CNN 中，一系列头部提供了前一个头部的边界框估计，而不是 RPN 的 RoI，解释为迭代地改进边界框估计（图 b、d）。理论上，下一个头部的输出应该逐步改进边界框位置，但是训练具有小 IoU 阈值的边界框精炼器不会将 IoU 提高到一定值（上图 c）。因此，Cascade R-CNN 被设计为不同专业回归器的级联（图 d）。因此，更深的阶段能够逐步提高到更高的 IoU 阈值，如下面的 IoU 直方图所述。

![](https://s2.loli.net/2023/05/22/PLz8c5y6sMUQYlb.png)



- 概括

1. 指出 IoU 阈值对物体检测的影响，以及简单修改阈值的问题。
2. 观察到不同的模型在不同的 IoU 范围内表现最好。
3. 级联边界框回归器可确保高置信度边界框输出，而不会引入其他问题。



## 总结

我们回顾了多阶段目标检测的主要方法。这些算法的进步速度真是惊人。普通的 R-CNN 算法速度慢且效率低下。高级算法的许多关键见解都基于共享特征（例如 SPP-Net、Fast R-CNN、Mask R-CNN）并支持对先前固定的管道组件（例如 Fast R-CNN、Faster R）进行梯度训练-CNN, Cascade R-CNN) 来有效地学习更丰富的特征。目标检测是计算机视觉中的一个重要领域，多阶段目标检测是目标检测的主流方法。

最近在多阶段目标检测方面的一项工作是 DetectoRS，它提出通过提出递归特征金字塔来改进网络的主干。虽然最近对对象检测的关注已经转向基于 Transformer 的方法，但这些关于多阶段对象检测的论文总体上提供了对深度学习的深刻见解。

