# 机器学习模型的生命周期

![](https://s2.loli.net/2023/05/20/VLo5f8IQHbFlnCZ.png)



您的模型如何变化？[Source](https://towardsdatascience.com/the-lifetime-of-a-machine-learning-model-392e1fadf84a "Source")



## 诞生

当我们构建、训练、拟合或估计我们的模型时，这些数字工具就诞生了。这个阶段几乎从拥有分析目标、数据、计算机、算法以及数据科学家现在已经非常了解的其他一切开始。无论您收集什么其他工具，永远不要忘记分析或科学目标，以便您的最终模型有意义并满足特定需求。你的模特什么时候出生的？当您完成培训并将其保存以供就业/部署时，该工具的生命周期就开始了。

![](https://s2.loli.net/2023/05/20/GbOlJxPaAHyQ5cW.png)



这个新生儿有什么前途？这将取决于分析目标，因此我们在构建它时不能忘记这部分。该模型可以服务于预测任务、指标解释或假设情景模拟，以及许多其他选择。这个工具将用于某些事情。可以是简单快速的事情，也可以是复杂、耗时和长期的事情。这种使用将决定该模型的剩余寿命。如果该模型用于参数的一次性解释，那么生命就没有多少了。但是，如果该模型用于预测，并且旨在为具有在线数据收集的系统提供服务，那么生命就在这个新生儿面前。接下来是什么？



## 维护

随着我们继续使用该模型，支持模型训练的数据条件将开始发生变化。就在那个变化的时候，模型也开始经历变化。如果我们在训练时建立一个预测精度高的流失预测模型，那么在不久或遥远的将来，要预测的客户的条件或行为就会开始发生变化。这种变化挑战了我们学习模型的预测性能。当这些变化发生时，我们的模型进入一个新阶段，我们称之为维护。

在维护阶段，我们可能需要新数据。有了它，我们可以更新模型的规格。这与使用另一台机器（例如汽车）并在机器无法正常工作时调整零件没有什么不同。我们不会深入探讨执行模型维护的策略或解决方案，但一般来说，我们的模型需要经历一个调整过程才能使它们恢复到令人满意的性能。

![](https://s2.loli.net/2023/05/20/2SF6iOxdvKmNGUP.png)



维护机器学习模型与重新训练模型并不完全相同。有些模型可能非常简单，以至于用更新的数据重新训练它们也同样简单。这可能是线性架构或具有很少层和神经元的网络的情况。但是当模型如此复杂且具有大型和深层架构时，维护阶段需要比重新训练模型的费用简单得多。这是当今机器学习领域最重要的主题之一，因为这些工具非常强大，但从长远来看维护起来却非常昂贵。

调整或更新模型后，它就可以恢复使用了，因此模型正在服务的任何进程都可以继续使用更新后的版本。我们的机器可以继续使用。尽管如此，这台机器已经发生了变化。如果你愿意的话，它已经被使用、消耗，并且已经转变为与原始状态略有不同的东西。就像铅笔一样，我们的模型会遇到我们需要削尖它们的尖端以保护它们以便我们可以继续使用它们的时刻。



## 迁移

在机器学习的道路上，我们可能需要走一个出口：转移。当我第一次看到有人调换他们的汽车轮胎在结冰的路上行驶时，我曾经访问过令人惊叹的冰岛。然后当他们回到城市时，他们又换回了普通轮胎。当我开始研究迁移学习时，这个概念变得如此清晰，同时还记得冰岛汽车轮胎的转换。当新的环境/领域开始发挥作用时，我们的模型会进入一个称为迁移的新阶段。

![](https://s2.loli.net/2023/05/20/3lowxuNWryaBkqs.png)



正如同一辆汽车可以通过更换轮胎而无需购买另一辆单独的汽车来适应不同的地面一样，我们可以添加或调整我们模型的某些部分以服务于新领域的新目的，而无需构建新模型。迁移学习是机器学习文献中的另一个研究子领域，旨在优化模型的调整以简化新环境下训练模型的工作。流行的例子是图像识别模型。我们用某些类别的图像训练它们，然后其他人转移这些模型以识别新类别的图像。许多企业现在使用 RegNet、VGG、Inception 或 AlexNet 等模型来调整它们以满足自己的需求。

当我们转移一个模型时，在某种程度上，一个新模型诞生了，它有自己的生命周期，与原来的模型分开。它将像原始模型一样需要维护。有了这个，我们已经从拥有一个初始实体到可能创建一整套模型。毫无疑问，这些数字工具背后确实存在生命周期。

![](https://s2.loli.net/2023/05/20/M3pQghBEfSlzn5O.png)



## 我们的模型会死吗？

简短的回答是：是的。例如，当它们的分析性能在系统上不尽如人意，或者当它们变得如此庞大和如此不同以至于原始模型已成为过去时，它们确实可以停止存在。正如我们在开始时所说，岩石、铅笔和汽车在某个时候都会停止存在。在这方面，模型与这些东西没有什么不同。

尽管该模型可能会灭绝，但直到今天，对它们何时达到这一点的问题的答案是我们在机器学习研究社区中想要回答的最大问题。监控机器学习和模型维护性能的许多发展都与模型何时不再起作用的问题有关。

这个答案不是微不足道的原因之一是因为我们不断需要标签来量化性能的满意度。但机器和统计学习最大的悖论恰恰是标签不可用，而我们构建这些工具来预测它们。另一个原因是，定义性能变化的接受限度可能非常主观。虽然科学家可以提出一些限制，但企业可能有不同的容忍度。

以下是数据科学家在回答这个问题（当前未解决的问题）时也可以考虑的一些要点：

- 训练数据是否过时？ （什么是“太过时了”）
- 当前版本与模型的原始版本有多相似？ （什么是“相似”？）
- 输入特征的可变性和与目标变量的关系是否完全漂移了？ （协变量和概念漂移，机器学习维护研究中的两个最大课题）。
- 部署模型的物理进程是否还在使用？如果物理基础设施不再支持模型的部署，这无疑标志着其生命周期的结束。

不再为模特而活并不一定是消极的事情，更像是她们的一条进化之路。我们需要了解它的生命周期，以使我们的物理和数字系统保持最新状态并具有令人满意的性能。