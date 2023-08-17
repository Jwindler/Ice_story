# PyTorch模型性能分析与优化

训练深度学习模型，尤其是大型模型，可能是一项昂贵的支出。我们可以使用的管理这些成本的主要方法之一是性能优化。性能优化是一个迭代过程，我们不断寻找提高应用程序性能的机会，然后利用这些机会。在之前的文章中（例如此处），我们强调了拥有适当工具来进行此分析的重要性。工具的选择可能取决于许多因素，包括训练加速器的类型（例如 GPU、HPU 或其他）和训练框架。



![](https://s2.loli.net/2023/08/16/EYsZPSRHXdQUxir.png)



[本文](https://towardsdatascience.com/pytorch-model-performance-analysis-and-optimization-10c3c5822869 "Source")的重点是在 GPU 上使用 PyTorch 进行训练。更具体地说，我们将重点关注 PyTorch 的内置性能分析器 PyTorch Profiler，以及查看其结果的方法之一，PyTorch Profiler TensorBoard 插件。

这篇文章并不是要取代有关 PyTorch Profiler 的官方 PyTorch 文档或使用 TensorBoard 插件来分析分析器结果。我们的目的是展示如何在日常开发过程中使用这些工具。

一段时间以来，我对 TensorBoard 插件教程的一个部分特别感兴趣。本教程介绍了一个在流行的 Cifar10 数据集上训练的分类模型（基于 Resnet 架构）。接下来演示如何使用 PyTorch Profiler 和 TensorBoard 插件来识别和修复数据加载器中的瓶颈。

![](https://s2.loli.net/2023/08/16/71wVvPUEhTKyW84.png)



如果仔细观察，你会发现优化后的GPU利用率为40.46%。现在没有办法粉饰这一点：这些结果绝对是糟糕的，应该让你彻夜难眠。正如我们过去所扩展的，GPU 是我们训练机器中最昂贵的资源，我们的目标应该是最大化其利用率。 40.46% 的利用率结果通常代表着加速训练和节省成本的重要机会。当然，我们可以做得更好！在这篇博文中，我们将努力做得更好。我们将首先尝试重现官方教程中提供的结果，看看我们是否可以使用相同的工具来进一步提高训练性能。



## 玩具示例

下面的代码块包含 TensorBoard 插件教程定义的训练循环，并进行了两处小修改：

1. 我们使用与本教程中使用的 CIFAR10 数据集具有相同属性和行为的假数据集。
2. 我们初始化 torch.profiler.schedule，将预热标志设置为 3，将重复标志设置为 1。我们发现，预热步骤数量的轻微增加提高了分析结果的稳定性。

```python
import numpy as np
import torch
import torch.nn
import torch.optim
import torch.profiler
import torch.utils.data
import torchvision.datasets
import torchvision.models
import torchvision.transforms as T
from torchvision.datasets.vision import VisionDataset
from PIL import Image

class FakeCIFAR(VisionDataset):
    def __init__(self, transform):
        super().__init__(root=None, transform=transform)
        self.data = np.random.randint(low=0,high=256,size=(10000,32,32,3),dtype=np.uint8)
        self.targets = np.random.randint(low=0,high=10,size=(10000),dtype=np.uint8).tolist()

    def __getitem__(self, index):
        img, target = self.data[index], self.targets[index]
        img = Image.fromarray(img)
        if self.transform is not None:
            img = self.transform(img)
        return img, target

    def __len__(self) -> int:
        return len(self.data)

transform = T.Compose(
    [T.Resize(224),
     T.ToTensor(),
     T.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))])

train_set = FakeCIFAR(transform=transform)
train_loader = torch.utils.data.DataLoader(train_set, batch_size=32, 
                                           shuffle=True)

device = torch.device("cuda:0")
model = torchvision.models.resnet18(weights='IMAGENET1K_V1').cuda(device)
criterion = torch.nn.CrossEntropyLoss().cuda(device)
optimizer = torch.optim.SGD(model.parameters(), lr=0.001, momentum=0.9)
model.train()

# train step
def train(data):
    inputs, labels = data[0].to(device=device), data[1].to(device=device)
    outputs = model(inputs)
    loss = criterion(outputs, labels)
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

# training loop wrapped with profiler object
with torch.profiler.profile(
        schedule=torch.profiler.schedule(wait=1, warmup=4, active=3, repeat=1),
        on_trace_ready=torch.profiler.tensorboard_trace_handler('./log/resnet18'),
        record_shapes=True,
        profile_memory=True,
        with_stack=True
) as prof:
    for step, batch_data in enumerate(train_loader):
        if step >= (1 + 4 + 3) * 1:
            break
        train(batch_data)
        prof.step()  # Need to call this at the end of each step
```



本教程中使用的 GPU 是 Tesla V100-DGXS-32GB。在这篇文章中，我们尝试使用包含 Tesla V100-SXM2–16GB GPU 的 Amazon EC2 p3.2xlarge 实例重现本教程的性能结果并进行改进。尽管它们共享相同的架构，但这两种 GPU 之间存在一些差异。我们使用 AWS PyTorch 2.0 Docker 映像运行训练脚本。 TensorBoard 查看器概述页面中显示的训练脚本的性能结果如下图所示：

![](https://s2.loli.net/2023/08/16/1g8xopOtMlke7A9.png)



我们首先注意到，与教程相反，我们实验中的概述页面（torch-tb-profiler 版本 0.4.1）将三个分析步骤合并为一个。因此，平均总步时间为 80 毫秒，而不是报告的 240 毫秒。这可以在“跟踪”选项卡中清楚地看到（根据我们的经验，该选项卡几乎总是提供更准确的报告），其中每个步骤大约需要 80 毫秒。

![](https://s2.loli.net/2023/08/16/B8ePFStachIUwD6.png)



请注意，我们的起始点为 31.65% GPU 利用率和 80 毫秒的步长时间，与教程中分别介绍的 23.54% 和 132 毫秒的起始点不同。这可能是由于训练环境（包括 GPU 类型和 PyTorch 版本）的差异造成的。我们还注意到，虽然教程基线结果清楚地将性能问题诊断为 DataLoader 中的瓶颈，但我们的结果却并非如此。我们经常发现数据加载瓶颈会在“概览”选项卡中将自己伪装成高比例的“CPU Exec”或“其他”。



## 优化1：多进程数据加载

让我们首先应用本教程中所述的多进程数据加载。由于 Amazon EC2 p3.2xlarge 实例有 8 个 vCPU，我们将 DataLoader 工作线程的数量设置为 8 以获得最大性能：

```python
train_loader = torch.utils.data.DataLoader(train_set, batch_size=32, 
                               shuffle=True, num_workers=8)
```

本次优化的结果如下所示：

![](https://s2.loli.net/2023/08/16/UmL4z5I1SrbDkJW.png)



对单行代码的更改使 GPU 利用率提高了 200% 以上（从 31.65% 增加到 72.81%），并将训练步骤时间减少了一半以上（从 80 毫秒减少到 37 毫秒）。

本教程中的优化过程到此结束。虽然我们的 GPU 利用率 (72.81%) 比教程中的结果 (40.46%) 高很多，但我毫不怀疑，像我们一样，您会发现这些结果仍然非常不令人满意。

个人评论，您可以随意跳过：想象一下，如果 PyTorch 在 GPU 上训练时默认应用多进程数据加载，可以节省多少全球资金！确实，使用多重处理可能会产生一些不需要的副作用。尽管如此，必须有某种形式的自动检测算法可以运行，以排除识别潜在问题场景的存在，并相应地应用此优化。



## 优化2：内存固定

如果我们分析上次实验的 Trace 视图，我们可以看到大量时间（37 毫秒中的 10 毫秒）仍然花费在将训练数据加载到 GPU 上。

![](https://s2.loli.net/2023/08/16/swTBXZVfClGJQ81.png)



为了解决这个问题，我们将应用 PyTorch 推荐的另一个优化来简化数据输入流，即内存固定。使用固定内存可以提高主机到 GPU 数据复制的速度，更重要的是，允许我们使它们异步。这意味着我们可以在 GPU 中准备下一个训练批次，同时在当前批次上运行训练步骤。有关更多详细信息以及内存固定的潜在副作用，请参阅 PyTorch 文档。

此优化需要更改两行代码。首先，我们将 DataLoader 的 pin_memory 标志设置为 True。

```python
train_loader = torch.utils.data.DataLoader(train_set, batch_size=32, 
                          shuffle=True, num_workers=8, pin_memory=True)
```

然后我们将主机到设备的内存传输（在训练函数中）修改为非阻塞：

```python
inputs, labels = data[0].to(device=device, non_blocking=True), \
                 data[1].to(device=device, non_blocking=True)
```

内存固定优化的结果如下所示：

![](https://s2.loli.net/2023/08/16/hCqZGuwV58PoRvU.png)



我们的 GPU 利用率现在达到了可观的 92.37%，并且我们的步数时间进一步减少。但我们仍然可以做得更好。请注意，尽管进行了这种优化，性能报告仍然表明我们花费了大量时间将数据复制到 GPU 中。我们将在下面的步骤 4 中再次讨论这一点。



## 优化3：增加批量大小

对于我们的下一个优化，我们将注意力集中在上一个实验的内存视图上：

![](https://s2.loli.net/2023/08/16/Nxqd7VsLZWhpkTo.png)



该图表显示，在 16 GB 的 GPU 内存中，我们的利用率峰值低于 1 GB。这是资源利用不足的一个极端例子，通常（尽管并非总是）表明有提高性能的机会。控制内存利用率的一种方法是增加批处理大小。在下图中，我们显示了将批处理大小增加到 512（内存利用率增加到 11.3 GB）时的性能结果。

![](https://s2.loli.net/2023/08/16/eHKEsrfWhF64kt2.png)



虽然 GPU 利用率指标没有太大变化，但我们的训练速度显着提高，从每秒 1200 个样本（批量大小 32 为 46 毫秒）到每秒 1584 个样本（批量大小 512 为 324 毫秒）。

注意：与我们之前的优化相反，增加批量大小可能会对训练应用程序的行为产生影响。不同的模型对批量大小的变化表现出不同程度的敏感度。有些可能只需要对优化器设置进行一些调整即可。对于其他人来说，调整到大批量可能会更困难甚至不可能。请参阅上一篇文章，了解大批量训练中涉及的一些挑战。



## 优化4：减少主机到设备的复制

您可能注意到了我们之前的结果中饼图中代表主机到设备数据副本的红色大碍眼。解决这种瓶颈最直接的方法就是看看是否可以减少每批的数据量。请注意，在图像输入的情况下，我们将数据类型从 8 位无符号整数转换为 32 位浮点数，并在执行数据复制之前应用归一化。在下面的代码块中，我们建议对输入数据流进行更改，其中我们延迟数据类型转换和规范化，直到数据位于 GPU 上：

```python
# maintain the image input as an 8-bit uint8 tensor
transform = T.Compose(
    [T.Resize(224),
     T.PILToTensor()
     ])
train_set = FakeCIFAR(transform=transform)
train_loader = torch.utils.data.DataLoader(train_set, batch_size=1024, shuffle=True, num_workers=8, pin_memory=True)

device = torch.device("cuda:0")
model = torch.compile(torchvision.models.resnet18(weights='IMAGENET1K_V1').cuda(device), fullgraph=True)
criterion = torch.nn.CrossEntropyLoss().cuda(device)
optimizer = torch.optim.SGD(model.parameters(), lr=0.001, momentum=0.9)
model.train()

# train step
def train(data):
    inputs, labels = data[0].to(device=device, non_blocking=True), \
                     data[1].to(device=device, non_blocking=True)
    # convert to float32 and normalize
    inputs = (inputs.to(torch.float32) / 255. - 0.5) / 0.5
    outputs = model(inputs)
    loss = criterion(outputs, labels)
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
```

由于这一变化，从 CPU 复制到 GPU 的数据量减少了 4 倍，并且红色碍眼的现象几乎消失了：

![](https://s2.loli.net/2023/08/16/yi6FzNZXehx38O1.png)



我们现在的 GPU 利用率达到新高，达到 97.51%（！！），训练速度达到每秒 1670 个样本！让我们看看我们还能做什么。



## 优化5：将渐变设置为“无”

在这个阶段，我们似乎充分利用了 GPU，但这并不意味着我们不能更有效地利用它。一种流行的优化据说可以减少 GPU 中的内存操作，即在每个训练步骤中将模型参数梯度设置为 None 而不是零。有关此优化的更多详细信息，请参阅 PyTorch 文档。实现此优化所需要做的就是将optimizer.zero_grad调用的set_to_none设置为True：

```python
optimizer.zero_grad(set_to_none=True)
```

在我们的例子中，这种优化并没有以任何有意义的方式提高我们的性能。



## 优化6：自动混合精度

GPU 内核视图显示 GPU 内核处于活动状态的时间量，并且可以成为提高 GPU 利用率的有用资源：

![](https://s2.loli.net/2023/08/16/UAu4CLijtN8O6VS.png)



该报告中最引人注目的细节之一是未使用 GPU Tensor Core。 Tensor Core 可在相对较新的 GPU 架构上使用，是用于矩阵乘法的专用处理单元，可以显着提高 AI 应用程序性能。它们的缺乏使用可能代表着优化的主要机会。

由于 Tensor Core 是专门为混合精度计算而设计的，因此提高其利用率的一种直接方法是修改我们的模型以使用自动混合精度（AMP）。在 AMP 模式下，模型的部分会自动转换为较低精度的 16 位浮点并在 GPU TensorCore 上运行。

重要的是，请注意，AMP 的完整实现可能需要梯度缩放，但我们的演示中并未包含该梯度缩放。在进行调整之前，请务必查看有关混合精度训练的文档。

下面的代码块演示了启用 AMP 所需的训练步骤的修改。

```python
def train(data):
    inputs, labels = data[0].to(device=device, non_blocking=True), \
                     data[1].to(device=device, non_blocking=True)
    inputs = (inputs.to(torch.float32) / 255. - 0.5) / 0.5
    with torch.autocast(device_type='cuda', dtype=torch.float16):
        outputs = model(inputs)
        loss = criterion(outputs, labels)
    # Note - torch.cuda.amp.GradScaler() may be required  
    optimizer.zero_grad(set_to_none=True)
    loss.backward()
    optimizer.step()
```

对 Tensor Core 利用率的影响如下图所示。尽管它继续表明有进一步改进的机会，但仅用一行代码，利用率就从 0% 跃升至 26.3%。

![](https://s2.loli.net/2023/08/16/lXc9tZCNKj5BbY8.png)



除了提高 Tensor Core 利用率之外，使用 AMP 还可以降低 GPU 内存利用率，从而释放更多空间来增加批处理大小。下图捕获了 AMP 优化且批量大小设置为 1024 后的训练性能结果：

![](https://s2.loli.net/2023/08/16/bW3ZQPcxREC4HS1.png)



尽管 GPU 利用率略有下降，但我们的主要吞吐量指标进一步增加了近 50%，从每秒 1670 个样本增加到 2477 个样本。我们正在发挥作用！

注意：降低模型部分的精度可能对其收敛产生有意义的影响。与增加批量大小（见上文）的情况一样，使用混合精度的影响会因模型而异。在某些情况下，AMP 会毫不费力地工作。其他时候，您可能需要更加努力地调整自动缩放器。还有一些时候，您可能需要显式设置模型不同部分的精度类型（即手动混合精度）。



## 优化7：在图形模式下训练

我们将应用的最终优化是模型编译。与默认的 PyTorch 急切执行模式相反，其中每个 PyTorch 操作都“急切”运行，编译 API 将模型转换为中间计算图，然后以最适合底层的方式编译为低级计算内核。

以下代码块演示了应用模型编译所需的更改：

```python
model = torchvision.models.resnet18(weights='IMAGENET1K_V1').cuda(device)
model = torch.compile(model)
```

模型编译优化结果如下所示：

![](https://s2.loli.net/2023/08/16/J3ylANiqazPQR4g.png)



与之前实验中的 2477 个样本相比，模型编译进一步将我们的吞吐量提高到每秒 3268 个样本，性能额外提升了 32% (!!)。

图编译改变训练步骤的方式在 TensorBoard 插件的不同视图中非常明显。例如，内核视图表明使用了新的（融合的）GPU 内核，而跟踪视图（如下所示）显示了与我们之前看到的完全不同的模式。

![](https://s2.loli.net/2023/08/16/WrlSngL2e1vd5Ay.png)



## 总结

在这篇文章中，我们展示了玩具分类模型性能优化的巨大潜力。尽管还有其他性能分析器可供您使用，每种分析器都有其优点和缺点，但我们选择了 PyTorch Profiler 和 TensorBoard 插件，因为它们易于集成。

我们应该强调的是，成功优化的路径将根据训练项目的细节（包括模型架构和训练环境）而有很大差异。在实践中，实现您的目标可能比我们在此介绍的示例更困难。我们描述的一些技术可能对您的表现影响不大，甚至可能使情况变得更糟。我们还注意到，我们选择的精确优化以及我们选择应用它们的顺序有些随意。强烈鼓励您根据项目的具体细节开发自己的工具和技术来实现优化目标。

机器学习工作负载的性能优化有时被视为次要的、非关键的和令人讨厌的。我希望我们已经成功地让您相信，节省开发时间和成本的潜力值得在性能分析和优化方面进行有意义的投资。而且，嘿，您甚至可能会发现它很有趣:)。