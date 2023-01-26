# Python 异步: 什么是事件循环 ?（6）

asyncio 程序的核心是事件循环。在本节中，我们将花点时间看一下 asyncio 事件循环。



## 1. 什么是 Asyncio 事件循环

事件循环是用于在单个线程中执行协程的环境。事件循环是异步程序的核心。

它做了很多事情，例如：

1. 执行协程。
2. 执行回调。
3. 执行网络输入/输出。
4. 运行子进程。



事件循环是一种常见的设计模式，并且由于在 JavaScript 中的使用而在最近变得非常流行。

事件循环，顾名思义，就是一个循环。它管理一个任务列表（协同程序）并尝试在循环的每次迭代中按顺序推进每个任务，以及执行其他任务，如执行回调和处理 I/O。

“asyncio”模块提供了访问事件循环并与之交互的功能。这不是典型应用程序开发所必需的。

相反，为框架开发人员提供了对事件循环的访问，这些开发人员希望在 asyncio 模块之上构建或为其库启用 asyncio。

asyncio 模块提供了一个用于访问当前事件循环对象的低级 API，以及一套可用于与事件循环交互的方法。

低级 API 适用于将 asyncio 扩展、补充和集成到第三方库中的框架开发人员。我们很少需要与 asyncio 程序中的事件循环交互，而是使用高级 API。

尽管如此，我们还是可以简单地探讨一下如何获取事件循环。



## 2. 如何启动和获取事件循环

我们在 asyncio 应用程序中创建事件循环的典型方法是通过 asyncio.run() 函数。该函数接受一个协程并将执行它直到完成。我们通常将它传递给我们的主协程并从那里运行我们的程序。有用于创建和访问事件循环的低级函数。asyncio.new_event_loop() 函数将创建一个新的事件循环并返回对它的访问。

```python
...
# create and access a new asyncio event loop
loop = asyncio.new_event_loop()
```

我们可以用一个有效的例子来证明这一点。

在下面的示例中，我们将创建一个新的事件循环，然后报告其详细信息。

```python
# SuperFastPython.com
# example of creating an event loop
import asyncio
 
# create and access a new asyncio event loop
loop = asyncio.new_event_loop()
# report defaults of the loop
print(loop)
```

运行示例创建事件循环，然后报告对象的详细信息。我们可以看到，在这种情况下，事件循环的类型为 _UnixSelectorEventLoop 并且没有运行，但也没有关闭。

```python
<_UnixSelectorEventLoop running=False closed=False debug=False>
```

如果 asyncio 事件循环已经在运行，我们可以通过 asyncio.get_running_loop() 函数访问它。

```python
...
# access he running event loop
loop = asyncio.get_running_loop()
```

还有一个用于获取或启动事件循环的函数，称为 asyncio.get_event_loop()，但它在 Python 3.10 中已弃用，不应使用。



## 3. 什么是事件循环对象

事件循环作为 Python 对象实现。事件循环对象定义了事件循环的实现方式，并提供了与循环交互的通用 API，定义在 AbstractEventLoop 类中。不同平台的事件循环有不同的实现。例如，Windows 和基于 Unix 的操作系统将以不同的方式实现事件循环，因为在这些平台上实现非阻塞 I/O 的底层方式不同。

SelectorEventLoop 类型的事件循环是基于 Unix 的操作系统（如 Linux 和 macOS）的默认设置。

ProactorEventLoop 类型的事件循环是 Windows 上的默认设置。

第三方库可能会实现自己的事件循环以针对特定功能进行优化。



## 4. 为什么要访问事件循环

为什么我们要访问 asyncio 程序之外的事件循环？

我们可能希望从正在运行的 asyncio 程序外部访问事件循环的原因有很多。

1. 监控任务的进度。
2. 发布任务并从中获取结果。
3. 解雇并忘记一次性任务。



asyncio 事件循环可以在程序中用作基于协程任务的线程池的替代方案。事件循环也可以嵌入到普通的 asyncio 程序中并根据需要访问。