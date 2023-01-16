# Python: 什么是异步？ （2）

广义上，asyncio 是指使用协程在 Python 中实现异步编程的能力。

具体来说，它指的是两个要素：

1. 在 Python 3.4 中将“asyncio”模块添加到 Python 标准库中。
2. 在 Python 3.5 中向 Python 语言添加了 async/await 表达式。

模块和语言的变化共同促进了支持基于协程的并发、非阻塞 I/O 和异步编程的 Python 程序的开发。

让我们仔细看看 asyncio 的这两个方面，从语言的变化开始。



## 1. 异步支持

Python 语言已更改为通过添加表达式和类型来适应 asyncio。更具体地说，它被更改为支持协程作为一流的概念。反过来，协程是 asyncio 程序中使用的并发单元。协程是一个可以挂起和恢复的函数。

协程可以通过“async def”表达式定义。它可以接受参数并返回一个值，就像函数一样。

```python
# define a coroutine
async def custom_coro():
	# ...
```

调用协程函数会创建一个协程对象，这是一个新的类。它不执行协程函数。

```python
...
# create a coroutine object
coro = custom_coro()
```

协程可以通过 await 表达式执行另一个协程。这会暂停调用者并安排目标执行。

```python
...
# suspend and schedule the target
await custom_coro()
```

异步迭代器是产生可等待对象的迭代器。可以使用“async for”表达式遍历异步迭代器。

```python
...
# traverse an asynchronous iterator
async for item in async_iterator:
	print(item)
```

这不会并行执行 for 循环。相反，执行 for 循环的调用协程将挂起并在内部等待迭代器产生的每个可等待对象。

异步上下文管理器是可以等待进入和退出方法的上下文管理器。“async with”表达式用于创建和使用异步上下文管理器。调用协程将在进入上下文管理器块之前挂起并等待上下文管理器，在离开上下文管理器块时也是如此。

这些是为支持协程而对 Python 语言进行的主要更改的总结。



## 2. 异步模块

“asyncio”模块提供函数和对象，用于使用异步编程范例开发基于协程的程序。具体来说，它支持带有子进程（用于执行命令）和流（用于 TCP 套接字编程）的非阻塞 I/O。

asyncio 模块的核心是事件循环。这是运行基于协程的程序并实现协程之间协作多任务处理的机制。该模块同时提供高级和低级 API。高级 API 是为我们 Python 应用程序开发人员准备的。在大多数情况下，低级 API 适用于框架开发人员，而不是我们。大多数用例都可以使用高级 API 来满足，这些 API 提供实用程序来处理协程、流、同步原语、子进程和队列，以便在协程之间共享数据。较低级别的 API 为高级 API 提供了基础，包括事件循环的内部结构、传输协议、策略等。

现在我们大致了解了 asyncio 是什么，它用于异步编程。