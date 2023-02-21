# Python 异步: 常见错误（22）

本节举例说明开发人员在 Python 中使用 asyncio 时遇到的一般错误。



## 1. 尝试通过调用协程来运行协程

asyncio 初学者遇到的最常见错误是像调用函数一样调用协程。

例如，我们可以使用“async def”表达式定义协程：

```python
# custom coroutine
async def custom_coro():
	print('hi there')
```

然后初学者将尝试像函数一样调用这个协程，并期望打印消息被报告。

```python
...
# error attempt at calling a coroutine like a function
custom_coro()
```

像调用函数一样调用协程不会执行协程体。相反，它将创建一个协程对象。然后可以在 asyncio 运行时中等待该对象，例如事件循环。

我们可以使用 asyncio.run() 函数启动事件循环来运行协程。

```python
...
# run a coroutine
asyncio.run(custom_coro())
```

或者，我们可以暂停当前协程并使用“await”表达式调度另一个协程。

```python
...
# schedule a coroutine
await custom_coro()
```



## 2. 不要让协程在事件循环中运行

如果协程未运行，您将收到如下运行时警告：

```python
sys:1: RuntimeWarning: coroutine 'custom_coro' was never awaited
```

如果您创建协程对象但不安排它在 asyncio 事件循环中执行，就会发生这种情况。

例如，您可能会尝试从常规 Python 程序中调用协程：

```python
...
# attempt to call the coroutine
custom_coro()
```

这不会调用协程。相反，它将创建一个协程对象。

```python
...
# create a coroutine object
coro = custom_coro()
```

如果您不允许此协程运行，您将收到运行时错误。正如我们在上一节中看到的，您可以通过启动异步事件循环并将协程对象传递给它来让协程运行。

```python
...
# create a coroutine object
coro = custom_coro()
# run a coroutine
asyncio.run(coro)
```

或者，在复合语句的一行中：

```python
...
# run a coroutine
asyncio.run(custom_coro())
```

如果您在 asyncio 程序中遇到此错误，那是因为您已经创建了一个协程并且没有安排它执行。

这可以使用 await 表达式来实现。

```python
...
# create a coroutine object
coro = custom_coro()
# suspend and allow the other coroutine to run
await coro
```

或者，您可以安排它作为任务独立运行。

```python
...
# create a coroutine object
coro = custom_coro()
# schedule the coro to run as a task interdependently
task = asyncio.create_task(coro)
```



## 3. 使用低级 Asyncio API

初学者的一个大问题是他们使用了错误的 asyncio API。这很常见，原因有很多。

- API 在最新版本的 Python 中发生了很大变化。
- API 文档页面让事情变得混乱，显示了两个 API。
- Web 上其他地方的示例混合使用不同的 API。

使用错误的 API 会使事情变得更冗长（例如更多代码）、更困难并且更难理解。

Asyncio 提供了两个 API。

- 面向应用程序开发人员的高级 API（我们）
- 面向框架和库开发人员（不是我们）的低级 API

较低级别的 API 为高级 API 提供了基础，包括事件循环的内部结构、传输协议、策略等。我们应该几乎总是坚持使用高级 API。我们在入门时绝对必须坚持使用高级 API。我们有时可能会深入研究低级 API 以实现特定结果。

如果您开始获取事件循环的句柄或使用“循环”变量来做事，那您就错了。通过高级 API 驱动 asyncio 一段时间。开发一些程序。熟悉异步编程和随意运行协程。



## 4. 过早退出主协程

asyncio 程序的一个主要混淆点是没有给任务足够的时间来完成。我们可以通过 asyncio.create_task() 方法安排许多协程在 asyncio 程序中独立运行。

主协程是 asyncio 程序的入口点，然后可以继续执行其他活动。如果主协程退出，则 asyncio 程序将终止。

即使有一个或多个协程作为任务独立运行，程序也会终止。这会让你措手不及。

您可以发出许多任务，然后让主协程恢复，期望所有发出的任务都在自己的时间内完成。相反，如果主协程没有其他事情可做，它应该等待剩余的任务。

这可以通过首先通过 asyncio.all_tasks() 函数获取一组所有正在运行的任务，从该集合中删除自身，然后通过 asyncio.wait() 函数等待剩余的任务来实现。

```python
...
# get a set of all running tasks
all_tasks = asyncio.all_tasks()
# get the current tasks
current_task = asyncio.current_task()
# remove the current task from the list of all tasks
all_tasks.remove(current_task)
# suspend until all tasks are completed
await asyncio.wait(all_tasks)
```



## 5. 假设竞争条件和死锁是不可能的

并发编程具有特定于并发的故障模式的风险。这包括竞争条件和死锁等问题。

竞争条件涉及两个或多个并发单元同时执行相同的临界区并使资源或数据处于不一致或意外状态。这可能会导致数据损坏和数据丢失。

死锁是指并发单元等待永远不会发生的情况，例如资源变得可用。许多 Python 开发人员认为，使用 asyncio 中的协程不可能出现这些问题。

原因是任何时候在事件循环中只能运行一个协程。的确，一次只能运行一个协程。

问题是，协程可以挂起和恢复，并且可以在使用共享资源或共享变量时这样做。如果不保护关键部分，异步程序中可能会出现竞争条件。如果不仔细管理同步原语，就会发生死锁。

因此，创建 asyncio 程序以确保协程安全是很重要的，协程安全是一个类似于线程安全和进程安全的概念，适用于协程。