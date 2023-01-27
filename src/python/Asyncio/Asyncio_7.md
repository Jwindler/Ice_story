# Python 异步: 创建和运行异步任务（7）

您可以从 asyncio 程序中的协程创建任务对象。任务提供独立调度和运行的协程的句柄，并允许查询、取消任务，以及稍后检索结果和异常。异步事件循环管理任务。因此，所有协程都成为事件循环中的任务并作为任务进行管理。

让我们仔细看看 asyncio 任务。



## 1.  什么是异步任务

异步任务是一个调度并独立运行 asyncio 协程的对象。它提供了一个调度协程的句柄，asyncio 程序可以查询并使用它来与协程交互。

任务是从协程创建的。它需要一个协程对象，包装协程，安排它执行，并提供与之交互的方法。任务独立执行。这意味着它被安排在 asyncio 事件循环中，并且无论创建它的协程中发生了什么，它都会执行。这与直接执行协程不同，后者调用者必须等待它完成。

asyncio.Task 类扩展了 asyncio.Future 类，一个实例是可等待的。Future 是一个较低级别的类，代表最终会到达的结果。扩展 Future 类的类通常被称为 Future-like。

因为异步任务是可等待的，这意味着协程可以使用 await 表达式等待任务完成。

```python
...
# wait for a task to be done
await task
```

现在我们知道什么是 asyncio 任务，让我们看看如何创建一个。



## 2. 如何创建任务

使用提供的协程实例创建任务。回想一下协程是使用 async def 表达式定义的，看起来像一个函数。

```python
# define a coroutine
async def task_coroutine():
	# ...
```

任务只能在协程中创建和调度。创建和调度任务有两种主要方式，它们是：

1. 使用高级 API 创建任务（首选）
2. 使用低级 API 创建任务



### 2.1. 高级 API

可以使用 asyncio.create_task() 函数创建任务。asyncio.create_task() 函数接受一个协程实例和一个可选的任务名称，并返回一个 asyncio.Task 实例。

```python
...
# create a coroutine
coro = task_coroutine()
# create a task from a coroutine
task = asyncio.create_task(coro)
```

这可以通过在一行中使用复合语句来实现。

```python
...
# create a task from a coroutine
task = asyncio.create_task(task_coroutine())
```

这将做几件事：

1. 将协程包装在异步任务实例中。
2. 安排任务在当前事件循环中执行。
3. 返回一个任务实例



任务实例可以被丢弃，通过方法与之交互，并由协程等待。这是从 asyncio 程序中的协程创建任务的首选方法。



### 2.2. 低级 API

也可以使用较低级别的 asyncio API 从协程创建任务。

第一种方法是使用 asyncio.ensure_future() 函数。此函数采用任务、未来或类似未来的对象，例如协程，以及可选的用于调度它的循环。如果没有提供循环，它将被安排在当前事件循环中。

如果为这个函数提供了协程，它会为我们包装在一个实例中，然后返回。

```python
...
# create and schedule the task
task = asyncio.ensure_future(task_coroutine())
```

我们可以用来创建和调度任务的另一个低级函数是 loop.create_task() 方法。此函数需要访问特定的事件循环，在该事件循环中将协程作为任务执行。

我们可以通过 asyncio.get_event_loop() 函数获取 asyncio 程序中当前事件循环的实例。然后可以使用它来调用 create_task() 方法来创建一个 Task 实例并安排它执行。

```python
...
# get the current event loop
loop = asyncio.get_event_loop()
# create and schedule the task
task = loop.create_task(task_coroutine())
```



## 3. 任务何时运行？

创建任务后的一个常见问题是它什么时候运行？

虽然我们可以通过 create_task() 函数调度协程作为任务独立运行，但它可能不会立即运行。事实上，直到事件循环有机会运行，任务才会执行。

直到所有其他协程都没有运行并且轮到任务运行时才会发生这种情况。

例如，如果我们有一个 asyncio 程序，其中有一个创建和调度任务的协程，则调度的任务将不会运行，直到创建任务的调用协程被挂起。

如果调用协程选择休眠，选择等待另一个协程或任务，或者选择等待已安排的新任务，则可能会发生这种情况。

```python
...
# create a task from a coroutine
task = asyncio.create_task(task_coroutine())
# await the task, allowing it to run
await task
```

现在我们知道什么是任务以及如何安排它们。