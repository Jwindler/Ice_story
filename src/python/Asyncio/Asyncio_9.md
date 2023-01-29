# Python 异步: 当前和正在运行的任务（9）

我们可以反省在 asyncio 事件循环中运行的任务。这可以通过为当前运行的任务和所有正在运行的任务获取一个 asyncio.Task 对象来实现。



## 1. 如何获取当前任务

我们可以通过 asyncio.current_task() 函数获取当前任务。此函数将为当前正在运行的任务返回一个任务对象。

```python
...
# get the current task
task = asyncio.current_task()
```

这将为当前正在运行的任务返回一个任务对象。这可能是：

1. 传递给 asyncio.run() 的主协程。
2. 通过 asyncio.create_task() 在 asyncio 程序中创建和调度的任务。



一个任务可以创建并运行另一个协程（例如，不包含在任务中）。从协程中获取当前任务将为正在运行的任务返回一个 Task 对象，但不会返回当前正在运行的协程。

如果协程或任务需要有关自身的详细信息，例如用于日志记录的任务名称，则获取当前任务会很有帮助。

我们可以探索如何为用于启动 asyncio 程序的主协程获取 Task 实例。下面的示例定义了一个用作程序入口点的协程。它报告一条消息，然后获取当前任务并报告其详细信息。

这是第一个重要的示例，因为它强调所有协程都可以作为异步事件循环中的任务进行访问。

下面列出了完整的示例。

```python
# SuperFastPython.com
# example of getting the current task from the main coroutine
import asyncio
 
# define a main coroutine
async def main():
    # report a message
    print('main coroutine started')
    # get the current task
    task = asyncio.current_task()
    # report its details
    print(task)
 
# start the asyncio program
asyncio.run(main())
```

运行该示例首先创建主协程并使用它来启动 asyncio 程序。main() 协程运行并首先报告一条消息。

然后它检索当前任务，这是一个代表自身的任务对象，即当前正在运行的协程。然后它会报告当前正在运行的任务的详细信息。

我们可以看到该任务具有第一个任务的默认名称“Task-1”，并且正在执行 main() 协程，即当前正在运行的协程。

这突出表明我们可以使用 asyncio.current_task() 函数来访问当前正在运行的协程的任务对象，该对象自动包装在任务对象中。

```python
main coroutine started
<Task pending name='Task-1' coro=<main() running at ...> cb=[_run_until_complete_cb() at ...]>
```



## 2. 如何获取所有任务

我们可能需要访问异步程序中的所有任务。这可能有很多原因，例如：

- 反省程序的当前状态或复杂性。
- 记录所有正在运行的任务的详细信息。
- 查找可以查询或取消的任务。

我们可以通过 asyncio.all_tasks() 函数在 asyncio 程序中获取一组所有已计划和正在运行（尚未完成）的任务。

```python
...
# get all tasks
tasks = asyncio.all_tasks()
```

这将返回 asyncio 程序中所有任务的集合。它是一个集合，因此每个任务只代表一次。

如果出现以下情况，将包括一项任务：

1. 任务已安排但尚未运行。
2. 该任务当前正在运行（例如，但当前已暂停）



该集合还将包括当前正在运行的任务的任务，例如正在执行调用 asyncio.all_tasks() 函数的协程的任务。

另外，回想一下用于启动 asyncio 程序的 asyncio.run() 方法会将提供的协程包装在任务中。这意味着所有任务的集合将包括程序入口点的任务。

我们可以探索在一个 asyncio 程序中有很多任务的情况，然后得到一组所有任务。

在此示例中，我们首先创建 10 个任务，每个任务包装并运行相同的协程。主协程然后获取程序中计划或运行的所有任务的集合并报告它们的详细信息。

下面列出了完整的示例。

```python
# SuperFastPython.com
# example of starting many tasks and getting access to all tasks
import asyncio
 
# coroutine for a task
async def task_coroutine(value):
    # report a message
    print(f'task {value} is running')
    # block for a moment
    await asyncio.sleep(1)
 
# define a main coroutine
async def main():
    # report a message
    print('main coroutine started')
    # start many tasks
    started_tasks = [asyncio.create_task(task_coroutine(i)) for i in range(10)]
    # allow some of the tasks time to start
    await asyncio.sleep(0.1)
    # get all tasks
    tasks = asyncio.all_tasks()
    # report all tasks
    for task in tasks:
        print(f'> {task.get_name()}, {task.get_coro()}')
    # wait for all tasks to complete
    for task in started_tasks:
        await task
 
# start the asyncio program
asyncio.run(main())
```

运行该示例首先创建主协程并使用它来启动 asyncio 程序。main() 协程运行并首先报告一条消息。然后它创建并安排 10 个包装自定义协程的任务。然后 main() 协程会阻塞片刻以允许任务开始运行。任务开始运行，每个任务报告一条消息，然后休眠。

main() 协程恢复并获取程序中所有任务的列表。然后它报告每个的名称和协程。最后，它枚举已创建的任务列表并等待每个任务完成。

这突出表明我们可以获得 asyncio 程序中所有任务的集合，其中包括创建的任务以及代表程序入口点的任务。

```python
main coroutine started
task 0 is running
task 1 is running
task 2 is running
task 3 is running
task 4 is running
task 5 is running
task 6 is running
task 7 is running
task 8 is running
task 9 is running
> Task-9, <coroutine object task_coroutine at 0x10e186e30>
> Task-2, <coroutine object task_coroutine at 0x10e184e40>
> Task-11, <coroutine object task_coroutine at 0x10e186f10>
> Task-7, <coroutine object task_coroutine at 0x10e186d50>
> Task-4, <coroutine object task_coroutine at 0x10e185700>
> Task-10, <coroutine object task_coroutine at 0x10e186ea0>
> Task-8, <coroutine object task_coroutine at 0x10e186dc0>
> Task-5, <coroutine object task_coroutine at 0x10e186ab0>
> Task-1, <coroutine object main at 0x10e1847b0>
> Task-3, <coroutine object task_coroutine at 0x10e184f90>
> Task-6, <coroutine object task_coroutine at 0x10e186ce0>
```

接下来，我们将探讨如何同时运行多个协程。