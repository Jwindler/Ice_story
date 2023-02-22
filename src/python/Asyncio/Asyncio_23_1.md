# Python 异步: 常见问题 Part_1（23）

本节回答开发人员在 Python 中使用 asyncio 时提出的常见问题。



## 1. 如何停止任务？

我们可以通过 asyncio.Task 对象上的 cancel() 方法取消任务。如果任务被取消，cancel() 方法返回 True，否则返回 False。

```python
...
# cancel the task
was_cancelled = task.cancel()
```

如果任务已经完成，则无法取消，cancel() 方法将返回 False，任务不会处于已取消状态。

下次任务有机会运行时，它将引发 CancelledError 异常。如果 CancelledError 异常未在包装协程内处理，任务将被取消。

否则，如果在包装协程内处理了 CancelledError 异常，任务将不会被取消。cancel() 方法还可以接受一个消息参数，该参数将在 CancelledError 的内容中使用。

我们可以探索如何取消正在运行的任务。

在这个例子中，我们定义了一个任务协程，它报告一条消息然后阻塞片刻。

然后我们定义用作 asyncio 程序入口点的主协程。它报告一条消息，创建并安排任务，然后等待片刻。

然后主协程在运行时恢复和取消任务。它再等一会儿，让任务响应取消请求。然后主协程报告取消任务的请求是否成功。

任务被取消，然后完成。主协程然后在关闭程序之前报告任务的状态是否已取消。

```python
# SuperFastPython.com
# example of canceling a running task
import asyncio
 
# define a coroutine for a task
async def task_coroutine():
    # report a message
    print('executing the task')
    # block for a moment
    await asyncio.sleep(1)
 
# custom coroutine
async def main():
    # report a message
    print('main coroutine started')
    # create and schedule the task
    task = asyncio.create_task(task_coroutine())
    # wait a moment
    await asyncio.sleep(0.1)
    # cancel the task
    was_cancelled = task.cancel()
    # report whether the cancel request was successful
    print(f'was canceled: {was_cancelled}')
    # wait a moment
    await asyncio.sleep(0.1)
    # check the status of the task
    print(f'canceled: {task.cancelled()}')
    # report a final message
    print('main coroutine done')
 
# start the asyncio program
asyncio.run(main())
```

运行该示例会启动 asyncio 事件循环并执行 main() 协程。main() 协程报告一条消息，然后创建并调度任务协程。然后它暂停并等待片刻以允许任务协程开始运行。任务运行，报告消息并休眠一段时间。

main() 协程恢复和取消任务。它报告取消任务的请求已成功。然后它会休眠片刻，让任务响应要取消的请求。

task_coroutine() 恢复并引发 CancelledError 异常，导致任务失败并完成。main()协程恢复并报告任务是否处于取消状态。在这种情况下，确实如此。

此示例突出显示了取消正在运行的任务的正常情况。

```python
main coroutine started
executing the task
was canceled: True
canceled: True
main coroutine done
```



## 2. 如何等待任务完成？

我们可以通过直接等待 asyncio.Task 对象来等待任务完成。

```python
...
# wait for the task to finish
await task
```

我们可以在一行中创建和等待任务。

```python
...
# create and wait for the task to finish
await asyncio.create_task(custom_coro())
```



## 3. 如何从任务中获取返回值？

我们可能需要将协程的值返回给调用者。我们可以通过等待从协程中检索返回值。它假定正在等待的另一个协程返回一个值。

```python
# coroutine that returns a value
async def other_coro():
	return 100
```

等待其他协程将挂起调用协程并安排其他协程执行。一旦另一个协程完成，调用协程将恢复。返回值将从另一个协程传递给调用者。

```python

...
# execute coroutine and retrieve return value
value = await other_coro()
```

协程可以包装在 asyncio.Task 对象中。这有助于独立执行协程而无需当前协程等待它。

这可以使用 asyncio.create_task() 函数来实现。

```python
...
# wrap coroutine in a task and schedule it for execution
task = asyncio.create_task(other_coro())
```

有两种方法可以从 asyncio.Task 中检索返回值，它们是：

- 等待任务。
- 调用结果() 方法。

我们可以等待任务来检索返回值。如果任务已安排或正在运行，则调用者将挂起，直到任务完成并提供返回值。如果任务完成，将立即提供返回值。

```python
...
# get the return value from a task
value = await task
```

与协程不同，我们可以多次等待任务而不会引发错误。

```python
...
# get the return value from a task
value = await task
# get the return value from a task
value = await task
```

我们还可以通过调用 asyncio.Task 对象的 result() 方法来获取任务的返回值。

```python
...
# get the return value from a task
value = task.result()
```

这需要完成任务。如果不是，将引发 InvalidStateError 异常。如果任务被取消，将引发 CancelledError 异常。



## 4. 如何在后台运行任务？

我们可以通过将协程包装在 asyncio.Task 对象中来在后台运行协程。这可以通过调用 asyncio.create_task() 函数并将其传递给协程来实现。

协程将被包装在一个 Task 对象中，并被安排执行。将返回任务对象，调用者不会挂起。

```python
...
# schedule the task for execution
task = asyncio.create_task(other_coroutine())
```

至少在当前协程出于任何原因挂起之前，任务不会开始执行。我们可以通过暂停片刻让任务开始运行来帮助解决问题。这可以通过休眠零秒来实现。

```python
...
# suspend for a moment to allow the task to start running
await asyncio.sleep(0)
```

这将暂停调用者一小会儿，并允许请求运行的机会。这不是必需的，因为调用者可能会在未来某个时间暂停或作为正常执行的一部分终止。一旦调用者没有事情要做，我们也可以直接等待任务。

```python
...
# wait for the task to complete
await task
```



## 5. 如何等待所有后台任务？

我们可以等待 asyncio 程序中的所有独立任务。这可以通过首先通过 asyncio.all_tasks() 函数获取一组所有当前正在运行的任务来实现。

```python
...
# get a set of all running tasks
all_tasks = asyncio.all_tasks()
```

这将返回一个集合，其中包含一个 asyncio.Task 对象，用于当前正在运行的每个任务，包括 main() 协程。

我们不能直接等待这个集合，因为它会永远阻塞，因为它包含当前任务。因此，我们可以获取当前正在运行的任务的 asyncio.Task 对象并将其从集合中删除。

这可以通过首先调用 asyncio.current_task() 方法来获取当前协程的任务，然后通过 remove() 方法将其从集合中删除来实现。

```python
...
# get the current tasks
current_task = asyncio.current_task()
# remove the current task from the list of all tasks
all_tasks.remove(current_task)
```

最后，我们可以等待剩余的任务集。这将挂起调用者，直到集合中的所有任务都完成。

```python
...
# suspend until all tasks are completed
await asyncio.wait(all_tasks)
```

将它们结合在一起，下面添加到 main() 协程末尾的代码片段将等待所有后台任务完成。

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

