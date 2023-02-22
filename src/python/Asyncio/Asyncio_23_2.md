# Python 异步: 常见问题 Part_2（23）

本节回答开发人员在 Python 中使用 asyncio 时提出的常见问题。



## 6. 正在运行的任务是否会阻止事件循环退出？

不会！

独立调度和运行的任务不会阻止事件循环退出。如果你的主协程没有其他活动要完成并且有独立的任务在后台运行，你应该检索正在运行的任务并等待它们



## 7. 如何显示正在运行的任务的进度？

我们可以在每个任务上使用 done 回调函数来显示进度。完成回调是我们可以在 asyncio.Task 上注册的函数。

一旦任务完成，它就会被调用，无论是正常还是失败。done 回调函数是一个常规函数，而不是协程，并将与其关联的 asyncio.Task 作为参数。我们可以对所有任务使用相同的回调函数并以通用方式报告进度，例如通过报告消息。

```python
# callback function to show progress of tasks
def progress(task):
    # report progress of the task
    print('.', end='')
```

我们可以在我们发出的每个 asyncio.Task 上注册一个回调函数。这可以通过在每个任务上使用 add_done_callback() 方法并将回调函数的名称传递给它来实现。

```python
...
# add a done callback to a task
task.add_done_callback(progress)
```



## 8. 如何在延迟后运行任务？

我们可以开发一个自定义包装器协程来在延迟后执行目标协程。包装协程可能有两个参数，一个协程和一个以秒为单位的时间。它将在给定的延迟间隔内休眠（以秒为单位），然后等待提供的协程。下面的 delay() 协程实现了这一点。

```python
# coroutine that will start another coroutine after a delay in seconds
async def delay(coro, seconds):
    # suspend for a time limit in seconds
    await asyncio.sleep(seconds)
    # execute the other coroutine
    await coro
```

要使用包装协程，可以创建协程对象并直接等待或作为任务独立执行。例如，调用者可以挂起并调度延迟协程并等待它完成：

```python
...
# execute a coroutine after a delay
await delay(coro, 10)
```

或者，调用者可以安排延迟协程独立运行：

```python
...
# execute a coroutine after a delay independently
_ = asyncio.create_task(delay(coro, 10))
```



## 9. 如何运行后续任务？

asyncio中主要有3种方式来发布后续任务:

1. 从已完成的任务本身安排后续任务。
2. 安排呼叫者的后续任务。
3. 使用完成回调自动安排后续任务。

完成的任务可以发布自己的后续任务。这可能需要检查一些状态以确定是否应该发出后续任务。然后可以通过调用 asyncio.create_task() 来安排任务。

```python
...
# schedule a follow-up task
task = asyncio.create_task(followup_task())
```

任务本身可以选择等待后续任务，也可以让其在后台独立完成。

```python
...
# wait for the follow-up task to complete
await task
```

发出任务的调用者可以选择发出后续任务。例如，当调用者发出第一个任务时，它可能会保留 asyncio.Task 对象。然后它可以检查任务的结果或任务是否成功完成。然后调用者可以决定发出后续任务。它可能会也可能不会直接等待后续任务。

```python
...
# issue and await the first task
task = await asyncio.create_task(task())
# check the result of the task
if task.result():
	# issue the follow-up task
	followup = await asyncio.create_task(followup_task())
```

我们可以使用 done 回调函数自动执行后续任务。例如，发出任务的调用者可以在任务本身上注册一个完成的回调函数。done 回调函数必须将 asyncio.Task 对象作为参数，并且只有在任务完成后才会被调用。然后它可以选择发布后续任务。done回调函数是一个普通的Python函数，不是协程，所以不能等待后续任务

例如，回调函数可能如下所示：

```python
# callback function
def callback(task):
    # schedule and await the follow-up task
    _ = asyncio.create_task(followup())
```

调用者可以发出第一个任务并注册完成的回调函数。

```python
...
# schedule and the task
task = asyncio.create_task(work())
# add the done callback function
task.add_done_callback(callback)
```



## 10. 如何执行阻塞 I/O 或 CPU 绑定函数？

asyncio 模块提供了两种在 asyncio 程序中执行阻塞调用的方法。第一种是使用 asyncio.to_thread() 函数。这是在高级 API 中，供应用程序开发人员使用。asyncio.to_thread() 函数采用要执行的函数名和任何参数。该函数在单独的线程中执行。它返回一个可以作为独立任务等待或安排的协程。

```python
...
# execute a function in a separate thread
await asyncio.to_thread(task)
```

在返回的协程有机会在事件循环中运行之前，任务不会开始执行。asyncio.to_thread() 函数在后台创建一个 ThreadPoolExecutor 来执行阻塞调用。因此，asyncio.to_thread() 函数仅适用于 IO 绑定任务。

另一种方法是使用 loop.run_in_executor() 函数。这是在低级异步 API 中，首先需要访问事件循环，例如通过 asyncio.get_running_loop() 函数。loop.run_in_executor() 函数接受一个执行器和一个要执行的函数。

如果没有为执行器提供，则使用默认执行器，即 ThreadPoolExecutor。loop.run_in_executor() 函数返回一个可等待对象，如果需要可以等待它。任务将立即开始执行，因此返回的可等待对象不需要等待或安排阻塞调用开始执行。

```python
...
# get the event loop
loop = asyncio.get_running_loop()
# execute a function in a separate thread
await loop.run_in_executor(None, task)
```

或者，可以创建一个执行器并将其传递给 loop.run_in_executor() 函数，该函数将在执行器中执行异步调用。

在这种情况下，调用者必须管理执行器，一旦调用者完成它就将其关闭。

```python
...
# create a process pool
with ProcessPoolExecutor as exe:
	# get the event loop
	loop = asyncio.get_running_loop()
	# execute a function in a separate thread
	await loop.run_in_executor(exe, task)
	# process pool is shutdown automatically...
```

这两种方法允许阻塞调用作为异步任务在 asyncio 程序中执行。