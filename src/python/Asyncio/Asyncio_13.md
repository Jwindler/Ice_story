# Python 异步: 保护任务免于取消（13）

Asyncio 任务可以通过调用它们的 cancel() 方法来取消。我们可以通过将任务包装在对 asyncio.shield() 的调用中来保护任务不被取消。

让我们仔细看看。



##  1. 什么是 Asyncio shield()

asyncio.shield() 函数在 Future 中包装了一个可等待对象，它将吸收要取消的请求。

这意味着被屏蔽的未来可以传递给可能尝试取消它的任务，并且取消请求看起来像是成功的，除了被屏蔽的任务或协程将继续运行。

它可能在 asyncio 程序中很有用，其中某些任务可以取消，但其他任务（可能具有更高优先级）则不能。

它也可能在某些任务可以安全取消的程序中很有用，例如那些在设计时考虑了 asyncio 的任务，而其他任务则不能安全终止，因此必须避免取消。

现在我们知道了 asyncio.shield() 是什么，让我们看看如何使用它。



## 2. 如何使用 Asyncio shield()

asyncio.shield() 函数将保护另一个任务或协程不被取消。它以一个可等待对象作为参数并返回一个 asyncio.Future 对象。

然后可以直接等待 Future 对象或将其传递给另一个任务或协程。

```python
...
# shield a task from cancellation
shielded = asyncio.shield(task)
# await the shielded task
await shielded
```

返回的 Future 可以通过调用 cancel() 方法取消。

如果内部任务正在运行，请求将被报告为成功。

```python
...
# cancel a shielded task
was_canceld = shielded.cancel()
```

任何等待 Future 对象的协程都会引发 asyncio.CancelledError，这可能需要处理。

```python
...
try:
	# await the shielded task
	await asyncio.shield(task)
except asyncio.CancelledError:
	# ...
```

重要的是，对 Future 对象的取消请求不会传播到内部任务。这意味着取消请求被护盾吸收了。

```python
...
# create a task
task = asyncio.create_task(coro())
# create a shield
shield = asyncio.shield(task)
# cancel the shield (does not cancel the task)
shield.cancel()
```

如果协程被提供给 asyncio.shield() 函数，它将被包装在 asyncio.Task() 中并立即调度。

这意味着不需要等待屏蔽来让内部协程运行。

如果被屏蔽的任务被取消，取消请求将向上传播到屏蔽，屏蔽也将被取消。

```python
...
# create a task
task = asyncio.create_task(coro())
# create a shield
shield = asyncio.shield(task)
# cancel the task (also cancels the shield)
task.cancel()
```

现在我们知道如何使用 asyncio.shield() 函数，让我们看一些有效的例子。



## 3. 示例

我们可以探索如何使用 asyncio.shield() 来保护任务不被取消。

在这个例子中，我们定义了一个简单的协程任务，它接受一个整数参数，休眠一秒钟，然后返回参数。然后可以创建协程并将其安排为任务。

我们可以定义第二个协程，它接受一个任务，休眠几分之一秒，然后取消提供的任务。

在主协程中，我们可以屏蔽第一个任务，然后将其传递给第二个任务，然后等待被屏蔽的任务。

期望是屏蔽将被取消并保持内部任务完好无损。取消将中断主协程。我们可以在程序结束时检查内部任务的状态，我们希望它已经正常完成，而不管屏蔽上的取消请求如何。

```python
# SuperFastPython.com
# example of using asyncio shield to protect a task from cancellation
import asyncio
 
# define a simple asynchronous
async def simple_task(number):
    # block for a moment
    await asyncio.sleep(1)
    # return the argument
    return number
 
# cancel the given task after a moment
async def cancel_task(task):
    # block for a moment
    await asyncio.sleep(0.2)
    # cancel the task
    was_cancelled = task.cancel()
    print(f'cancelled: {was_cancelled}')
 
# define a simple coroutine
async def main():
    # create the coroutine
    coro = simple_task(1)
    # create a task
    task = asyncio.create_task(coro)
    # created the shielded task
    shielded = asyncio.shield(task)
    # create the task to cancel the previous task
    asyncio.create_task(cancel_task(shielded))
    # handle cancellation
    try:
        # await the shielded task
        result = await shielded
        # report the result
        print(f'>got: {result}')
    except asyncio.CancelledError:
        print('shielded was cancelled')
    # wait a moment
    await asyncio.sleep(1)
    # report the details of the tasks
    print(f'shielded: {shielded}')
    print(f'task: {task}')
 
# start
asyncio.run(main())
```

运行示例首先创建 main() 协程并将其用作应用程序的入口点。创建任务协程，然后将其包装并安排在任务中。然后该任务就不会被取消。

然后将屏蔽的任务传递给 cancel_task() 协程，该协程包装在任务中并进行调度。主协程然后等待受保护的任务，该任务需要 CancelledError 异常。

任务运行片刻然后休眠。取消任务运行片刻，休眠，恢复然后取消屏蔽任务。取消请求报告它已成功。

这会在受保护的 Future 中引发 CancelledError 异常，但不会在内部任务中引发。

main() 协程恢复并响应 CancelledError 异常，报告一条消息。然后它会睡一会儿。任务恢复、完成并返回一个值。

最后，main() 协程恢复，并报告被屏蔽的未来和内部任务的状态。我们可以看到屏蔽的未来被标记为已取消，而内部任务被标记为正常完成并提供返回值。

此示例突出显示了如何使用防护罩来成功保护内部任务不被取消。

```python
cancelled: True
shielded was cancelled
shielded: <Future cancelled>
task: <Task finished name='Task-2' coro=<simple_task() done, defined at ...> result=1>
```

