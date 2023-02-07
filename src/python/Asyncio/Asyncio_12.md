# Python 异步: 等待有时间限制的协程（12）

我们可以使用 asyncio.wait_for() 函数等待 asyncio 任务或协程超时完成。如果在任务完成之前超时已过，任务将被取消。



## 1. 什么是 Asyncio wait_for()

asyncio.wait_for() 函数允许调用者等待 asyncio 任务或协程超时完成。如果没有指定超时，wait_for() 函数将等待直到任务完成。如果在任务完成之前指定了超时并超时，那么任务将被取消。

这允许调用者既可以设置他们愿意等待任务完成的时间，又可以通过在超时结束时取消任务来强制执行超时。

现在我们知道了 asyncio.wait_for() 函数是什么，让我们看看如何使用它。



## 2. 如何使用 Asyncio wait_for()

asyncio.wait_for() 函数接受一个等待和超时。等待对象可能是协程或任务。必须指定超时，并且可以是无超时、整数或浮点秒数。wait_for() 函数返回一个协程，该协程在明确等待或作为任务调度之前不会执行。

```python
...
# wait for a task to complete
await asyncio.wait_for(coro, timeout=10)
```

如果提供协程，则在执行 wait_for() 协程时将其转换为任务。如果在任务完成之前超时已过，任务将被取消，并引发 asyncio.TimeoutError，这可能需要处理。

```python
...
# execute a task with a timeout
try:
	# wait for a task to complete
	await asyncio.wait_for(coro, timeout=1)
except asyncio.TimeoutError:
	# ...
```

如果等待的任务因未处理的异常而失败，则该异常将传播回等待 wait_for() 协程的调用者，在这种情况下可能需要处理它。

```python
...
# execute a task that may fail
try:
	# wait for a task to complete
	await asyncio.wait_for(coro, timeout=1)
except asyncio.TimeoutError:
	# ...
except Exception:
	# ...
```

接下来，让我们看看如何在超时时调用 wait_for()。



## 3. 带有超时的 Asyncio wait_for() 示例

我们可以探索如何在任务完成之前等待具有超时的协程。在此示例中，我们执行上述协程，但调用方等待 0.2 秒或 200 毫秒的固定超时。回想一下，一秒等于 1,000 毫秒。

任务协程被修改，使其休眠一秒以上，确保超时总是在任务完成之前到期。

```python
# SuperFastPython.com
# example of waiting for a coroutine with a timeout
from random import random
import asyncio
 
# coroutine to execute in a new task
async def task_coro(arg):
    # generate a random value between 0 and 1
    value = 1 + random()
    # report message
    print(f'>task got {value}')
    # block for a moment
    await asyncio.sleep(value)
    # report all done
    print('>task done')
 
# main coroutine
async def main():
    # create a task
    task = task_coro(1)
    # execute and wait for the task without a timeout
    try:
        await asyncio.wait_for(task, timeout=0.2)
    except asyncio.TimeoutError:
        print('Gave up waiting, task canceled')
 
# start the asyncio program
asyncio.run(main())
```

运行示例首先创建 main() 协程并将其用作 asyncio 程序的入口点。main() 协程创建任务协程。然后它调用 wait_for() 并传递任务协程并将超时设置为 0.2 秒。

main()协程被挂起，执行task_coro()。它报告一条消息并休眠片刻。main() 协程在超时结束后恢复。 wait_for()协程取消task_coro()协程，main()协程挂起。

task_coro() 再次运行并响应要终止的请求。它引发 TimeoutError 异常并终止。main() 协程恢复并处理由 task_coro() 引发的 TimeoutError。

这突出显示了我们如何调用带超时的 wait_for() 函数，并在任务未在超时内完成时取消任务。

由于使用了随机数，程序每次运行时的输出都会不同。

```python
>task got 0.685375224799321
Gave up waiting, task canceled
```



