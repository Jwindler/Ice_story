# Python 异步: 在 Asyncio 中运行阻塞任务（14）

阻塞任务是阻止当前线程继续进行的任务。

如果在 asyncio 程序中执行阻塞任务，它会停止整个事件循环，从而阻止任何其他协程继续进行。

我们可以通过 asyncio.to_thread() 和 loop.run_in_executor() 函数在 asyncio 程序中异步运行阻塞调用。



## 1. 阻塞任务

asyncio的重点是异步编程和非阻塞IO。然而，我们经常需要在 asyncio 应用程序中执行阻塞函数调用。

这可能有很多原因，例如：

- 执行 CPU 密集型任务，例如计算某事。
- 执行阻塞 IO 绑定任务，如从文件读取或写入。
- 调用不支持 asyncio 的第三方库。

直接在 asyncio 程序中进行阻塞调用将导致事件循环在执行阻塞调用时停止。它不允许其他协程在后台运行。

我们如何在 asyncio 程序中异步执行阻塞调用？



## 2. 如何运行阻塞任务

asyncio 模块提供了两种在 asyncio 程序中执行阻塞调用的方法。

第一种是使用 asyncio.to_thread() 函数。这是在高级 API 中，供应用程序开发人员使用。

asyncio.to_thread() 函数采用要执行的函数名和任何参数。

该函数在单独的线程中执行。它返回一个可以作为独立任务等待或安排的协程。

```python
...
# execute a function in a separate thread
await asyncio.to_thread(task)
```

在返回的协程有机会在事件循环中运行之前，任务不会开始执行。asyncio.to_thread() 函数在后台创建一个 ThreadPoolExecutor 来执行阻塞调用。因此，asyncio.to_thread() 函数仅适用于 IO 绑定任务。

另一种方法是使用 loop.run_in_executor() 函数。

这是在低级异步 API 中，首先需要访问事件循环，例如通过 asyncio.get_running_loop() 函数。

loop.run_in_executor() 函数接受一个执行器和一个要执行的函数。

如果没有为执行器提供，则使用默认执行器，即 ThreadPoolExecutor。

loop.run_in_executor() 函数返回一个可等待对象，如果需要可以等待它。任务将立即开始执行，因此返回的可等待对象不需要等待或安排阻塞调用开始执行。

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

现在我们知道如何在 asyncio 程序中执行阻塞调用，让我们看一些有效的例子。



## 3. 实例

我们可以探索如何使用 asyncio.to_thread() 在 asyncio 程序中执行阻塞 IO 绑定调用。

在这个例子中，我们将定义一个函数来阻塞调用者几秒钟。然后，我们将使用 asyncio.to_thread() 函数在 asyncio 的线程池中异步执行此函数。

这将使呼叫者腾出时间继续其他活动。

```python
# SuperFastPython.com
# example of running a blocking io-bound task in asyncio
import asyncio
import time
 
# a blocking io-bound task
def blocking_task():
    # report a message
    print('Task starting')
    # block for a while
    time.sleep(2)
    # report a message
    print('Task done')
 
# main coroutine
async def main():
    # report a message
    print('Main running the blocking task')
    # create a coroutine for  the blocking task
    coro = asyncio.to_thread(blocking_task)
    # schedule the task
    task = asyncio.create_task(coro)
    # report a message
    print('Main doing other things')
    # allow the scheduled task to start
    await asyncio.sleep(0)
    # await the task
    await task
 
# run the asyncio program
asyncio.run(main())
```

运行示例首先创建 main() 协程并将其作为 asyncio 程序的入口点运行。main() 协程运行并报告一条消息。然后它发出对线程池的阻塞函数调用的调用。然后将协程包装在任务中并独立执行。

main() 协程可以自由地继续其他活动。在这种情况下，它会休眠片刻以允许计划任务开始执行。这使得目标函数可以在后台下发给 ThreadPoolExecutor 并开始运行。

然后 main() 协程挂起并等待任务完成。阻塞函数报告一条消息，休眠 2 秒，然后报告最后一条消息。

这突出了我们如何在一个单独的线程中与 asyncio 程序异步执行阻塞 IO 绑定任务。

```python
Main running the blocking task
Main doing other things
Task starting
Task done
```

