# Python 异步: 等待任务集合（11）

我们可以通过 asyncio.wait() 函数等待异步任务完成。可以等待不同的条件，例如所有任务完成、第一个任务完成以及第一个任务因异常而失败。

让我们仔细看看。



## 1. 什么是 asyncio.wait()

asyncio.wait() 函数可用于等待一组异步任务完成。回想一下，asyncio 任务是包装协程的 asyncio.Task 类的一个实例。它允许独立调度和执行协程，Task 实例提供任务句柄以查询状态和获取结果。

wait() 函数允许我们等待一组任务完成。等待调用可以配置为等待不同的条件，例如所有任务完成、第一个任务完成以及第一个任务因错误而失败。

接下来，让我们看看如何使用 wait() 函数。



## 2. 如何使用 asyncio.wait()

asyncio.wait() 函数采用可等待对象的集合，通常是 Task 对象。

这可以是我们创建的列表、字典或任务对象集，例如通过在列表理解中调用 asyncio.create_task() 函数。

```python
...
# create many tasks
tasks = [asyncio.create_task(task_coro(i)) for i in range(10)]
```

asyncio.wait() 在满足任务集合的某些条件之前不会返回。默认情况下，条件是所有任务都已完成。

wait() 函数返回两个集合的元组。第一个集合包含所有满足条件的任务对象，第二个集合包含所有其他尚未满足条件的任务对象。

这些集被称为“完成”集和“待定”集。

```python
...
# wait for all tasks to complete
done, pending = await asyncio.wait(tasks)
```

从技术上讲，asyncio.wait() 是一个返回协程的协程函数。然后我们可以等待这个协程，它将返回集合的元组。

```python
...
# create the wait coroutine
wait_coro = asyncio.wait(tasks)
# await the wait coroutine
tuple = await wait_coro
```

等待的条件可以由默认设置为 asyncio.ALL_COMPLETED 的“return_when”参数指定。

```python
...
# wait for all tasks to complete
done, pending = await asyncio.wait(tasks, return_when=asyncio.ALL_COMPLETED)
```

我们可以通过将 return_when 设置为 FIRST_COMPLETED 来等待第一个任务完成

```python
...
# wait for the first task to be completed
done, pending = await asyncio.wait(tasks, return_when=asyncio.FIRST_COMPLETED)
```

当第一个任务完成并在完成集中返回时，其余任务不会被取消并继续并发执行。

我们可以通过将 return_when 设置为 FIRST_EXCEPTION 来等待第一个任务因异常而失败。

```python
...
# wait for the first task to fail
done, pending = await asyncio.wait(tasks, return_when=asyncio.FIRST_EXCEPTION)
```

在这种情况下，完成集将包含第一个因异常而失败的任务。如果没有任务因异常而失败，则完成集将包含所有任务，只有在所有任务完成后 wait() 才会返回。

我们可以通过以秒为单位的“超时”参数指定我们愿意等待给定条件的时间。

如果在满足条件之前超时到期，则返回任务元组以及当时满足条件的任何任务子集，例如如果等待所有任务完成，则完成的任务子集。

```python
...
# wait for all tasks to complete with a timeout
done, pending = await asyncio.wait(tasks, timeout=3)
```

如果在满足条件之前达到超时，则不会引发异常并且不会取消剩余任务。

现在我们知道如何使用 asyncio.wait() 函数，让我们看一些有效的例子。



## 3. 等待所有任务的示例

我们可以探索如何使用 asyncio.wait() 等待所有任务。在这个例子中，我们将定义一个简单的任务协程，它生成一个随机值，休眠几分之一秒，然后用生成的值报告一条消息。

然后，主协程将与协程一起在列表理解中创建许多任务，然后等待所有任务完成。

```python
# SuperFastPython.com
# example of waiting for all tasks to complete
from random import random
import asyncio
 
# coroutine to execute in a new task
async def task_coro(arg):
    # generate a random value between 0 and 1
    value = random()
    # block for a moment
    await asyncio.sleep(value)
    # report the value
    print(f'>task {arg} done with {value}')
 
# main coroutine
async def main():
    # create many tasks
    tasks = [asyncio.create_task(task_coro(i)) for i in range(10)]
    # wait for all tasks to complete
    done,pending = await asyncio.wait(tasks)
    # report results
    print('All done')
 
# start the asyncio program
asyncio.run(main())
```

运行示例首先创建 main() 协程并将其用作 asyncio 程序的入口点。

然后 main() 协程在列表理解中创建一个包含十个任务的列表，每个任务提供一个从 0 到 9 的唯一整数参数。

然后 main() 协程被挂起并等待所有任务完成。任务执行。每个生成一个随机值，休眠片刻，然后报告其生成的值。

所有任务完成后，main() 协程恢复并报告最终消息。这个例子强调了我们如何使用 wait() 函数来等待一组任务完成。

这可能是该函数最常见的用法。请注意，由于使用了随机数，每次运行程序时结果都会不同。

```python
>task 5 done with 0.0591009105682192
>task 8 done with 0.10453715687017351
>task 0 done with 0.15462838864295925
>task 6 done with 0.4103492027393125
>task 9 done with 0.45567100006991623
>task 2 done with 0.6984682905809402
>task 7 done with 0.7785363531316224
>task 3 done with 0.827386088873161
>task 4 done with 0.9481344994700972
>task 1 done with 0.9577302665040541
All done
```

