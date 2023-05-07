# 使用这些方法让你的 Python 并发任务执行得更好



## 问题

一直以来，Python的多线程性能因为GIL而一直没有达到预期。

所以从 3.4 版本开始，Python 引入了 asyncio 包，通过并发的方式并发执行 IO-bound 任务。经过多次迭代，asyncio API 的效果非常好，并发任务的性能相比多线程版本有了很大的提升。

但是，程序员在使用asyncio时还是会犯很多错误：

一个错误如下图所示，直接使用await协程方法，将对并发任务的调用从异步变为同步，最终失去并发特性。

```python
async def main():
    result_1 = await some_coro("name-1")
    result_2 = await some_coro("name-2")
```

另一个错误如下图所示，虽然程序员意识到他需要使用create_task创建一个任务在后台执行。而下面这种一个一个等待任务的方式，将不同时序的任务变成了有序的等待。

```python
async def main():
    task_1 = asyncio.create_task(some_coro("name-1"))
    task_2 = asyncio.create_task(some_coro("name-2"))
    
    result_1 = await task_1
    result_2 = await task_2
```

此代码将等待 task_1 先完成，而不管 task_2 是否先完成。



## 什么是并发任务执行？

那么，什么是真正的并发任务呢？我们用一张图来说明：

![](https://s2.loli.net/2023/05/07/2PeoCzWlNyEZtjV.png)



如图所示，一个并发流程应该由两部分组成：启动后台任务，将后台任务重新加入主函数，并获取结果。

大多数读者已经知道如何使用 create_task 启动后台任务。今天，我将介绍几种等待后台任务完成的方法以及每种方法的最佳实践。



## 开始

在开始介绍今天的主角之前，我们需要准备一个示例async方法来模拟IO绑定的方法调用，以及一个自定义的AsyncException，可以用来在测试抛出异常时友好地提示异常信息：

```python
from random import random, randint
import asyncio


class AsyncException(Exception):
    def __init__(self, message, *args, **kwargs):
        self.message = message
        super(*args, **kwargs)

    def __str__(self):
        return self.message


async def some_coro(name):
    print(f"Coroutine {name} begin to run")
    value = random()

    delay = randint(1, 4)
    await asyncio.sleep(delay)
    if value > 0.5:
        raise AsyncException(f"Something bad happen after delay {delay} second(s)")
    print(f"Coro {name} is Done. with delay {delay} second(s)")
    return value
```



## 并发执行方法比较

### 1. asyncio.gather

asyncio.gather 可用于启动一组后台任务，等待它们完成执行，并获取结果列表：

```python
async def main():
    aws, results = [], []
    for i in range(3):
        aws.append(asyncio.create_task(some_coro(f'name-{i}')))

    results = await asyncio.gather(*aws)  # need to unpack the list
    for result in results:
        print(f">got : {result}")

asyncio.run(main())
```

asyncio.gather 虽然组成了一组后台任务，但不能直接接受一个列表或集合作为参数。如果需要传入包含后台任务的列表，请解包。

asyncio.gather 接受一个 return_exceptions 参数。当return_exception的值为False时，任何后台任务抛出异常，都会抛给gather方法的调用者。而 gather 方法的结果列表是空的。

```python
async def main():
    aws, results = [], []
    for i in range(3):
        aws.append(asyncio.create_task(some_coro(f'name-{i}')))

    try:
        results = await asyncio.gather(*aws, return_exceptions=False)  # need to unpack the list
    except AsyncException as e:
        print(e)
    for result in results:
        print(f">got : {result}")

asyncio.run(main())
```

![](https://s2.loli.net/2023/05/07/KilkqbY69mr2wJE.png)



当return_exception的值为True时，后台任务抛出的异常不会影响其他任务的执行，最终会合并到结果列表中一起返回。

```python
results = await asyncio.gather(*aws, return_exceptions=True)
```

![](https://s2.loli.net/2023/05/07/nTGK41vM763IbmH.png)



接下来我们看看为什么gather方法不能直接接受一个列表，而是要对列表进行解包。因为当一个列表被填满并执行时，我们很难在等待任务完成时向列表中添加新任务。但是 gather 方法可以使用嵌套组将现有任务与新任务混合，解决了中间无法添加新任务的问题：

```python
async def main():
    aws, results = [], []
    for i in range(3):
        aws.append(asyncio.create_task(some_coro(f'name-{i}')))
    group_1 = asyncio.gather(*aws)  # note we don't use await now
    # when some situation happen, we may add a new task
    group_2 = asyncio.gather(group_1, asyncio.create_task(some_coro("a new task")))
    results = await group_2
    for result in results:
        print(f">got : {result}")

asyncio.run(main())
```

但是gather不能直接设置timeout参数。如果需要为所有正在运行的任务设置超时时间，就用这个姿势，不够优雅。

```python
async def main():
    aws, results = [], []
    for i in range(3):
        aws.append(asyncio.create_task(some_coro(f'name-{i}')))

    results = await asyncio.wait_for(asyncio.gather(*aws), timeout=2)
    for result in results:
        print(f">got : {result}")

asyncio.run(main())
```



### 2. asyncio.as_completed

有时，我们必须在完成一个后台任务后立即开始下面的动作。比如我们爬取一些数据，马上调用机器学习模型进行计算，gather方法不能满足我们的需求，但是我们可以使用as_completed方法。

在使用 asyncio.as_completed 方法之前，我们先看一下这个方法的源码。

```python
# This is *not* a @coroutine!  It is just an iterator (yielding Futures).
def as_completed(fs, *, timeout=None):
  # ...
  for f in todo:
      f.add_done_callback(_on_completion)
  if todo and timeout is not None:
      timeout_handle = loop.call_later(timeout, _on_timeout)
  for _ in range(len(todo)):
      yield _wait_for_one()
```

源码显示as_completed不是并发方法，返回一个带有yield语句的迭代器。所以我们可以直接遍历每个完成的后台任务，我们可以对每个任务单独处理异常，而不影响其他任务的执行：

```python
async def main():
    aws = []
    for i in range(5):
        aws.append(asyncio.create_task(some_coro(f"name-{i}")))

    for done in asyncio.as_completed(aws):  # we don't need to unpack the list
        try:
            result = await done
            print(f">got : {result}")
        except AsyncException as e:
            print(e)

asyncio.run(main())
```

as_completed 接受超时参数，超时后当前迭代的任务会抛出asyncio.TimeoutError：

```python
async def main():
    aws = []
    for i in range(5):
        aws.append(asyncio.create_task(some_coro(f"name-{i}")))

    for done in asyncio.as_completed(aws, timeout=2):  # we don't need to unpack the list
        try:
            result = await done
            print(f">got : {result}")
        except AsyncException as e:
            print(e)
        except asyncio.TimeoutError: # we need to handle the TimeoutError
            print("time out.")

asyncio.run(main())
```

![](https://s2.loli.net/2023/05/07/PrdloIN3ksBMmVh.png)



as_complete在处理任务执行的结果方面比gather灵活很多，但是在等待的时候很难往原来的任务列表中添加新的任务。



### 3. asyncio.wait

asyncio.wait 的调用方式与 as_completed 相同，但返回一个包含两个集合的元组：done 和 pending。 done 保存已完成执行的任务，而 pending 保存仍在运行的任务。

asyncio.wait 接受一个 return_when 参数，它可以取三个枚举值：

- 当return_when为asyncio.ALL_COMPLETED时，done存放所有完成的任务，pending为空。
- 当 return_when 为 asyncio.FIRST_COMPLETED 时，done 持有所有已完成的任务，而 pending 持有仍在运行的任务。

```python
async def main():
    aws = set()
    for i in range(5):
        aws.add(asyncio.create_task(some_coro(f"name-{i}")))

    done, pending = await asyncio.wait(aws, return_when=asyncio.FIRST_COMPLETED)
    for task in done:
        try:
            result = await task
            print(f">got : {result}")
        except AsyncException as e:
            print(e)
    print(f"the length of pending is {len(pending)}")

asyncio.run(main())
```

![](https://s2.loli.net/2023/05/07/CTQ4d9EGNwcf6qm.png)



- 当return_when为asyncio.FIRST_EXCEPTION时，done存放抛出异常并执行完毕的任务，pending存放仍在运行的任务。



当 return_when 为 asyncio.FIRST_COMPLETED 或 asyncio.FIRST_EXECEPTION 时，我们可以递归调用 asyncio.wait，这样我们就可以添加新的任务，并根据情况一直等待所有任务完成。

```python
async def main():
    pending = set()
    for i in range(5):
        pending.add(asyncio.create_task(some_coro(f"name-{i}")))  # note the type and name of the task list

    while pending:
        done, pending = await asyncio.wait(pending, return_when=asyncio.FIRST_EXCEPTION)
        for task in done:
            try:
                result = await task
                print(f">got : {result}")
            except AsyncException as e:
                print(e)
                pending.add(asyncio.create_task(some_coro("a new task")))
    print(f"the length of pending is {len(pending)}")

asyncio.run(main())
```

![](https://s2.loli.net/2023/05/07/avRzVqucl6Bd5GQ.png)



### 4. asyncio.TaskGroup

在 Python 3.11 中，asyncio 引入了新的 TaskGroup API，正式让 Python 支持结构化并发。此功能允许您以更 Pythonic 的方式管理并发任务的生命周期。



## 总结

[本文](https://towardsdatascience.com/use-these-methods-to-make-your-python-concurrent-tasks-perform-better-b693b7a633e1 "Source")介绍了 asyncio.gather、asyncio.as_completed 和 asyncio.wait API，还回顾了 Python 3.11 中引入的新 asyncio.TaskGroup 特性。

根据实际需要使用这些后台任务管理方式可以让我们的asyncio并发编程更加灵活。

