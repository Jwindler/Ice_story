# Python 异步: 同时运行多个协程（10）

asyncio 的一个好处是我们可以同时运行许多协程。这些协同程序可以在一个组中创建并存储，然后同时一起执行。这可以使用 asyncio.gather() 函数来实现。

让我们仔细看看。



## 1. 什么是 Asyncio gather()

asyncio.gather() 模块函数允许调用者将多个可等待对象组合在一起。分组后，可等待对象可以并发执行、等待和取消。

它是一个有用的实用函数，可用于分组和执行多个协程或多个任务。

```python
...
# run a collection of awaitables
results = await asyncio.gather(coro1(), asyncio.create_task(coro2()))
```

在我们可能预先创建许多任务或协程然后希望一次执行它们并等待它们全部完成后再继续的情况下，我们可以使用 asyncio.gather() 函数。

这是一种可能的情况，其中需要许多类似任务的结果，例如具有不同数据的相同任务或协程。

可等待对象可以并发执行，返回结果，并且主程序可以通过使用它所依赖的结果来恢复。

gather() 函数比简单地等待任务完成更强大。它允许将一组可等待对象视为单个可等待对象。

- 通过 await 表达式执行并等待组中的所有可等待对象完成。
- 从所有分组的等待对象中获取结果，稍后通过 result() 方法检索。
- 要通过 cancel() 方法取消的一组等待对象。
- 通过 done() 方法检查组中的所有可等待对象是否已完成。
- 仅当组中的所有任务完成时才执行回调函数。



## 2. 如何使用 Asyncio gather()

在本节中，我们将仔细研究如何使用 asyncio.gather() 函数。

asyncio.gather() 函数将一个或多个可等待对象作为参数。回想一下，可等待对象可能是协程、Future 或 Task。

因此，我们可以调用 gather() 函数：

- 多项任务
- 多个协程
- 任务和协程的混合

```python
...
# execute multiple coroutines
asyncio.gather(coro1(), coro2())
```

如果 Task 对象被提供给 gather()，它们将已经在运行，因为 Tasks 被安排为创建的一部分。asyncio.gather() 函数将可等待对象作为位置参数。

我们不能创建可等待对象的列表或集合并将其提供给收集，因为这会导致错误。

```python
...
# cannot provide a list of awaitables directly
asyncio.gather([coro1(), coro2()])
```

如果首先使用星号运算符 (*) 将其解压缩到单独的表达式中，则可以提供等待列表。

```python
...
# gather with an unpacked list of awaitables
asyncio.gather(*[coro1(), coro2()])
```

如果协程提供给 gather()，它们会自动包装在 Task 对象中。gather() 函数不会阻塞。

相反，它返回一个代表可等待对象组的 asyncio.Future 对象。

```python
...
# get a future that represents multiple awaitables
group = asyncio.gather(coro1(), coro2())
```

一旦创建了 Future 对象，它就会在事件循环中自动调度。awaitable 代表组，组中的所有 awaitable 都会尽快执行。这意味着如果调用者什么都不做，那么预定的可等待对象组将运行（假设调用者挂起）。

这也意味着您不必等待从 gather() 返回的 Future。

```python
...
# get a future that represents multiple awaitables
group = asyncio.gather(coro1(), coro2())
# suspend and wait a while, the group may be executing..
await asyncio.sleep(10)
```

可以等待返回的 Future 对象，它将等待组中的所有可等待对象完成。

```python
...
# run the group of awaitables
await group
```

等待从 gather() 返回的 Future 将返回可等待对象的返回值列表。

如果可等待对象没有返回值，则此列表将包含默认的“无”返回值。

```python
...
# run the group of awaitables and get return values
results = await group
```

这通常在一行中执行。

```python
...
# run tasks and get results on one line
results = await asyncio.gather(coro1(), coro2())
```



## 3. 列表中多个协程的 gather() 示例

预先创建多个协程然后再收集它们是很常见的。这允许程序准备要并发执行的任务，然后立即触发它们的并发执行并等待它们完成。

我们可以手动或使用列表理解将许多协程收集到一个列表中。

```python
...
# create many coroutines
coros = [task_coro(i) for i in range(10)]
```

然后我们可以用列表中的所有协程调用 gather()。协程列表不能直接提供给 gather() 函数，因为这会导致错误。相反，gather() 函数要求将每个可等待对象作为单独的位置参数提供。

这可以通过将列表展开为单独的表达式并将它们传递给 gather() 函数来实现。星号运算符 (*) 将为我们执行此操作。

```python
...
# run the tasks
await asyncio.gather(*coros)
```

将它们结合在一起，下面列出了使用 gather() 运行预先准备好的协程列表的完整示例。

```python
# SuperFastPython.com
# example of gather for many coroutines in a list
import asyncio
 
# coroutine used for a task
async def task_coro(value):
    # report a message
    print(f'>task {value} executing')
    # sleep for a moment
    await asyncio.sleep(1)
 
# coroutine used for the entry point
async def main():
    # report a message
    print('main starting')
    # create many coroutines
    coros = [task_coro(i) for i in range(10)]
    # run the tasks
    await asyncio.gather(*coros)
    # report a message
    print('main done')
 
# start the asyncio program
asyncio.run(main())
```

运行该示例会执行 main() 协程作为程序的入口点。main() 协程然后使用列表理解创建一个包含 10 个协程对象的列表。然后将此列表提供给 gather() 函数，并使用星号运算符将其解压缩为 10 个单独的表达式。

然后 main() 协程等待从调用 gather() 返回的 Future 对象，暂停并等待所有调度的协程完成它们的执行。协程会尽快运行，报告它们独特的消息并在终止前休眠。

只有在组中的所有协程都完成后，main() 协程才会恢复并报告其最终消息。这突出了我们如何准备协程集合并将它们作为单独的表达式提供给 gather() 函数。

```python
main starting
>task 0 executing
>task 1 executing
>task 2 executing
>task 3 executing
>task 4 executing
>task 5 executing
>task 6 executing
>task 7 executing
>task 8 executing
>task 9 executing
main done
```

