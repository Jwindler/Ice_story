# Python 异步: 定义、创建和运行协程（5）

我们可以在我们的 Python 程序中定义协程，就像定义新的子例程（函数）一样。一旦定义，协程函数可用于创建协程对象。“asyncio”模块提供了在事件循环中运行协程对象的工具，事件循环是协程的运行时。



## 1. 如何定义协程

协程可以通过“async def”表达式定义。这是用于定义子例程的“def”表达式的扩展。它定义了一个可以创建的协程，并返回一个协程对象。

```python
# define a coroutine
async def custom_coro():
	# ...
```

用“async def”表达式定义的协程被称为“协程函数”。

然后协程可以在其中使用特定于协程的表达式，例如 await、async for 和 async with。

```python
# define a coroutine
async def custom_coro():
	# await another coroutine
	await asyncio.sleep(1)
```



## 2. 如何创建协程

一旦定义了协程，就可以创建它。这看起来像是在调用一个子程序。

```python
...
# create a coroutine
coro = custom_coro()
```

这不会执行协程。它返回一个“协程”对象。“协程”Python 对象具有方法，例如 send() 和 close()。它是一种类型。

我们可以通过创建协程实例并调用 type() 内置函数来报告其类型来证明这一点。

```python
# SuperFastPython.com
# check the type of a coroutine
 
# define a coroutine
async def custom_coro():
    # await another coroutine
    await asyncio.sleep(1)
 
# create the coroutine
coro = custom_coro()
# check the type of the coroutine
print(type(coro))
```

运行示例报告创建的协程是一个“协程”类。我们还会得到一个 RuntimeError，因为协程已创建但从未执行过，我们将在下一节中探讨它。

```python
<class 'coroutine'>
sys:1: RuntimeWarning: coroutine 'custom_coro' was never awaited
```

协程对象是可等待的。这意味着它是一个实现了 __await__() 方法的 Python 类型。



## 3. 如何从 Python 运行协程

可以定义和创建协程，但它们只能在事件循环中执行。执行协程的事件循环，管理协程之间的协作多任务处理。

启动协程事件循环的典型方法是通过 asyncio.run() 函数。此函数接受一个协程并返回协程的值。提供的协程可以用作基于协程的程序的入口点。

```python
# SuperFastPython.com
# example of running a coroutine
import asyncio
# define a coroutine
async def custom_coro():
    # await another coroutine
    await asyncio.sleep(1)
 
# main coroutine
async def main():
    # execute my custom coroutine
    await custom_coro()
 
# start the coroutine program
asyncio.run(main())
```

现在我们知道如何定义、创建和运行协程，让我们花点时间了解事件循环。

