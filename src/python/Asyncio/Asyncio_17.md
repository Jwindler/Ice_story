# Python 异步: 异步上下文管理器（17）

上下文管理器是一种 Python 构造，它提供了一个类似 try-finally 的环境，具有一致的接口和方便的语法，例如通过“with”表达。

它通常与资源一起使用，确保在我们完成资源后始终关闭或释放资源，无论资源的使用是成功还是因异常而失败。

Asyncio 允许我们开发异步上下文管理器。

我们可以通过定义一个将 __aenter__() 和 __aexit__() 方法实现为协程的对象来在 asyncio 程序中创建和使用异步上下文管理器。



## 1. 什么是异步上下文管理器

异步上下文管理器是一个实现了 __aenter__() 和 __aexit__() 方法的 Python 对象。

在我们深入了解异步上下文管理器的细节之前，让我们回顾一下经典的上下文管理器。



### 1.1. Context Manager

上下文管理器是一个 Python 对象，它实现了 __enter__() 和 __exit__() 方法。

- __enter__() 方法定义了块开头发生的事情，例如打开或准备资源，如文件、套接字或线程池。
- __exit__() 方法定义退出块时发生的情况，例如关闭准备好的资源。

通过“with”表达式使用上下文管理器。通常，上下文管理器对象是在“with”表达式的开头创建的，并且会自动调用 __enter__() 方法。内容的主体通过命名的上下文管理器对象使用资源，然后 __aexit__() 方法在块退出时自动调用，通常或通过异常。

```python
...
# open a context manager
with ContextManager() as manager:
	# ...
# closed automatically
```

这反映了 try-finally 表达式。

```python
...
# create the object
manager = ContextManager()
try:
	manager.__enter__()
	# ...
finally:
	manager.__exit__()
```



### 1.2. Asynchronous Context Manager

“PEP 492 – Coroutines with async and await syntax”引入了异步上下文管理器。

它们提供了一个上下文管理器，可以在进入和退出时挂起。

__aenter__ 和 __aexit__ 方法被定义为协同程序，由调用者等待。这是使用“async with”表达式实现的。

因此，异步上下文管理器只能在 asyncio 程序中使用，例如在调用协程中。

- 什么是“async with”

“async with”表达式用于创建和使用异步上下文管理器。它是“with”表达式的扩展，用于异步程序中的协程。

“async with”表达式就像用于上下文管理器的“with”表达式，除了它允许在协同程序中使用异步上下文管理器。

为了更好地理解“async with”，让我们仔细看看异步上下文管理器。async with 表达式允许协程创建和使用上下文管理器的异步版本。

```python
...
# create and use an asynchronous context manager
async with AsyncContextManager() as manager:
	# ...
```

这相当于：

```python
...
# create or enter the async context manager
manager = await AsyncContextManager()
try:
	# ...
finally:
	# close or exit the context manager
	await manager.close()
```

请注意，我们正在实现与传统上下文管理器大致相同的模式，只是创建和关闭上下文管理器涉及等待协程。

这会暂停当前协程的执行，调度一个新的协程并等待它完成。因此，异步上下文管理器必须实现必须通过 async def 表达式定义的 __aenter__() 和 __aexit__() 方法。这使得它们自己协程也可能等待。



## 2. 如何使用异步上下文管理器

在本节中，我们将探讨如何在我们的 asyncio 程序中定义、创建和使用异步上下文管理器。



### 2.1. 定义

我们可以将异步上下文管理器定义为实现 __aenter__() 和 __aexit__() 方法的 Python 对象。

重要的是，这两种方法都必须使用“async def”定义为协程，因此必须返回可等待对象。

```python
# define an asynchronous context manager
class AsyncContextManager:
    # enter the async context manager
    async def __aenter__(self):
        # report a message
        print('>entering the context manager')
 
    # exit the async context manager
    async def __aexit__(self, exc_type, exc, tb):
        # report a message
        print('>exiting the context manager')
```

因为每个方法都是协程，所以它们本身可能等待协程或任务。

```python
# define an asynchronous context manager
class AsyncContextManager:
    # enter the async context manager
    async def __aenter__(self):
        # report a message
        print('>entering the context manager')
        # block for a moment
        await asyncio.sleep(0.5)
 
    # exit the async context manager
    async def __aexit__(self, exc_type, exc, tb):
        # report a message
        print('>exiting the context manager')
        # block for a moment
        await asyncio.sleep(0.5)
```



### 2.2. 使用

通过“async with”表达式使用异步上下文管理器。这将自动等待进入和退出协程，根据需要暂停调用协程。

```python
...
# use an asynchronous context manager
async with AsyncContextManager() as manager:
	# ...
```

因此，“async with”表达式和异步上下文管理器更普遍地只能在 asyncio 程序中使用，例如在协程中。

现在我们知道如何使用异步上下文管理器，让我们看一个有效的例子。



## 3. 异步上下文管理器和“异步”示例

我们可以探索如何通过“async with”表达式使用异步上下文管理器。

在这个例子中，我们将更新上面的例子，以正常方式使用上下文管理器。

我们将使用“async with”表达式，并在一行中创建并进入上下文管理器。这将自动等待 enter 方法。

然后我们可以在内部块中使用管理器。在这种情况下，我们将只报告一条消息。

退出内部块将自动等待上下文管理器的退出方法。将这个例子与前面的例子进行对比，可以看出“async with”表达式在 asyncio 程序中为我们做了多少繁重的工作。

```python
# SuperFastPython.com
# example of an asynchronous context manager via async with
import asyncio
 
# define an asynchronous context manager
class AsyncContextManager:
    # enter the async context manager
    async def __aenter__(self):
        # report a message
        print('>entering the context manager')
        # block for a moment
        await asyncio.sleep(0.5)
 
    # exit the async context manager
    async def __aexit__(self, exc_type, exc, tb):
        # report a message
        print('>exiting the context manager')
        # block for a moment
        await asyncio.sleep(0.5)
 
# define a simple coroutine
async def custom_coroutine():
    # create and use the asynchronous context manager
    async with AsyncContextManager() as manager:
        # report the result
        print(f'within the manager')
 
# start the asyncio program
asyncio.run(custom_coroutine())
```

运行示例首先创建 main() 协程并将其用作 asyncio 程序的入口点。

main() 协程运行并在“async with”表达式中创建我们的 AsyncContextManager 类的实例。

该表达式自动调用 enter 方法并等待协程。报告一条消息，协程暂时阻塞。

main() 协程恢复并执行上下文管理器的主体，打印一条消息。

块退出，自动等待上下文管理器的退出方法，报告消息并休眠片刻。

这突出了 asyncio 程序中异步上下文管理器的正常使用模式。

```python
>entering the context manager
within the manager
>exiting the context manager
```

