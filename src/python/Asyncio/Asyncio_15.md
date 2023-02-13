# Python 异步: 异步迭代器（15）

迭代是 Python 中的基本操作。我们可以迭代列表、字符串和所有其他结构。

Asyncio 允许我们开发异步迭代器。我们可以通过定义一个实现 __aiter__() 和 __anext__() 方法的对象来在 asyncio 程序中创建和使用异步迭代器。



## 1. 什么是异步迭代器

异步迭代器是一个实现了 __aiter__() 和 __anext__() 方法的对象。在我们仔细研究异步迭代器之前，让我们回顾一下经典迭代器。

### 1.1. Iterators

迭代器是实现特定接口的 Python 对象。具体来说，返回迭代器实例的 __iter__() 方法和使迭代器步进一个循环并返回值的 __next__() 方法。可以使用内置函数 next() 步进迭代器或使用 for 循环遍历迭代器。许多 Python 对象是可迭代的，最值得注意的是列表等容器。



### 1.2. Asynchronous Iterators

异步迭代器是实现特定接口的 Python 对象。异步迭代器必须实现 __aiter__() 和 __anext__() 方法。

- __aiter__() 方法必须返回迭代器的一个实例。
- __anext__() 方法必须返回一个步进迭代器的可等待对象。

异步迭代器只能在 asyncio 程序中步进或遍历，例如在协程中。

可以使用 anext() 内置函数步进异步迭代器，该函数返回执行迭代器一步的可等待对象，例如一次调用 __anext__() 方法。

可以使用“async for”表达式遍历异步迭代器，该表达式将在每次迭代时自动调用 anext() 并等待返回的 awaitable 以检索返回值。



## 2. 什么是“async for”循环？

async for 表达式用于遍历异步迭代器。它是一个异步的 for 循环语句。异步迭代器是产生可等待对象的迭代器。您可能还记得 awaitable 是可以等待的对象，例如协程或任务。

异步生成器将自动实现异步迭代器方法，允许它像异步迭代器一样被迭代。await for 表达式允许调用者遍历 awaitable 的异步迭代器并从每个迭代器中检索结果。

这与遍历集合或等待对象列表（例如协程对象）不同，相反，必须使用预期的异步迭代器方法提供返回的等待对象。在内部，async for 循环将根据需要自动解析或等待每个可等待的调度协程。

因为它是一个 for 循环，所以它假定（尽管不要求）每个被遍历的等待对象都会产生一个返回值。async for 循环必须在协程内使用，因为它在内部会使用只能在协程内使用的 await 表达式。async for 表达式可用于在协程中遍历异步迭代器。

```python
...
# traverse an asynchronous iterator
async for item in async_iterator:
	print(item)
```

这不会并行执行 for 循环。 asyncio 无法在一个 Python 线程中一次执行多个协程。

相反，这是一个异步 for 循环。不同的是，执行 for 循环的协程会暂停并在内部等待每个 awaitable。在幕后，这可能需要安排和等待协程，或者等待任务。我们也可以在列表理解中使用 async for 表达式。

```python
...
# build a list of results
results = [item async for item async_iterator]
```

这将构建异步迭代器的返回值列表。



## 3. 如何使用异步迭代器

在本节中，我们将仔细研究如何在 asyncio 程序中定义、创建、步进和遍历异步迭代器。让我们从如何定义异步迭代器开始。

- 定义异步迭代器

我们可以通过定义一个实现了 __aiter__() 和 __anext__() 方法的类来定义一个异步迭代器。这些方法通常在 Python 对象上定义。重要的是，因为 __anext__() 函数必须返回一个可等待对象，所以它必须使用“async def”表达式定义。迭代完成后，__anext__() 方法必须引发 StopAsyncIteration 异常。

```python
# define an asynchronous iterator
class AsyncIterator():
    # constructor, define some state
    def __init__(self):
        self.counter = 0
 
    # create an instance of the iterator
    def __aiter__(self):
        return self
 
    # return the next awaitable
    async def __anext__(self):
        # check for no further items
        if self.counter >= 10:
            raise StopAsyncIteration
        # increment the counter
        self.counter += 1
        # return the counter value
        return self.counter
```

因为异步迭代器是一个协程，并且每个迭代器返回一个在 asyncio 事件循环中调度和执行的等待对象，所以我们可以在迭代器的主体内执行和等待等待对象。

```python
...
# return the next awaitable
async def __anext__(self):
    # check for no further items
    if self.counter >= 10:
        raise StopAsyncIteration
    # increment the counter
    self.counter += 1
    # simulate work
    await asyncio.sleep(1)
    # return the counter value
    return self.counter
```



- 创建异步迭代器

要使用异步迭代器，我们必须创建迭代器。这涉及正常创建 Python 对象。

```python
...
# create the iterator
it = AsyncIterator()
```

这将返回一个“异步迭代器”，它是“异步迭代器”的一个实例。



- 迭代一个异步迭代器

可以使用 anext() 内置函数遍历迭代器的一步，就像使用 next() 函数的经典迭代器一样。结果是等待的可等待对象。

```python
...
# get an awaitable for one step of the iterator
awaitable = anext(it)
# execute the one step of the iterator and get the result
result = await awaitable
```

这可以一步实现。

```python
...
# step the async iterator
result = await anext(it)
```



- 遍历异步迭代器

异步迭代器也可以使用“async for”表达式在循环中遍历，该表达式将自动等待循环的每次迭代。

```python
...
# traverse an asynchronous iterator
async for result in AsyncIterator():
	print(result)
```

我们还可以使用带有“async for”表达式的异步列表理解来收集迭代器的结果。

```python
...
# async list comprehension with async iterator
results = [item async for item in AsyncIterator()]
```



## 4. 异步迭代器示例

我们可以探索如何使用“async for”表达式遍历异步迭代器。在此示例中，我们将更新之前的示例，以使用“async for”循环遍历迭代器直至完成。

此循环将自动等待从迭代器返回的每个可等待对象，检索返回值，并使其在循环体内可用，以便在这种情况下可以报告它。这可能是异步迭代器最常见的使用模式。

```python
# SuperFastPython.com
# example of an asynchronous iterator with async for loop
import asyncio
 
# define an asynchronous iterator
class AsyncIterator():
    # constructor, define some state
    def __init__(self):
        self.counter = 0
 
    # create an instance of the iterator
    def __aiter__(self):
        return self
 
    # return the next awaitable
    async def __anext__(self):
        # check for no further items
        if self.counter >= 10:
            raise StopAsyncIteration
        # increment the counter
        self.counter += 1
        # simulate work
        await asyncio.sleep(1)
        # return the counter value
        return self.counter
 
# main coroutine
async def main():
    # loop over async iterator with async for loop
    async for item in AsyncIterator():
        print(item)
 
# execute the asyncio program
asyncio.run(main())
```

运行示例首先创建 main() 协程并将其用作 asyncio 程序的入口点。main() 协程运行并启动 for 循环。



异步迭代器的一个实例被创建，循环使用 anext() 函数自动单步执行它以返回一个可等待对象。然后循环等待可等待对象并检索一个值，该值可用于报告它的循环体。然后重复这个过程，挂起 main() 协程，执行迭代器和挂起的一个步骤，然后恢复 main() 协程，直到迭代器耗尽。

一旦迭代器的内部计数器达到 10，就会引发 StopAsyncIteration。这不会终止程序。相反，它由“async for”表达式预期和处理并中断循环。

这突出显示了如何使用 async for 表达式遍历异步迭代器。

```python
1
2
3
4
5
6
7
8
9
10
```

