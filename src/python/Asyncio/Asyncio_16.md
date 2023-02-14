# Python 异步: 异步生成器（16）

生成器是 Python 的基本组成部分。生成器是一个至少有一个“yield”表达式的函数。它们是可以暂停和恢复的函数，就像协程一样。

实际上，Python 协程是 Python 生成器的扩展。Asyncio 允许我们开发异步生成器。我们可以通过定义一个使用“yield”表达式的协程来创建一个异步生成器。



## 1. 什么是异步生成器

异步生成器是使用 yield 表达式的协程。在我们深入了解异步生成器的细节之前，让我们先回顾一下经典的 Python 生成器。

### 1.1. Generators

生成器是一个 Python 函数，它通过 yield 表达式返回一个值。

```python
# define a generator
def generator():
	for i in range(10):
		yield i
```

生成器执行到 yield 表达式，之后返回一个值。这会在该点暂停生成器。下一次执行生成器时，它将从恢复点恢复并运行直到下一个 yield 表达式。

从技术上讲，生成器函数创建并返回一个生成器迭代器。生成器迭代器执行生成器函数的内容，根据需要产生和恢复。

可以使用内置函数 next() 分步执行生成器。

```python
...
# create the generator
gen = generator()
# step the generator
result = next(gen)
```

虽然，更常见的是迭代生成器直到完成，例如使用 for 循环或列表理解。

```python
...
# traverse the generator and collect results
results = [item for item in generator()]
```

接下来，让我们仔细看看异步生成器。



### 1.2. Asynchronous Generators

异步生成器是使用 yield 表达式的协程。与函数生成器不同，协程可以调度和等待其他协程和任务。

与经典生成器一样，异步生成器函数可用于创建可使用内置的 anext() 函数而不是 next() 函数遍历的异步生成器迭代器。

这意味着异步生成器迭代器实现了 __anext__() 方法并且可以与 async for 表达式一起使用。

这意味着生成器的每次迭代都被安排并执行为可等待的。 “async for”表达式将调度并执行生成器的每次迭代，暂停调用协程并等待结果。



## 2. 如何使用异步生成器

在本节中，我们将仔细研究如何在 asyncio 程序中定义、创建、步进和遍历异步生成器。

让我们从如何定义异步生成器开始。



### 2.1. 定义

我们可以通过定义一个至少有一个 yield 表达式的协程来定义一个异步生成器。

这意味着该函数是使用“async def”表达式定义的。

```python
# define an asynchronous generator
async def async_generator():
	for i in range(10)
		yield i
```

因为异步生成器是一个协程，并且每个迭代器返回一个在 asyncio 事件循环中调度和执行的等待对象，所以我们可以在生成器主体内执行和等待等待对象。

```python
# define an asynchronous generator that awaits
async def async_generator():
	for i in range(10)
		# suspend and sleep a moment
		await asyncio.sleep(1)
		# yield a value to the caller
		yield i
```

接下来，让我们看看如何使用异步生成器。



### 2.2. 创建

要使用异步生成器，我们必须创建生成器。这看起来像是调用它，而是创建并返回一个迭代器对象。

```python
...
# create the iterator
it = async_generator()
```

这将返回一种称为异步生成器迭代器的异步迭代器。



### 2.3. 一步

可以使用 anext() 内置函数遍历生成器的一个步骤，就像使用 next() 函数的经典生成器一样。

结果是等待的可等待对象。

```python
...
# get an awaitable for one step of the generator
awaitable = anext(gen)
# execute the one step of the generator and get the result
result = await awaitable
```

这可以一步实现。

```python
...
# step the async generator
result = await anext(gen)
```



### 2.4. 遍历

还可以使用“async for”表达式在循环中遍历异步生成器，该表达式将自动等待循环的每次迭代。

```python
...
# traverse an asynchronous generator
async for result in async_generator():
	print(result)
```

我们还可以使用带有“async for”表达式的异步列表理解来收集生成器的结果。

```python
...
# async list comprehension with async generator
results = [item async for item in async_generator()]
```



## 3. 异步生成器示例

我们可以探索如何使用“async for”表达式遍历异步生成器。

在此示例中，我们将更新之前的示例以使用“async for”循环遍历生成器直至完成。

此循环将自动等待从生成器返回的每个可等待对象，检索产生的值，并使其在循环体内可用，以便在这种情况下可以报告它。

这可能是异步生成器最常见的使用模式。

```python
# SuperFastPython.com
# example of asynchronous generator with async for loop
import asyncio
 
# define an asynchronous generator
async def async_generator():
    # normal loop
    for i in range(10):
        # block to simulate doing work
        await asyncio.sleep(1)
        # yield the result
        yield i
 
# main coroutine
async def main():
    # loop over async generator with async for loop
    async for item in async_generator():
        print(item)
 
# execute the asyncio program
asyncio.run(main())
```

运行示例首先创建 main() 协程并将其用作 asyncio 程序的入口点。main() 协程运行并启动 for 循环。

异步生成器的一个实例被创建，循环使用 anext() 函数自动单步执行它以返回一个可等待对象。然后循环等待可等待对象并检索一个值，该值可用于报告它的循环体。

然后重复此过程，挂起 main() 协程，执行生成器的迭代，然后挂起和恢复 main() 协程，直到生成器耗尽。

这突出显示了如何使用 async for 表达式遍历异步生成器。

```python
0
1
2
3
4
5
6
7
8
9
```

