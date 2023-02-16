# Python 异步: 异步推导式（18）

当我们想到“pythonic”时，理解，如列表和字典理解是 Python 的一个特性。

这是我们执行循环的一种方式，与许多其他语言不同。

Asyncio 允许我们使用异步推导式。

我们可以通过“async for”表达式使用异步推导式来遍历异步生成器和异步迭代器。



## 1. 什么是异步推导式

异步推导式是经典推导式的异步版本。Asyncio 支持两种类型的异步推导式，它们是“async for”推导式和“await”推导式。

在我们看每之前，让我们首先回顾一下经典的推导式。



## 2. 推导式

推导式允许以简洁的方式创建列表、字典和集合等数据集合。列表推导式允许从新列表表达式中的 for 表达式创建列表。

```python
...
# create a list using a list comprehension
result = [a*2 for a in range(100)]
```

还支持推导式来创建字典和集合。

```python
...
# create a dict using a comprehension
result = {a:i for a,i in zip(['a','b','c'],range(3))}
# create a set using a comprehension
result = {a for a in [1, 2, 3, 2, 3, 1, 5, 4]}
```



## 3. 异步推导式

异步推导式允许使用带有异步可迭代对象的“async for”表达式来创建列表、集合或字典。

```python
...
# async list comprehension with an async iterator
result = [a async for a in aiterable]
```

这将根据需要创建和安排协程或任务，并将其结果放入列表中。

回想一下，“async for”表达式只能在协程和任务中使用。

另外，回想一下异步迭代器是一个产生可等待对象的迭代器。

“async for”表达式允许调用者遍历等待对象的异步迭代器并从每个对象中检索结果。

在内部，async for 循环将根据需要自动解析或等待每个可等待的调度协程。

异步生成器自动实现异步迭代器的方法，也可用于异步推导式。

```python
...
# async list comprehension with an async generator
result = [a async for a in agenerator]
```





## 4. Await 推导式

“等待”表达式也可以在列表、集合或字典理解中使用，称为等待推导式。

与异步推导式一样，它只能在异步协程或任务中使用。

这允许通过挂起和等待一系列可等待对象来创建数据结构，如列表。

```python
...
# await list compression with a collection of awaitables
results = [await a for a in awaitables]
```

这将通过依次等待每个可等待对象来创建结果列表。

当前协程将被挂起以顺序执行可等待对象，这与使用 asyncio.gather() 并发执行它们不同，而且可能更慢。