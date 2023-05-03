# 在 Python 中将 Tqdm 与 Asyncio 结合使用



## 简介

### 困扰

在 Python 中使用并发编程来提高效率对于数据科学家来说并不罕见。在后台观察各种子进程或并发线程以保持我的计算或 IO 绑定任务的顺序总是令人满意的。

但是还有一点困扰我的是，当我在后台并发处理成百上千个文件或者执行成百上千个进程时，我总是担心会不会有几个任务偷偷挂了，整个代码永远跑不完。我也很难知道代码现在在哪里执行。

最糟糕的是，当我看着一个空白屏幕时，很难说出我的代码需要多长时间才能执行或 ETA 是多少。这对我安排工作日程的能力非常不利。

因此，我想要一种方法让我知道代码执行到了哪里。



### 已有方法

比较传统的做法是任务之间共享一块内存区域，在这块内存区域放一个计数器，当一个任务结束的时候让这个计数器+1，然后用一个线程不停的打印这个计数器的值。

这从来都不是一个好的解决方案：一方面，我需要在你现有的业务逻辑中添加一段用于计数的代码，这违反了“低耦合，高内聚”的原则。另一方面，由于线程安全问题，我必须非常小心锁定机制，这会导致不必要的性能问题。



### tqdm

![](https://s2.loli.net/2023/05/03/jocTbw8DSQm2pHs.png)



有一天，我发现了 tqdm 库，它使用进度条来可视化我的代码进度。我可以使用进度条来可视化我的 asyncio 任务的完成和预计到达时间吗？

那么[本文](https://towardsdatascience.com/using-tqdm-with-asyncio-in-python-5c0f6e747d55 "Source")我把这个方法分享给大家，让每个程序员都有机会监控自己并发任务的进度。



## 异步

在我们开始之前，我希望您了解一些 Python asyncio 的背景知识。我的文章描述了[asyncio](https://github.com/Jwindler/Ice_story "异步")的一些常用API的用法，这将有助于我们更好地理解tqdm的设计：

![](https://s2.loli.net/2023/05/03/TSYzCoKX7RI9Wj5.png)



## tqdm 概述

如官方网站所述，tqdm 是一个显示循环进度条的工具。它使用简单、高度可定制并且占用资源少。

一个典型的用法是将一个可迭代对象传递给 tqdm 构造函数，然后你会得到一个如下所示的进度条：

```python
from time import sleep
from tqdm import tqdm


def main():
    for _ in tqdm(range(100)):
        # do something in the loop
        sleep(0.1)


if __name__ == "__main__":
    main()
```

或者您可以在读取文件时手动浏览并更新进度条的进度：

```python
import os
from tqdm import tqdm


def main():
    filename = "../data/large-dataset"
    with (tqdm(total=os.path.getsize(filename)) as bar,
            open(filename, "r", encoding="utf-8") as f):
        for line in f:
            bar.update(len(line))


if __name__ == "__main__":
    main()
```

![](https://s2.loli.net/2023/05/03/LujznZ9q6wMhdmS.png)





## 将 tqdm 与异步集成

总体而言，tqdm 非常易于使用。但是，GitHub 上需要更多关于将 tqdm 与 asyncio 集成的信息。所以我深入研究了源代码，看看 tqdm 是否支持 asyncio。

幸运的是，最新版本的 tqdm 提供了包 tqdm.asyncio，它提供了类 tqdm_asyncio。

tqdm_asyncio 类有两个相关的方法。一个是 tqdm_asyncio.as_completed。从源码可以看出，它是对asyncio.as_completed的包装：

```python
@classmethod
    def as_completed(cls, fs, *, loop=None, timeout=None, total=None, **tqdm_kwargs):
        """
        Wrapper for `asyncio.as_completed`.
        """
        if total is None:
            total = len(fs)
        kwargs = {}
        if version_info[:2] < (3, 10):
            kwargs['loop'] = loop
        yield from cls(asyncio.as_completed(fs, timeout=timeout, **kwargs),
                       total=total, **tqdm_kwargs)
```

另一个是 tqdm_asyncio.gather ，从源代码可以看出，它基于模拟 asyncio.gather 功能的 tqdm_asyncio.as_completed 的实现：

```python
@classmethod
    async def gather(cls, *fs, loop=None, timeout=None, total=None, **tqdm_kwargs):
        """
        Wrapper for `asyncio.gather`.
        """
        async def wrap_awaitable(i, f):
            return i, await f

        ifs = [wrap_awaitable(i, f) for i, f in enumerate(fs)]
        res = [await f for f in cls.as_completed(ifs, loop=loop, timeout=timeout,
                                                 total=total, **tqdm_kwargs)]
        return [i for _, i in sorted(res)]
```

所以，接下来，我将描述这两个API的用法。在开始之前，我们还需要做一些准备工作。在这里，我写了一个简单的方法来模拟一个随机休眠时间的并发任务：

```python
import asyncio
import random

from tqdm.asyncio import tqdm_asyncio


class AsyncException(Exception):
    def __int__(self, message):
        super.__init__(self, message)


async def some_coro(simu_exception=False):
    delay = round(random.uniform(1.0, 5.0), 2)

    # We will simulate throwing an exception if simu_exception is True
    if delay > 4 and simu_exception:
        raise AsyncException("something wrong!")

    await asyncio.sleep(delay)

    return delay
```

紧接着，我们将创建 2000 个并发任务，然后使用 tqdm_asyncio.gather 而不是熟悉的 asyncio.gather 方法来查看进度条是否正常工作：

```python
async def main():
    tasks = []
    for _ in range(2000):
        tasks.append(some_coro())
    await tqdm_asyncio.gather(*tasks)

    print(f"All tasks done.")


if __name__ == "__main__":
    asyncio.run(main())
```

![](https://s2.loli.net/2023/05/03/SnXFiPkBgDVf3dr.png)



或者让我们用 tqdm_asyncio.as_completed 替换 tqdm_asyncio.gather 并重试：

```python
async def main():
    tasks = []
    for _ in range(2000):
        tasks.append(some_coro())

    for done in tqdm_asyncio.as_completed(tasks):
        await done

    print(f"The tqdm_asyncio.as_completed also works fine.")


if __name__ == "__main__":
    asyncio.run(main())
```

![](https://s2.loli.net/2023/05/03/R4wdkLYH78fjXgp.png)