# Python 异步: 使用和查询任务（8）

任务是异步程序的货币。在本节中，我们将仔细研究如何在我们的程序中与它们交互。



## 1. 任务生命周期

异步任务具有生命周期。首先，任务是从协程创建的。然后安排在事件循环中独立执行。在某个时候，它会运行。

在运行时它可能会被挂起，例如等待另一个协程或任务。它可能正常完成并返回结果或因异常而失败。

另一个协程可能会介入并取消任务。最终，它将完成并且无法再次执行。

我们可以将这个生命周期总结如下：

1. 创建
2. 预定
   1. 取消
3. 运行
   1. 暂停
   2. 结果
   3. Exception
   4. 取消
4. 完成

请注意，Suspended、Result、Exception 和 Canceled 本身并不是状态，它们是正在运行的任务的重要转换点。

下图总结了此生命周期，显示了每个阶段之间的转换。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230128170843424.png)



现在我们已经从高层次上熟悉了任务的生命周期，让我们仔细看看每个阶段。



## 2. 如何检查任务状态

创建任务后，我们可以检查任务的状态。我们可能要检查两种状态，它们是：

- 任务是否完成
- 任务是否取消

让我们依次仔细看看每一个。

### 2.1. 检查任务是否完成

我们可以通过 done() 方法检查任务是否完成。如果任务完成，该方法返回 True，否则返回 False。

```python
# check if a task is done
if task.done():
	# ...
```

如果任务有机会运行但现在不再运行，则该任务已完成。已安排的任务未完成。同样，正在运行的任务未完成。

如果出现以下情况，则完成任务：

1. 协程正常结束。
2. 协程显式返回。
3. 协程中出现意外错误或异常
4. 任务被取消。



### 2.2. 检查任务是否取消

我们可以通过 cancelled() 方法检查任务是否被取消。如果任务被取消，该方法返回 True，否则返回 False。

```python
...
# check if a task was canceled
if task.cancelled():
	# ...
```

如果在任务上调用 cancel() 方法并成功完成，则任务被取消，例如 cancel() 返回 True。

如果未调用 cancel() 方法，或者调用了 cancel() 方法但未能取消任务，则不会取消任务。



## 3. 如何获取任务结果

我们可以通过 result() 方法获取任务的结果。这将返回由 Task 包装的协程的返回值，如果包装的协程没有显式返回值，则返回 None 。

```python
...
# get the return value from the wrapped coroutine
value = task.result()
```

如果协程引发未处理的错误或异常，则在调用 result() 方法时会重新引发，并且可能需要处理。

```python
...
try:
	# get the return value from the wrapped coroutine
	value = task.result()
except Exception:
	# task failed and there is no result
```

如果任务被取消，则在调用 result() 方法时会引发 CancelledError 异常，可能需要进行处理。

```python
...
try:
	# get the return value from the wrapped coroutine
	value = task.result()
except asyncio.CancelledError:
	# task was canceled
```

因此，最好先检查任务是否已取消。

```python
...
# check if the task was not canceled
if not task.cancelled():
	# get the return value from the wrapped coroutine
	value = task.result()
else:
	# task was canceled
```

如果任务尚未完成，则在调用 result() 方法时会引发 InvalidStateError 异常，可能需要进行处理。

```python
...
try:
	# get the return value from the wrapped coroutine
	value = task.result()
except asyncio.InvalidStateError:
	# task is not yet done
```

因此，最好先检查任务是否已完成。

```python
...
# check if the task is not done
if not task.done():
	await task
# get the return value from the wrapped coroutine
value = task.result()
```



## 4. 如何获取任务异常

任务包装的协程可能会引发未处理的异常。这实际上会取消任务。

我们可以通过 exception() 方法在任务包装的协程中检索未处理的异常。

```python
...
# get the exception raised by a task
exception = task.exception()
```

如果包装协程中未引发未处理的异常，则返回 None 值。

如果任务被取消，则在调用 exception() 方法时会引发 CancelledError 异常，可能需要对其进行处理。

```python
...
try:
	# get the exception raised by a task
	exception = task.exception()
except asyncio.CancelledError:
	# task was canceled
```

因此，最好先检查任务是否已取消。

```python
...
# check if the task was not canceled
if not task.cancelled():
	# get the exception raised by a task
	exception = task.exception()
else:
	# task was canceled
```

如果任务尚未完成，则在调用 exception() 方法时会引发 InvalidStateError 异常，可能需要进行处理。

```python
...
try:
	# get the exception raised by a task
	exception = task.exception()
except asyncio.InvalidStateError:
	# task is not yet done
```

因此，最好先检查任务是否已完成。

```python
...
# check if the task is not done
if not task.done():
	await task
# get the exception raised by a task
exception = task.exception()
```



## 5. 如何取消任务

我们可以通过 cancel() 方法取消计划任务。如果任务被取消，则 cancel 方法返回 True，否则返回 False。

```python
...
# cancel the task
was_cancelled = task.cancel()
```

如果任务已经完成，则无法取消，cancel() 方法将返回 False，任务不会处于已取消状态。

下次任务有机会运行时，它将引发 CancelledError 异常。如果 CancelledError 异常未在包装协程内处理，任务将被取消。否则，如果在包装协程内处理了 CancelledError 异常，任务将不会被取消。

cancel() 方法还可以接受一个消息参数，该参数将在 CancelledError 的内容中使用。



## 6. 如何在任务中使用回调

我们可以通过 add_done_callback() 方法向任务添加完成回调函数。此方法采用任务完成时要调用的函数的名称。回调函数必须将 Task 实例作为参数。

```python
# done callback function
def handle(task):
	print(task)
 
...
# register a done callback function
task.add_done_callback(handle)
```

回想一下，当包装的协程返回时正常完成、引发未处理的异常或取消任务时，任务可能会完成。add_done_callback() 方法可用于添加或注册任意数量的 done 回调函数。

我们还可以通过 remove_done_callback() 函数删除或注销回调函数。

```python
...
# remove a done callback function
task.remove_done_callback(handle)
```



## 7. 如何设置任务名称

一个任务可能有一个名字。如果多个任务是从同一个协程创建的，那么这个名称会很有用，我们需要一些方法以编程方式区分它们。当通过“名称”参数从协程创建任务时，可以设置名称。

```python
...
# create a task from a coroutine
task = asyncio.create_task(task_coroutine(), name='MyTask')
```

任务的名称也可以通过 set_name() 方法设置。

```python
...
# set the name of the task
task.set_name('MyTask')
```

我们可以通过 get_name() 方法检索任务的名称。

```python
...
# get the name of a task
name = task.get_name()
```



