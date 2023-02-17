# Python 异步: 在非阻塞子进程中运行命令（19）

我们可以从 asyncio 执行命令。该命令将在我们可以使用非阻塞 I/O 写入和读取的子进程中运行。



## 1. 什么是 asyncio.subprocess.Process

asyncio.subprocess.Process 类提供了由 asyncio 运行的子进程的表示。它在 asyncio 程序中提供子进程的句柄，允许对其执行操作，例如等待和终止它。

该 API 与 multiprocessing.Process 类非常相似，可能与 subprocess.Popen 类更相似。具体来说，它与 subprocess.Popen 共享 wait()、communicate() 和 send_signal() 等方法以及 stdin、stdout 和 stderr 等属性。

现在我们知道了 asyncio.subprocess.Process 类是什么，让我们看看如何在我们的 asyncio 程序中使用它。

我们不直接创建 asyncio.subprocess.Process。相反，在 asyncio 程序中执行子进程时，会为我们创建一个类的实例。

有两种方法可以将外部程序作为子流程执行并获取 Process 实例，它们是：

- asyncio.create_subprocess_exec() 用于直接运行命令。
- asyncio.create_subprocess_shell() 用于通过 shell 运行命令。

让我们依次看一下每个示例。



## 2. 如何直接运行命令

命令是在命令行（终端或命令提示符）上执行的程序。这是另一个直接运行的程序。

Linux 和 macOS 上的常见示例可能是：

- ‘ls’ 列出目录的内容
- ‘cat’报告文件的内容
- “data”报告日期
- ‘echo’ 报告一个字符串
- ‘sleep’ 睡眠几秒钟

我们可以通过 create_subprocess_exec() 函数从 asyncio 程序执行命令。

asyncio.create_subprocess_exec() 函数接受一个命令并直接执行它。

这很有用，因为它允许命令在子进程中执行，并允许 asyncio 协程读取、写入和等待它。

与 asyncio.create_subprocess_shell() 函数不同，asyncio.create_subprocess_exec() 不会使用 shell 执行命令。

这意味着 shell 提供的功能，例如 shell 变量、脚本和通配符，在执行命令时不可用。

这也意味着执行命令可能更安全，因为没有机会进行 shell 注入。

现在我们知道了 asyncio.create_subprocess_exec() 的作用，让我们看看如何使用它。



### 2.1. 如何使用 Asyncio create_subprocess_exec()

asyncio.create_subprocess_exec() 函数将在子进程中执行给定的字符串命令。

它返回一个代表子进程的 asyncio.subprocess.Process 对象。

create_subprocess_exec() 函数是一个协程，这意味着我们必须等待它。它会在子流程启动后返回，而不是在子流程完成时返回。

```python
...
# execute a command in a subprocess
process = await asyncio.create_subprocess_exec('ls')
```

正在执行的命令的参数必须作为后续参数提供给 create_subprocess_exec() 函数。

```python
...
# execute a command with arguments in a subprocess
process = await asyncio.create_subprocess_exec('ls', '-l')
```

我们可以通过等待 wait() 方法来等待子进程完成。

```python
...
# wait for the subprocess to terminate
await process.wait()
```

我们可以通过调用 terminate() 或 kill() 方法直接停止子进程，这将在子进程中引发一个信号。

```python
...
# terminate the subprocess
process.terminate()
```

命令的输入和输出将由 stdin、stderr 和 stdout 处理。我们可以让 asyncio 程序处理子进程的输入或输出。

这可以通过指定输入或输出流并指定要重定向的常量来实现，例如 asyncio.subprocess.PIPE。

例如，我们可以将命令的输出重定向到 asyncio 程序：

```python
...
# start a subprocess and redirect output
process = await asyncio.create_subprocess_exec('ls', stdout=asyncio.subprocess.PIPE)
```

然后我们可以通过 asyncio.subprocess.Process 实例通过 communicate() 方法读取程序的输出。

此方法是协程，必须等待。它用于通过子流程发送和接收数据。

```python
...
# read data from the subprocess
line = process.communicate()
```

我们还可以通过以字节为单位设置“input”参数，通过 communicate() 方法将数据发送到子进程。

```python
...
# start a subprocess and redirect input
process = await asyncio.create_subprocess_exec('ls', stdin=asyncio.subprocess.PIPE)
# send data to the subprocess
process.communicate(input=b'Hello\n')
```

在后台，asyncio.subprocess.PIPE 将子进程配置为指向 StreamReader 或 StreamWriter，用于向子进程发送数据或从子进程发送数据，并且 communicate() 方法将从配置的读取器读取或写入字节。

我们可以通过子进程通过 stdin、stdout 和 stderr 属性直接与 StreamReader 或 StreamWriter 交互。

```python
...
# read a line from the subprocess output stream
line = await process.stdout.readline()
```

现在我们知道如何使用 create_subprocess_exec() 函数，让我们看一些工作示例。



### 2.2. Asyncio create_subprocess_exec() 示例

我们可以探索如何在 asyncio 的子进程中运行命令。在这个例子中，我们将执行“echo”命令来报告一个字符串。

echo 命令将直接在标准输出上报告提供的字符串。下面列出了完整的示例。

请注意，此示例假设您可以访问“echo”命令，我不确定它是否适用于 Windows。

```python
# SuperFastPython.com
# example of executing a command as a subprocess with asyncio
import asyncio
 
# main coroutine
async def main():
    # start executing a command in a subprocess
    process = await asyncio.create_subprocess_exec('echo', 'Hello World')
    # report the details of the subprocess
    print(f'subprocess: {process}')
 
# entry point
asyncio.run(main())
```

运行示例首先创建 main() 协程并将其作为 asyncio 程序的入口点执行。

main() 协程运行并调用 create_subprocess_exec() 函数来执行命令。

main() 协程在创建子进程时挂起。返回一个 Process 实例。

main() 协程恢复并报告子进程的详细信息。 main() 进程终止，asyncio 程序终止。

echo 命令的输出在命令行上报告。这突出了我们如何从 asyncio 程序执行命令。

```python
Hello World
subprocess: <Process 50249>
```



## 3. 如何通过 Shell 运行命令

我们可以使用 shell 执行命令。shell 是命令行的用户界面，称为命令行解释器 (CLI)。它将代表用户解释和执行命令。

它还提供诸如用于脚本、通配符、管道、shell 变量（例如 PATH）等的原始编程语言等功能。

例如，我们可以将一条命令的输出重定向为另一条命令的输入，比如将“/etc/services”文件的内容重定向到word count命令“wc”中，统计行数：

```sh
cat /etc/services | wc -l
```

基于 Unix 的操作系统中的 shell 示例包括：

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230217111714284.png)

shell 已经在运行，它被用来启动 Python 程序。您无需执行任何特殊操作即可获取或访问 shell。

我们可以通过 create_subprocess_shell() 函数从 asyncio 程序执行命令。

asyncio.create_subprocess_shell() 函数接受一个命令并使用当前用户 shell 执行它。

这很有用，因为它不仅允许执行命令，还允许使用 shell 的功能，例如重定向、通配符等。

该命令将在执行 asyncio 程序的进程的子进程中执行。重要的是，asyncio 程序能够与子进程异步交互，例如通过协程。

通过 shell 而不是直接执行命令时，可能会有安全考虑。

这是因为请求执行命令和正在执行的命令之间至少存在一层间接和解释，允许可能的恶意注入。

现在我们知道了 asyncio.create_subprocess_shell() 的作用，让我们看看如何使用它。



### 3.1. 如何使用 Asyncio create_subprocess_shell()

asyncio.create_subprocess_shell() 函数将通过当前 shell 执行给定的字符串命令。

它返回一个表示进程的 asyncio.subprocess.Process 对象。

它与我们在上一节中看到的 create_subprocess_shell() 函数非常相似。不过，我们将回顾如何使用该函数以及如何通过 Process 实例与流程交互（以防您直接跳到本节）。

create_subprocess_shell() 函数是一个协程，这意味着我们必须等待它。它会在子流程启动后返回，而不是在子流程完成时返回。

```python
...
# start a subprocess
process = await asyncio.create_subprocess_shell('ls')
```

我们可以通过等待 wait() 方法来等待子进程完成。

```python
...
# wait for the subprocess to terminate
await process.wait()
```

我们可以通过调用 terminate() 或 kill() 方法直接停止子进程，这将在子进程中引发一个信号。

命令的输入和输出将由 shell 处理，例如标准输入、标准错误和标准输出。

我们可以让 asyncio 程序处理子进程的输入或输出。

这可以通过指定输入或输出流并指定要重定向的常量来实现，例如 asyncio.subprocess.PIPE。

例如，我们可以将命令的输出重定向到 asyncio 程序：

```python
...
# start a subprocess and redirect output
process = await asyncio.create_subprocess_shell('ls', stdout=asyncio.subprocess.PIPE)
```

然后我们可以通过 asyncio.subprocess.Process 实例通过 communicate() 方法读取程序的输出。

此方法是协程，必须等待。它用于通过子流程发送和接收数据。

```python
...
# read data from the subprocess
line = process.communicate()
```

我们还可以通过以字节为单位设置“input”参数，通过 communicate() 方法将数据发送到子进程。

```python
...
# start a subprocess and redirect input
process = await asyncio.create_subprocess_shell('ls', stdin=asyncio.subprocess.PIPE)
# send data to the subprocess
process.communicate(input=b'Hello\n')
```

在后台，asyncio.subprocess.PIPE 将子进程配置为指向 StreamReader 或 StreamWriter，用于向子进程发送数据或从子进程发送数据，并且 communicate() 方法将从配置的读取器读取或写入字节。

我们可以通过子进程通过 stdin、stdout 和 stderr 属性直接与 StreamReader 或 StreamWriter 交互。

```python
...
# read a line from the subprocess output stream
line = await process.stdout.readline()
```

现在我们知道如何使用 create_subprocess_shell() 函数，让我们看一些工作示例。



### 3.2. Asyncio create_subprocess_shell() 示例

我们可以探索如何使用 shell 在 asyncio 的子进程中运行命令。在这个例子中，我们将执行“echo”命令来报告一个字符串。

echo 命令将直接在标准输出上报告提供的字符串。下面列出了完整的示例。

请注意，此示例假设您可以访问“echo”命令，我不确定它是否适用于 Windows。

```python
# SuperFastPython.com
# example of executing a shell command as a subprocess with asyncio
import asyncio
 
# main coroutine
async def main():
    # start executing a shell command in a subprocess
    process = await asyncio.create_subprocess_shell('echo Hello World')
    # report the details of the subprocess
    print(f'subprocess: {process}')
 
# entry point
asyncio.run(main())
```

运行示例首先创建 main() 协程并将其作为 asyncio 程序的入口点执行。main() 协程运行并调用 create_subprocess_shell() 函数来执行命令。

main() 协程运行并调用 create_subprocess_shell() 函数来执行命令。main() 协程恢复并报告子进程的详细信息。 main() 进程终止，asyncio 程序终止。

echo 命令的输出在命令行上报告。这突出显示了我们如何使用 shell 从 asyncio 程序执行命令。

```python
subprocess: <Process 43916>
Hello World
```

