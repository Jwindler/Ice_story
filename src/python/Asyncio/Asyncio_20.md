# Python 异步: 非阻塞流（20）

asyncio 的一个主要好处是能够使用非阻塞流。



## 1. 异步流

Asyncio 提供非阻塞 I/O 套接字编程。这是通过流提供的。

可以打开提供对流写入器和流写入器的访问的套接字。然后可以使用协同程序从流中写入和读取数据，并在适当的时候暂停。完成后，可以关闭套接字。

异步流功能是低级的，这意味着必须手动实现所需的任何协议。

这可能包括常见的 Web 协议，例如：

- 用于与 Web 服务器交互的 HTTP 或 HTTPS
- 用于与电子邮件服务器交互的 SMTP
- 用于与文件服务器交互的 FTP。

这些流还可用于创建服务器以使用标准协议处理请求，或开发您自己的特定于应用程序的协议。

现在我们知道什么是异步流，让我们看看如何使用它们。



## 2. 如何打开连接

可以使用 asyncio.open_connection() 函数打开 asyncio TCP 客户端套接字连接。

这是一个必须等待的协程，一旦套接字连接打开就会返回。

该函数返回用于与套接字交互的 StreamReader 和 StreamWriter 对象。

```python
...
# open a connection
reader, writer = await asyncio.open_connection(...)
```

asyncio.open_connection() 函数采用许多参数来配置套接字连接。两个必需的参数是主机和端口。

host 是一个字符串，指定要连接的服务器，例如域名或 IP 地址。

port为socket端口号，如HTTP服务器为80，HTTPS服务器为443，SMTP为23等。

```python
...
# open a connection to an http server
reader, writer = await asyncio.open_connection('www.google.com', 80)
```

通过 SSL 协议支持加密套接字连接。最常见的例子可能是正在取代 HTTP 的 HTTPS。这可以通过将“ssl”参数设置为 True 来实现。

```python

...
# open a connection to an https server
reader, writer = await asyncio.open_connection('www.google.com', 443, ssl=True)
```



## 3. 如何启动服务器

可以使用 asyncio.start_server() 函数打开 asyncio TCP 服务器套接字。这是一个必须等待的协程。

该函数返回一个代表正在运行的服务器的 asyncio.Server 对象。

```python
...
# start a tcp server
server = await asyncio.start_server(...)
```

三个必需的参数是回调函数、主机和端口。回调函数是一个由名称指定的自定义函数，每次客户端连接到服务器时都会调用该函数。

主机是客户端将指定连接的域名或 IP 地址。端口是接收连接的套接字端口号，例如 21 用于 FTP 或 80 用于 HTTP。

```python
# handle connections
async def handler(reader, writer):
	# ...
 
...
# start a server to receive http connections
server = await asyncio.start_server(handler, '127.0.0.1', 80)
```



## 4. 如何使用 StreamWriter 写入数据

我们可以使用 asyncio.StreamWriter 将数据写入套接字。数据以字节形式写入。可以使用 write() 方法将字节数据写入套接字。

```python
...
# write byte data
writer.write(byte_data)
```

或者，可以使用 writelines() 方法写入组织成列表或可迭代的多“行”字节数据。

```python
...
# write lines of byte data
writer.writelines(byte_lines)
```

写入数据块或挂起调用协程的方法都没有。写入字节数据后，最好通过 drain() 方法清空套接字。这是一个Coroutine，将暂停呼叫者，直到传输字节并准备好插座为止。

```python
...
# write byte data
writer.write(byte_data)
# wait for data to be transmitted
await writer.drain()
```



## 5. 如何使用 StreamReader 读取数据

我们可以使用 asyncio.StreamReader 从套接字读取数据。数据以字节格式读取，因此字符串在使用前可能需要进行编码。所有读取方法都是必须等待的协程。

可以通过 read() 方法读取任意数量的字节，该方法将一直读取到文件末尾 (EOF)。

```python
...
# read byte data
byte_data = await reader.read()
```

此外，可以通过“n”参数指定要读取的字节数。如果您知道下一个响应的预期字节数，这可能会有所帮助。

```python
...
# read byte data
byte_data = await reader.read(n=100)
```

可以使用 readline() 方法读取单行数据。这将返回字节，直到遇到换行符“\n”或 EOF。

这在阅读使用文本行操作的标准协议时很有用。

```python
...
# read a line data
byte_line = await reader.readline()
```

此外，还有一个 readexactly() 方法来读取确切数量的字节，否则会引发异常，还有一个 readuntil() 方法将读取字节，直到读取字节形式的指定字符。



## 6. 如何关闭连接

可以通过 asyncio.StreamWriter 关闭套接字。可以调用 close() 方法来关闭套接字。此方法不会阻塞。

```python
...
# close the socket
writer.close()
```

虽然 close() 方法不会阻塞，但我们可以等待套接字完全关闭后再继续。这可以通过 wait_closed() 方法来实现。这是一个可以等待的协程。

```python
...
# close the socket
writer.close()
# wait for the socket to close
await writer.wait_closed()
```

我们可以通过 is_closing() 方法检查套接字是否已经关闭或正在关闭。

```python
...
# check if the socket is closed or closing
if writer.is_closing():
	# ...
```

