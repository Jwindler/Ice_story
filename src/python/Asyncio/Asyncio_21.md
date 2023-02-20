# Python 异步: 检查网站状态示例（21）

我们可以通过打开流并写入和读取 HTTP 请求和响应来使用 asyncio 查询网站的 HTTP 状态。

然后我们可以使用 asyncio 并发查询多个网站的状态，甚至动态报告结果。



## 1. 如何使用 Asyncio 检查 HTTP 状态

asyncio 模块提供了对打开套接字连接和通过流读写数据的支持。我们可以使用此功能来检查网页的状态。

这可能涉及四个步骤，它们是：

- 打开一个连接
- 写一个请求
- 读取响应
- 关闭连接



## 2. 打开 HTTP 连接

可以使用 asyncio.open_connection() 函数在 asyncio 中打开连接。在众多参数中，该函数采用字符串主机名和整数端口号。

这是一个必须等待的协程，它返回一个 StreamReader 和一个 StreamWriter，用于使用套接字进行读写。

这可用于在端口 80 上打开 HTTP 连接。

```python
...
# open a socket connection
reader, writer = await asyncio.open_connection('www.google.com', 80)
```

我们还可以使用 ssl=True 参数打开 SSL 连接。这可用于在端口 443 上打开 HTTPS 连接。

```python
...
# open a socket connection
reader, writer = await asyncio.open_connection('www.google.com', 443)
```



## 3. 写入 HTTP 请求

打开后，我们可以向 StreamWriter 写入查询以发出 HTTP 请求。例如，HTTP 版本 1.1 请求是纯文本格式的。我们可以请求文件路径“/”，它可能如下所示：

```python
GET / HTTP/1.1
Host: www.google.com
```

重要的是，每行末尾必须有一个回车和一个换行符（\r\n），末尾有一个空行。

作为 Python 字符串，这可能如下所示：

```python
'GET / HTTP/1.1\r\n'
'Host: www.google.com\r\n'
'\r\n'
```

在写入 StreamWriter 之前，此字符串必须编码为字节。这可以通过对字符串本身使用 encode() 方法来实现。默认的“utf-8”编码可能就足够了。

```python
...
# encode string as bytes
byte_data = string.encode()
```

然后可以通过 StreamWriter 的 write() 方法将字节写入套接字。

```python
...
# write query to socket
writer.write(byte_data)
```

写入请求后，最好等待字节数据发送完毕并等待套接字准备就绪。这可以通过 drain() 方法来实现。这是一个必须等待的协程。

```python
...
# wait for the socket to be ready.
await writer.drain()
```



## 4. 读取 HTTP 响应

发出 HTTP 请求后，我们可以读取响应。这可以通过套接字的 StreamReader 来实现。可以使用读取一大块字节的 read() 方法或读取一行字节的 readline() 方法来读取响应。

我们可能更喜欢 readline() 方法，因为我们使用的是基于文本的 HTTP 协议，它一次发送一行 HTML 数据。readline() 方法是协程，必须等待。

```python
...
# read one line of response
line_bytes = await reader.readline()
```

HTTP 1.1 响应由两部分组成，一个由空行分隔的标头，然后是一个空行终止的主体。header 包含有关请求是否成功以及将发送什么类型的文件的信息，body 包含文件的内容，例如 HTML 网页。

HTTP 标头的第一行包含服务器上所请求页面的 HTTP 状态。每行都必须从字节解码为字符串。

这可以通过对字节数据使用 decode() 方法来实现。同样，默认编码为“utf_8”。

```python
...
# decode bytes into a string
line_data = line_bytes.decode()
```



## 5. 关闭 HTTP 连接

我们可以通过关闭 StreamWriter 来关闭套接字连接。这可以通过调用 close() 方法来实现。

```python
...
# close the connection
writer.close()
```

这不会阻塞并且可能不会立即关闭套接字。现在我们知道如何使用 asyncio 发出 HTTP 请求和读取响应，让我们看一些检查网页状态的示例。



## 6. 顺序检查 HTTP 状态的示例

我们可以开发一个示例来使用 asyncio 检查多个网站的 HTTP 状态。

在此示例中，我们将首先开发一个协程来检查给定 URL 的状态。然后我们将为排名前 10 的网站中的每一个调用一次这个协程。

首先，我们可以定义一个协程，它将接受一个 URL 字符串并返回 HTTP 状态。

```python
# get the HTTP/S status of a webpage
async def get_status(url):
	# ...
```

必须将 URL 解析为其组成部分。我们在发出 HTTP 请求时需要主机名和文件路径。我们还需要知道 URL 方案（HTTP 或 HTTPS）以确定是否需要 SSL。

这可以使用 urllib.parse.urlsplit() 函数来实现，该函数接受一个 URL 字符串并返回所有 URL 元素的命名元组。

```python
...
# split the url into components
url_parsed = urlsplit(url)
```

然后我们可以打开基于 URL 方案的 HTTP 连接并使用 URL 主机名。

```python
...
# open the connection
if url_parsed.scheme == 'https':
    reader, writer = await asyncio.open_connection(url_parsed.hostname, 443, ssl=True)
else:
    reader, writer = await asyncio.open_connection(url_parsed.hostname, 80)
```

接下来，我们可以使用主机名和文件路径创建 HTTP GET 请求，并使用 StreamWriter 将编码字节写入套接字。

```python
...
# send GET request
query = f'GET {url_parsed.path} HTTP/1.1\r\nHost: {url_parsed.hostname}\r\n\r\n'
# write query to socket
writer.write(query.encode())
# wait for the bytes to be written to the socket
await writer.drain()
```

接下来，我们可以读取 HTTP 响应。我们只需要包含 HTTP 状态的响应的第一行。

```python
...
# read the single line response
response = await reader.readline()
```

然后可以关闭连接。

```python
...
# close the connection
writer.close()
```

最后，我们可以解码从服务器读取的字节、远程尾随空白，并返回 HTTP 状态。

```python
...
# decode and strip white space
status = response.decode().strip()
# return the response
return status
```

将它们结合在一起，下面列出了完整的 get_status() 协程。它没有任何错误处理，例如无法访问主机或响应缓慢的情况。这些添加将为读者提供一个很好的扩展。

```python
# get the HTTP/S status of a webpage
async def get_status(url):
    # split the url into components
    url_parsed = urlsplit(url)
    # open the connection
    if url_parsed.scheme == 'https':
        reader, writer = await asyncio.open_connection(url_parsed.hostname, 443, ssl=True)
    else:
        reader, writer = await asyncio.open_connection(url_parsed.hostname, 80)
    # send GET request
    query = f'GET {url_parsed.path} HTTP/1.1\r\nHost: {url_parsed.hostname}\r\n\r\n'
    # write query to socket
    writer.write(query.encode())
    # wait for the bytes to be written to the socket
    await writer.drain()
    # read the single line response
    response = await reader.readline()
    # close the connection
    writer.close()
    # decode and strip white space
    status = response.decode().strip()
    # return the response
    return status
```

接下来，我们可以为我们要检查的多个网页或网站调用 get_status() 协程。在这种情况下，我们将定义一个世界排名前 10 的网页列表。

```python
...
# list of top 10 websites to check
sites = ['https://www.google.com/',
    'https://www.youtube.com/',
    'https://www.facebook.com/',
    'https://twitter.com/',
    'https://www.instagram.com/',
    'https://www.baidu.com/',
    'https://www.wikipedia.org/',
    'https://yandex.ru/',
    'https://yahoo.com/',
    'https://www.whatsapp.com/'
    ]
```

然后我们可以使用我们的 get_status() 协程依次查询每个。在这种情况下，我们将在一个循环中按顺序这样做，并依次报告每个状态。

```python
...
# check the status of all websites
for url in sites:
    # get the status for the url
    status = await get_status(url)
    # report the url and its status
    print(f'{url:30}:\t{status}')
```

在使用 asyncio 时，我们可以做得比顺序更好，但这提供了一个很好的起点，我们可以在以后进行改进。将它们结合在一起，main() 协程查询前 10 个网站的状态。

```python
# main coroutine
async def main():
    # list of top 10 websites to check
    sites = ['https://www.google.com/',
        'https://www.youtube.com/',
        'https://www.facebook.com/',
        'https://twitter.com/',
        'https://www.instagram.com/',
        'https://www.baidu.com/',
        'https://www.wikipedia.org/',
        'https://yandex.ru/',
        'https://yahoo.com/',
        'https://www.whatsapp.com/'
        ]
    # check the status of all websites
    for url in sites:
        # get the status for the url
        status = await get_status(url)
        # report the url and its status
        print(f'{url:30}:\t{status}')
```

最后，我们可以创建 main() 协程并将其用作 asyncio 程序的入口点。

```python
...
# run the asyncio program
asyncio.run(main())
```

将它们结合在一起，下面列出了完整的示例。

```python
# SuperFastPython.com
# check the status of many webpages
import asyncio
from urllib.parse import urlsplit
 
# get the HTTP/S status of a webpage
async def get_status(url):
    # split the url into components
    url_parsed = urlsplit(url)
    # open the connection
    if url_parsed.scheme == 'https':
        reader, writer = await asyncio.open_connection(url_parsed.hostname, 443, ssl=True)
    else:
        reader, writer = await asyncio.open_connection(url_parsed.hostname, 80)
    # send GET request
    query = f'GET {url_parsed.path} HTTP/1.1\r\nHost: {url_parsed.hostname}\r\n\r\n'
    # write query to socket
    writer.write(query.encode())
    # wait for the bytes to be written to the socket
    await writer.drain()
    # read the single line response
    response = await reader.readline()
    # close the connection
    writer.close()
    # decode and strip white space
    status = response.decode().strip()
    # return the response
    return status
 
# main coroutine
async def main():
    # list of top 10 websites to check
    sites = ['https://www.google.com/',
        'https://www.youtube.com/',
        'https://www.facebook.com/',
        'https://twitter.com/',
        'https://www.instagram.com/',
        'https://www.baidu.com/',
        'https://www.wikipedia.org/',
        'https://yandex.ru/',
        'https://yahoo.com/',
        'https://www.whatsapp.com/'
        ]
    # check the status of all websites
    for url in sites:
        # get the status for the url
        status = await get_status(url)
        # report the url and its status
        print(f'{url:30}:\t{status}')
 
# run the asyncio program
asyncio.run(main())
```

运行示例首先创建 main() 协程并将其用作程序的入口点。main() 协程运行，定义前 10 个网站的列表。然后顺序遍历网站列表。 main()协程挂起调用get_status()协程查询一个网站的状态。

get_status() 协程运行、解析 URL 并打开连接。它构造一个 HTTP GET 查询并将其写入主机。读取、解码并返回响应。main() 协程恢复并报告 URL 的 HTTP 状态。

对列表中的每个 URL 重复此操作。该程序大约需要 5.6 秒才能完成，或者平均每个 URL 大约需要半秒。这突出了我们如何使用 asyncio 来查询网页的 HTTP 状态。

尽管如此，它并没有充分利用 asyncio 来并发执行任务。

```python
https://www.google.com/       :	HTTP/1.1 200 OK
https://www.youtube.com/      :	HTTP/1.1 200 OK
https://www.facebook.com/     :	HTTP/1.1 302 Found
https://twitter.com/          :	HTTP/1.1 200 OK
https://www.instagram.com/    :	HTTP/1.1 200 OK
https://www.baidu.com/        :	HTTP/1.1 200 OK
https://www.wikipedia.org/    :	HTTP/1.1 200 OK
https://yandex.ru/            :	HTTP/1.1 302 Moved temporarily
https://yahoo.com/            :	HTTP/1.1 301 Moved Permanently
https://www.whatsapp.com/     :	HTTP/1.1 302 Found
```



## 7. 并发查看网站状态示例

asyncio 的一个好处是我们可以同时执行许多协程。我们可以使用 asyncio.gather() 函数在 asyncio 中并发查询网站的状态。

此函数采用一个或多个协程，暂停执行提供的协程，并将每个协程的结果作为可迭代对象返回。然后我们可以遍历 URL 列表和可迭代的协程返回值并报告结果。

这可能是比上述方法更简单的方法。首先，我们可以创建一个协程列表。

```python
...
# create all coroutine requests
coros = [get_status(url) for url in sites]
```

接下来，我们可以执行协程并使用 asyncio.gather() 获取可迭代的结果。

请注意，我们不能直接提供协程列表，而是必须将列表解压缩为单独的表达式，这些表达式作为位置参数提供给函数。

```python
...
# execute all coroutines and wait
results = await asyncio.gather(*coros)
```

这将同时执行所有协程并检索它们的结果。然后我们可以遍历 URL 列表和返回状态并依次报告每个。

```python
...
# process all results
for url, status in zip(sites, results):
    # report status
    print(f'{url:30}:\t{status}')
```

将它们结合在一起，下面列出了完整的示例。

```python
# SuperFastPython.com
# check the status of many webpages
import asyncio
from urllib.parse import urlsplit
 
# get the HTTP/S status of a webpage
async def get_status(url):
    # split the url into components
    url_parsed = urlsplit(url)
    # open the connection
    if url_parsed.scheme == 'https':
        reader, writer = await asyncio.open_connection(url_parsed.hostname, 443, ssl=True)
    else:
        reader, writer = await asyncio.open_connection(url_parsed.hostname, 80)
    # send GET request
    query = f'GET {url_parsed.path} HTTP/1.1\r\nHost: {url_parsed.hostname}\r\n\r\n'
    # write query to socket
    writer.write(query.encode())
    # wait for the bytes to be written to the socket
    await writer.drain()
    # read the single line response
    response = await reader.readline()
    # close the connection
    writer.close()
    # decode and strip white space
    status = response.decode().strip()
    # return the response
    return status
 
# main coroutine
async def main():
    # list of top 10 websites to check
    sites = ['https://www.google.com/',
        'https://www.youtube.com/',
        'https://www.facebook.com/',
        'https://twitter.com/',
        'https://www.instagram.com/',
        'https://www.baidu.com/',
        'https://www.wikipedia.org/',
        'https://yandex.ru/',
        'https://yahoo.com/',
        'https://www.whatsapp.com/'
        ]
    # create all coroutine requests
    coros = [get_status(url) for url in sites]
    # execute all coroutines and wait
    results = await asyncio.gather(*coros)
    # process all results
    for url, status in zip(sites, results):
        # report status
        print(f'{url:30}:\t{status}')
 
# run the asyncio program
asyncio.run(main())
```

运行该示例会像以前一样执行 main() 协程。在这种情况下，协程列表是在列表理解中创建的。

然后调用 asyncio.gather() 函数，传递协程并挂起 main() 协程，直到它们全部完成。协程执行，同时查询每个网站并返回它们的状态。

main() 协程恢复并接收可迭代的状态值。然后使用 zip() 内置函数遍历此可迭代对象和 URL 列表，并报告状态。

这突出了一种更简单的方法来同时执行协程并在所有任务完成后报告结果。它也比上面的顺序版本更快，在我的系统上完成大约 1.4 秒。

```python
https://www.google.com/       :	HTTP/1.1 200 OK
https://www.youtube.com/      :	HTTP/1.1 200 OK
https://www.facebook.com/     :	HTTP/1.1 302 Found
https://twitter.com/          :	HTTP/1.1 200 OK
https://www.instagram.com/    :	HTTP/1.1 200 OK
https://www.baidu.com/        :	HTTP/1.1 200 OK
https://www.wikipedia.org/    :	HTTP/1.1 200 OK
https://yandex.ru/            :	HTTP/1.1 302 Moved temporarily
https://yahoo.com/            :	HTTP/1.1 301 Moved Permanently
https://www.whatsapp.com/     :	HTTP/1.1 302 Found
```

