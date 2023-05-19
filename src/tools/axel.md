# Axel – 用于 Linux 的命令行文件下载加速器

如果您是那种喜欢下载和试用多个 Linux 发行版的人，我们相信您会张开双臂欢迎一个说到做到的下载加速器——一个按照其描述进行操作的下载加速器。

在[本指南](https://www.tecmint.com/axel-commandline-download-accelerator-for-linux/ "Source")中，我们将向您介绍 Axel，这是一个没有依赖项（除了 gcc 和 makeutils）的轻量级 wget 克隆。

![](https://s2.loli.net/2023/05/19/v3IHZUtbpO4jqyg.png)



虽然它的描述表明它特别适用于字节关键系统，但 axel 可以安装在任何地方并且不仅可以用于通过 HTTP/FTP 链接同时下载多个文件，还可以加快它们的速度。



## 安装

正如我们之前提到的，axel 不仅仅是另一个下载工具。它通过使用多个连接从目标检索文件来加速 HTTP 和 FTP 下载，也可以配置为使用多个镜像。

如果这还不足以让你有动力去尝试，让我们补充一下，axel 支持自动中止和恢复在给定时间段后无响应或不返回任何数据的连接。

此外，如果您有权这样做，您可以利用 axel 打开多个同时的 FTP 连接到一个服务器，以增加每个连接分配的带宽。

如果您不允许这样做或不确定，您可以改为打开多个连接到单独的服务器并同时从所有服务器下载。

最后但同样重要的是，axel 与其他 Linux 下载加速器的不同之处在于它在下载时将所有数据放在一个文件中，而不是将数据写入单独的文件并在稍后阶段加入它们。

在 CentOS/RHEL 8/7 中，您需要启用 EPEL 存储库才能安装 axel：

```sh
yum install epel-release
yum install axel
```

在 Fedora 中，它可以从默认存储库中获得。

```shell
yum install axel   
dnf install axel   [On Fedora 23+ releases]
```

在 Debian 及其衍生版如 Ubuntu、Linux Mint 中，可以直接使用 aptitude 安装 axel：

```sh
aptitude install axel
```

在 Arch Linux 和相关发行版（例如 Manjaro Linux 和 OpenSUSE Linux）上，您可以直接安装 axel：

```sh
sudo pacman -S axel       [On Arch/Manjaro]
sudo zypper install axel  [On OpenSUSE]
```



## 配置

您可以使用 /etc/axelrc 配置 axel 并在调用它时在命令行中传递更多所需的选项。配置文件有详细记录，但我们将在此处查看最有用的选项：

- reconnect_delay 是 axel 在再次尝试启动与服务器的新连接之前等待的秒数。
- max_speed 值以每秒字节数 (B/s) 为单位。考虑到可用带宽后，您可能希望将此变量设置为适当的值。这将帮助您防止 axel 在下载时消耗大量带宽。

> 重要提示：请注意，实际最大下载速率将取决于您的 Internet 连接

- num_connections 是 axel 将尝试启动的最大连接数。推荐值 (4) 对于大多数情况已经足够，主要是出于对其他 FTP 用户的尊重。请注意，某些服务器甚至可能不允许多个连接。
- connection_timeout 指示 axel 在尝试中止并自动恢复之前等待接收响应的秒数。
- http_proxy 允许您设置代理服务器，以防 HTTP_PROXY 环境变量尚未在系统范围内设置。此变量使用与 HTTP_PROXY (http://:PORT) 相同的格式。
- no_proxy 是本地域的列表，以逗号分隔，axel 不应尝试通过代理访问这些域。此设置是可选的。
- buffer_size 表示一次从所有当前连接读取的最大数量（以字节为单位）。
- verbose 允许您选择是否在屏幕上打印与下载相关的消息。如果您想禁用它，请将其设置为 0，如果您仍想看到消息，请将其设置为 1。
- 如果您有多个接口，interfaces 可以让您列出可以访问 Internet 的网络接口。如果未明确设置，axel 将使用路由表中的第一个接口。

如果仔细观察，您会发现大多数命令行选项与配置文件中的选项相似。此外，-o (–output) 选项允许您指定输出文件名。

如果使用，它将覆盖源文件名。如果您设置任何命令行选项，它们将覆盖配置文件中的设置。



## 使用

我们将使用配置文件中的以下设置（取消注释相应行）：

```sh
reconnect_delay = 20
max_speed = 500000
num_connections = 4
connection_timeout = 30
buffer_size = 10240
verbose = 1
```

![](https://s2.loli.net/2023/05/19/5g3zhWNVpQeZXfo.png)



我们现在将使用 wget 和 axel 比较 HTTP 和 FTP 链接的下载时间。您可以选择任何大小的任何文件，但为简单起见，我们将从以下位置下载 100 MB 的文件：

- ftp://speedtest:speedtest@ftp.otenet.gr/test100Mb.db
- http://speedtest.ftp.otenet.gr/files/test100Mb.db



### FTP

使用 wget 进行 FTP 下载（平均 459 KB/s）：

```sh
wget ftp://speedtest:speedtest@ftp.otenet.gr/test100Mb.db
```

![](https://s2.loli.net/2023/05/19/wvPexYQK4Bas8r2.png)



### axel

使用 axel 下载 FTP（平均 1181.43 KB/s）：

```sh
axel -n 10 --output=axel-test100Mb.db ftp://speedtest:speedtest@ftp.otenet.gr/test100Mb.db
```

![](https://s2.loli.net/2023/05/19/F4tmYbWfHZVQrBO.png)



正如您在我们上面执行的测试结果中看到的那样，axel 可以显着加速 FTP 或 HTTP 下载。



## 总结

在本文中，我们解释了如何使用 axel，一种 FTP/HTTP 下载加速器，并展示了它如何比 wget 等其他程序执行得更快，因为它能够同时打开多个到远程服务器的连接。