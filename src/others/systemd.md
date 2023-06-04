# 如何在 Linux 中列出 Systemd 下所有正在运行的服务

Linux系统提供多种系统服务（如进程管理、登录、syslog、cron等）和网络服务（如远程登录、电子邮件、打印机、虚拟主机、数据存储、文件传输、域名解析等） （使用 DNS）、动态 IP 地址分配（使用 DHCP）等等）。

从技术上讲，服务是一个进程或一组进程（通常称为守护进程）在后台持续运行，等待请求进入（尤其是来自客户端的请求）。

Linux 支持不同的方式来管理（启动、停止、重新启动、在系统启动时启用自动启动等）服务，通常是通过进程或服务管理器。大多数（如果不是全部）现代 Linux 发行版现在都使用相同的进程管理器：systemd。

Systemd 是 Linux 的系统和服务管理器； init 进程的替代品，与 SysV 和 LSB init 脚本兼容，systemctl 命令是管理 systemd 的主要工具。

在[本指南](https://www.tecmint.com/list-all-running-services-under-systemd-in-linux/ "Source")中，我们将演示如何在 Linux 中列出 systemd 下所有正在运行的服务。



## 在 Linux 中列出 SystemD 下正在运行的服务

当您运行不带任何参数的 systemctl 命令时，它将显示所有加载的 systemd 单元的列表（阅读 systemd 文档以获取有关 systemd 单元的更多信息），包括服务，显示它们的状态（无论是否处于活动状态）。

```sh
systemctl 
```

![](https://s2.loli.net/2023/06/04/8QXGlkHfNWsbzIC.png)



要列出系统上所有已加载的服务（无论是活动的、正在运行的、退出的还是失败的，请使用 list-units 子命令和带有服务值的 --type 开关。

```sh
# systemctl list-units --type=service
OR
# systemctl --type=service
```

![](https://s2.loli.net/2023/06/04/1U3nDFpazyHmeMk.png)



要列出所有已加载但处于活动状态的服务，包括正在运行的和已退出的服务，您可以添加值为 active 的 --state 选项，如下所示。

```sh
# systemctl list-units --type=service --state=active
OR
# systemctl --type=service --state=active
```

![](https://s2.loli.net/2023/06/04/MU4NhWJXvQ9cFk3.png)



但要快速浏览所有正在运行的服务（即所有已加载和正在运行的服务），请运行以下命令。

```sh
# systemctl list-units --type=service --state=running 
OR
# systemctl --type=service --state=running
```

![](https://s2.loli.net/2023/06/04/oEpFKiw2DAJum4a.png)



如果您经常使用前面的命令，您可以如图所示在 ~/.bashrc 文件中创建一个别名命令，以便轻松调用它。

```sh
vim ~/.bashrc
```

然后在别名列表下添加以下行，如屏幕截图所示。

```sh
alias running_services='systemctl list-units  --type=service  --state=running'
```

![](https://s2.loli.net/2023/06/04/tumaKNLzJqAhvnU.png)



保存文件中的更改并关闭它。从现在开始，使用“running_services”命令查看服务器上所有已加载、正在运行的服务的列表。

```sh
# running_services	#use the Tab completion 
```

![](https://s2.loli.net/2023/06/04/udGBaYFzXerK5xO.png)



此外，服务的一个重要方面是它们使用的端口。要确定守护进程正在侦听的端口，您可以使用 netstat 或 ss 命令，如图所示。

其中标志 -l 表示打印所有侦听套接字，-t 显示所有 TCP 连接，-u 显示所有 UDP 连接，-n 表示打印数字端口号（而不是应用程序名称），-p 表示显示应用程序名称。

```sh
# netstat -ltup | grep zabbix_agentd
OR
# ss -ltup | grep zabbix_agentd
```

第五列显示套接字：Local Address:Port。在这种情况下，进程 zabbix_agentd 正在侦听端口 10050。

![](https://s2.loli.net/2023/06/04/qalkWGd2ow5DyB4.png)



此外，如果您的服务器正在运行防火墙服务，该服务控制如何阻止或允许进出所选服务或端口的流量，您可以使用 firewall-cmd 或 ufw 命令列出已在防火墙中打开的服务或端口（取决于您使用的 Linux 发行版），如图所示。

```sh
# firewall-cmd --list-services   [FirewallD]
# firewall-cmd --list-ports

$ sudo ufw status     [UFW Firewall]
```

![](https://s2.loli.net/2023/06/04/L3oXGx2zcnECr4J.png)



目前为止就这样了！在本指南中，我们演示了如何在 Linux 中查看 systemd 下正在运行的服务。我们还介绍了如何检查正在侦听的端口服务以及如何查看在系统防火墙中打开的服务或端口。