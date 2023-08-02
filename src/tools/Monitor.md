# 如何一目了然地监控远程 Linux 系统

Glances 是一款免费的开源、现代、跨平台、实时 top 和类似 htop 的系统监控工具，与同类工具相比，它提供了先进的功能，并且可以在不同的模式下运行：作为独立模式、客户端/服务器模式，并在 Web 服务器模式下。

考虑到Web服务器模式，您不一定需要通过SSH登录远程服务器来运行glances，您可以在Web服务器模式下运行它并通过Web浏览器访问它来远程监控您的Linux服务器，如下所述。

要在 Web 服务器模式下运行 Glance，您需要使用适合您的 Linux 发行版的命令，将其与 Python Bottle 模块（一个快速、简单且轻量级的 WSGI 微型 Web 框架）一起安装。



```sh
$ sudo apt install glances python-bottle         [On Debian, Ubuntu and Mint]
$ sudo yum install glances python-bottle         [On RHEL/CentOS/Fedora and Rocky/AlmaLinux]
$ sudo emerge -a sys-apps/glances python-bottle  [On Gentoo Linux]
$ sudo apk add glances python-bottle             [On Alpine Linux]
$ sudo pacman -S glances python-bottle           [On Arch Linux]
$ sudo zypper install glances python-bottle      [On OpenSUSE]    
```

或者，使用所示的 PIP 命令安装它。

```sh
$ sudo pip install bottle
```

安装上述软件包后，使用 -w 标志启动 Glaces 以在 Web 服务器模式下运行它。默认情况下，它将侦听端口 61208。

```sh
$ glances -w 
OR
$ glances -w &
```

如果您正在运行firewalld服务，那么您应该打开端口61208以允许入站流量到达该端口。

```sh
$ sudo firewall-cmd --permanent --add-port=61208/tcp
$ sudo firewall-cmd --reload
```

对于 UFW 防火墙，运行以下命令。

```sh
$ sudo ufw allow 61208/tcp
$ sudo ufw reload
```

之后，从 Web 浏览器中使用 URL http://SERVER_IP:61208/ 访问 Glances UI。

如果您使用 systemd 系统和服务管理器，则可以在 Web 服务器模式下将 Glas 作为一项服务运行，以实现高效管理，如下一节所述。实际上我更喜欢这种方法作为后台进程运行。



## 在 Web 服务器模式下将 Glance 作为服务运行

首先在/usr/lib/systemd/system/glancesweb.service 下创建服务单元文件。

```sh
$ sudo vim /usr/lib/systemd/system/glancesweb.service
```

然后将下面的单元文件配置复制并粘贴到其中。

```sh
[Unit]
Description = Glances in Web Server Mode
After = network.target

[Service]
ExecStart = /usr/bin/glances  -w  -t  5

[Install]
WantedBy = multi-user.target
```

上面的配置告诉systemd这是一个unit-of-type服务，它应该在network.target之后加载。

一旦系统位于网络目标中，systemd 将调用命令“/usr/bin/glances -w -t 5”作为服务。 -t 指定实时更新的时间间隔（以秒为单位）。

[install] 部分通知 systemd “multi-user.target” 需要此服务。因此，当您启用它时，会创建一个从 

/etc/systemd/system/multi-user.target.wants/glancesweb.service 到 

/usr/lib/systemd/system/glancesweb.service 的符号链接。禁用它将删除该符号链接。

接下来，启用新的 systemd 服务，启动并查看其状态，如下所示。

```sh
$ sudo systemctl enable glancesweb.service
$ sudo systemctl start glancesweb.service
$ sudo systemctl status glancesweb.service
```

最后，在您的 Web 浏览器中，使用 URL http://SERVER_IP:61208/ 在任何设备（智能手机、平板电脑或计算机）上通过 Glances UI 远程监控您的 Linux 服务器。

![](https://s2.loli.net/2023/08/02/qYrPodlEcjUeLyf.png)



![](https://s2.loli.net/2023/08/02/jxhwRGTy4lHYzEv.png)



您可以更改页面的刷新率，只需在 URL 末尾添加以秒为单位的句点，这会将刷新率设置为 8 秒。

```sh
http://SERVERI_P:61208/8	
```

在 Web 服务器模式下运行 Glance 的一个缺点是，如果 Internet 连接较差，客户端很容易与服务器断开连接。

您可以从[本指南](https://www.tecmint.com/glances-monitor-remote-linux-systems/ "Source")中了解如何创建新的 systemd 服务：如何在 Linux 中创建 Systemd 单元文件

