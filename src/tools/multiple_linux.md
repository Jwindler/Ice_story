# 如何在多个 Linux 服务器上运行多个命令

如果你正在管理多台 Linux 服务器，并且你想在所有 Linux 服务器上运行多个命令，但你不知道该怎么做。不用担心，在这个简单的服务器管理[指南](https://www.tecmint.com/run-multiple-commands-on-multiple-linux-servers/ "Source")中，我们将向您展示如何在多个 Linux 服务器上同时运行多个命令。

为此，您可以使用 pssh（并行 ssh）程序，这是一个用于在多个主机上并行执行 ssh 的命令行实用程序。使用它，您可以从 shell 脚本向所有 ssh 进程发送输入。



## 创建 Shell 脚本

因此，您需要首先准备一个脚本，其中包含您要在不同服务器上执行的 Linux 命令。在此示例中，我们将编写一个脚本，该脚本将从多个服务器收集以下信息：

- 检查服务器的正常运行时间
- 检查谁登录以及他们在做什么
- 根据内存使用情况列出前 5 个正在运行的进程。

首先使用您喜欢的编辑器创建一个名为 commands.sh 的脚本。

```sh
# vi commands.sh
```

接下来，将以下命令添加到脚本中，如图所示。

```sh
#!/bin/bash 
###############################################################################
#Script Name    : commands.sh                       
#Description    : execute multiple commands on multiple servers                                                                     
#Author         : Aaron Kili Kisinga       
#Email          : aaronkilik@gmail.com 
################################################################################
echo
# show system uptime
uptime
echo
# show who is logged on and what they are doing
who
echo
# show top 5 processe by RAM usage 
ps -eo cmd,pid,ppid,%mem,%cpu --sort=-%mem | head -n 6

exit 0
```

保存文件并关闭它。然后如图所示使脚本可执行。

```sh
# chmod +x commands.sh
```



## 创建 PSSH 主机文件

接下来，在 hosts.txt 文件中添加要在其上运行命令的服务器列表，格式为 [user@]host[:port] 或仅提供服务器 IP 地址。

但我们建议您使用可以在 .ssh/config 文件中指定的 ssh 别名，如如何配置自定义 ssh 连接以简化远程访问中所述。

这种方法更有效和可靠，它允许您为每个远程服务器指定配置选项（如主机名、标识文件、端口、用户名等）。
以下是我们的示例 ssh 主机别名文件，也就是用户特定的 ssh 配置文件。

```sh
# vi ~/.ssh/config
```

![](https://s2.loli.net/2023/06/17/aQj9cABLoRhItbY.png)



接下来，创建一个 hosts.txt 文件，在这里您可以简单地指定别名（使用 .ssh/config 文件中的 Host 关键字定义的名称），如图所示。

```sh
# vi hosts.txt 
```

添加服务器别名。

```sh
server1
server2
server3
```



## 通过脚本在多个 Linux 服务器上运行命令

现在通过指定 hosts.txt 文件以及包含要在多个远程服务器上运行的多个命令的脚本来运行以下 pssh 命令。

```sh
# pssh -h hosts.txt -P -I<./commands.sh
```

上述命令中使用的标志的含义：

- `-h` – 读取主机文件。
- `-P` – 告诉 pssh 在输出到达时显示输出。
- `-I` – 读取输入并发送到每个 ssh 进程。

![](https://s2.loli.net/2023/06/17/DKAIfaNEBGyYPkS.png)