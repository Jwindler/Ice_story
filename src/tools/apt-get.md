# Ubuntu 包管理的 20 个“apt-get”命令

在引入 apt 命令之前，apt-get 命令是基于 Debian 的 Linux 发行版中使用的主要包管理命令。

使用 apt-get 命令，您可以在系统上安装、删除、升级、搜索和管理软件包。然而，从 Ubuntu 16.04 和 Debian 9 开始，apt 命令成为推荐的包管理命令行工具，尽管 apt-get 仍然可用且功能齐全。



## 什么是 apt-get 命令？

apt-get命令是一个功能强大且免费的包管理命令行程序，用于与Ubuntu的APT（高级打包工具）库配合执行新软件包的安装、删除现有软件包、升级现有软件包、甚至用于升级整个操作系统。

apt-get 命令的语法如下：

```sh
sudo apt-get <options> <command>
```

这里，<options> 表示您可以与该命令一起使用的任何其他标志或修饰符，<command> 指定您要执行的操作，例如安装、升级、删除或搜索包。



## 1. 更新Ubuntu系统包

“update”命令用于从 /etc/apt/sources.list 文件中指定的源重新同步包索引文件。更新命令从其位置获取包并将包更新到较新的版本。

```sh
sudo apt-get update
```

![](https://s2.loli.net/2023/07/07/5PF6f2OicA8gmSy.png)



## 2. 升级Ubuntu系统包

‘upgrade’命令用于升级系统上当前安装的所有软件包。在任何情况下，都不会删除当前安装的软件包，也不会检索或安装尚未安装的软件包来满足升级依赖性。

```sh
sudo apt-get upgrade
```

![](https://s2.loli.net/2023/07/07/GYlfzEdF2u1KS6B.png)



## 3. 安装软件包

“install”子命令由一个或多个希望从存储库安装或升级的包进行跟踪。例如，要安装或更新名为 wget 的包，您可以运行：

```sh
sudo apt-get install wget
```

![](https://s2.loli.net/2023/07/07/bKqJ7A1S2FgVGzx.png)



或者，您还可以使用 apt-cache 命令在安装之前根据给定的搜索词（例如名称或描述）在系统包缓存中搜索包。



## 4. 安装多个软件包

您可以在命令中添加多个软件包名称，以便同时安装多个软件包。例如，以下命令将安装软件包“nethogs”和“goaccess”。

```sh
sudo apt-get install nethogs goaccess
```

![](https://s2.loli.net/2023/07/07/TWYjpQ6tLCPAsoy.png)



## 5. 使用通配符安装多个软件包

借助正则表达式，您可以使用一个字符串添加多个包。例如，我们使用 * 通配符来安装多个包含“*name*”字符串的软件包，名称将为“package-name”。

```sh
sudo apt-get install '*name*'
```

![](https://s2.loli.net/2023/07/07/gwLvh23ziMSUaQc.png)



## 6. 安装包而不升级

使用子“--no-upgrade”命令将阻止已安装的软件包升级。

```sh
sudo apt-get install packageName --no-upgrade
```

![](https://s2.loli.net/2023/07/07/I7OpoWgqdY8RfrV.png)



## 7. 更新单个包

“--only-upgrade”命令不会安装新的软件包，而只会升级已安装的软件包并禁用新安装的软件包。

```sh
sudo apt-get install packageName --only-upgrade
```

![](https://s2.loli.net/2023/07/07/4KUHAzsXg2wRYf6.png)



## 8. 安装特定的软件包版本

假设您希望仅安装特定版本的软件包，只需将“=”与软件包名称一起使用并附加所需的版本即可。

```sh
sudo apt-get install vsftpd=3.0.5-0ubuntu1
```

![](https://s2.loli.net/2023/07/07/l8OCcGEUHjbY7Iy.png)



## 9. 卸载不带配置的包

要卸载软件包而不删除其配置文件（以便以后重新使用相同的配置），请使用删除命令，如下所示。

```sh
sudo apt-get remove vsftpd
```

![](https://s2.loli.net/2023/07/07/XqJOtE4F8RZYW3x.png)



## 10. 完全删除带有配置的包

要删除软件包及其配置文件，请使用“purge”子命令，如下所示。

```sh
sudo apt-get purge vsftpd
```

![](https://s2.loli.net/2023/07/07/lgaIvh3RXMTqbL9.png)



或者，您可以将这两个命令组合在一起，如下所示。



```sh
sudo apt-get remove --purge vsftpd
```



## 11. 清除 Apt 缓存以节省磁盘空间

“clean”命令用于通过清理从本地存储库检索（下载）的 .deb 文件（包）来释放磁盘空间。

```sh
sudo apt-get autoclean
```

![](https://s2.loli.net/2023/07/07/Qd7y52fZGT6VHpP.png)



## 12. 下载软件包的源代码

要仅下载特定包的源代码，请使用选项“--download-only source”和“package-name”，如图所示。

```sh
sudo apt-get --download-only source vsftpd
```

![](https://s2.loli.net/2023/07/07/aWwL4Hbe8uZ3AKf.png)



## 13. 下载并解压源码包

要将包的源代码下载并解压到特定目录，请键入以下命令。

```sh
sudo apt-get source vsftpd
```

![](https://s2.loli.net/2023/07/07/tweo4TPHxc6RM8J.png)



当尝试从存储库下载包的源代码时，您可能会遇到一个常见错误“E：您必须将一些‘deb-src’URI 放入您的sources.list 中”。



## 14. 从源代码编译 Ubuntu 软件包

您还可以使用选项“--compile”同时下载、解压和编译源代码，如下所示。

```sh
sudo apt-get --compile source goaccess
```

![](https://s2.loli.net/2023/07/07/XBE48D7N1VZx3Uu.png)



## 15. 下载包而不安装

使用“下载”选项，您可以下载任何给定的包而无需安装它。例如，以下命令只会将“nethogs”包下载到当前工作目录。

```sh
sudo apt-get download nethogs
```

![](https://s2.loli.net/2023/07/07/I1SkKqy6vMdPAYa.png)



## 16. 查看软件包变更日志

“changelog”标志下载软件包更改日志并显示已安装的软件包版本。

```sh
sudo apt-get changelog vsftpd
```

![](https://s2.loli.net/2023/07/07/blDpg8ufAGiKWEM.png)



## 17. 查看 Ubuntu 中损坏的依赖关系

“check”命令是一个诊断工具，用于更新包缓存并检查损坏的依赖项。

```sh
sudo apt-get check
```

![](https://s2.loli.net/2023/07/07/W6yKPJDQXnlEHmx.png)



## 18. 安装包的构建依赖项

‘build-dep’命令搜索系统中的本地存储库并安装curl包的构建依赖项。如果本地存储库中不存在该包，它将返回错误代码。

```sh
sudo apt-get build-dep curl
```

![](https://s2.loli.net/2023/07/07/EcMzg1PCLGiDJUk.png)



## 19. 自动删除已安装的软件包

“autoremove”子命令用于自动删除某些软件包，这些软件包本来是为了满足其他软件包的依赖关系而安装的，但现在不再需要了。例如，以下命令将删除已安装的软件包及其依赖项。

```sh
sudo apt-get autoremove vsftpd
```

![](https://s2.loli.net/2023/07/07/e4Z1SglCuAtPXnr.png)



## 20. apt-get 命令帮助

apt-get help 命令显示内置帮助文档，以及与 apt-get 命令一起使用的可用选项。

```sh
sudo apt-get help
```

![](https://s2.loli.net/2023/07/07/KFabyGP4J2hfcqC.png)



我已经使用 apt-get 命令介绍了大部分可用选项，但仍然有更多可用选项，您可以从终端使用“man apt-get”查看它们。

我希望您喜欢阅读[这篇文章](https://www.tecmint.com/apt-get-command/ "Source")，如果我遗漏了任何内容并且您希望我添加到列表中。请随时在下面的评论中提及这一点。