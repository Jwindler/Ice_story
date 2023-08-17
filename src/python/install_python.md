# 如何在 Ubuntu 中安装最新的 Python 版本

Python 是增长最快的主要通用编程语言。其原因有很多，例如其可读性和灵活性、易于学习和使用、可靠性和效率。

目前使用的 Python 有两个主要版本 – 2 和 3（Python 的现在和未来）；前者不会出现新的主要版本，而后者正在积极开发中，并且在过去几年中已经发布了许多稳定版本。 Python 3 的最新稳定版本是版本 3.11。

在较新的 Ubuntu 版本上，预安装了 Python 3.10 或 Python 3.8，而较旧的 Ubuntu 版本则不然。

在本文中，我们将解释如何使用 deadsnakes PPA 通过 apt 包管理器在所有 Ubuntu 版本上安装最新的 Python 3.11 版本。

要从所有主要 Linux 发行版中的源安装最新版本的 Python，请查看本指南：



## Install

要安装最新的 Python 3.11 版本，您可以使用“deadsnakes”团队 PPA，其中包含为 Ubuntu 打包的最新 Python 版本。

```sh
$ sudo add-apt-repository ppa:deadsnakes/ppa
$ sudo apt update
$ sudo apt install python3.11
```

如果您想在 Ubuntu 系统中安装特定的 Python 版本或多个版本的 Python，只需运行以下命令并输入所示的 Python 版本号即可。

```sh
$ sudo apt install python3.10
$ sudo apt install python3.9
$ sudo apt install python3.8
$ sudo apt install python3.7
$ sudo apt install python3.6
```

要查看系统上安装的所有 Python 二进制文件的列表，请运行以下 ls 命令。

```sh
$ ls -l /usr/bin/python*
```

列出 Python 二进制文件

```sh
lrwxrwxrwx 1 root root      10 Apr 22  2022 /usr/bin/python3 -> python3.10
-rwxr-xr-x 1 root root 5901416 Apr  2  2022 /usr/bin/python3.10
-rwxr-xr-x 1 root root 6705016 Oct 24 15:56 /usr/bin/python3.11
-rwxr-xr-x 1 root root     960 Dec 23  2020 /usr/bin/python3-futurize
-rwxr-xr-x 1 root root     964 Dec 23  2020 /usr/bin/python3-pasteurize
```

从上面截图的输出来看，测试系统默认的Python版本是3.10，您也可以使用以下命令检查Python版本。

```sh
$ python -V

Python 3.10.4
```

要使用 Python 11，请调用以下命令。

```sh
$ python3.11
```

访问 Python Shell

```sh
Python 3.11.0 (main, Oct 24 2022, 19:56:13) [GCC 11.2.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> print ("TecMint #1 Linux Blog");
TecMint #1 Linux Blog
>>> quit()
```

要退出 Python 解释器，请键入以下命令并按 Enter。

```sh
quit()
OR
exit()
```



## 设置默认版本

如果您在 Ubuntu 系统中安装了多个版本的 Python，并且只想将一个版本设置为默认版本，那么您需要执行一些额外的步骤，如图所示。

```sh
$ python3 --version
$ sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.10 1
$ sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 2
$ sudo update-alternatives --config python3
$ python3 --version
```

![](https://s2.loli.net/2023/08/16/4OxcFjytbDfGkNr.png)



就这样！在这篇短文中，我们解释了如何通过 apt 包管理器在 Ubuntu 中安装 Python 3.11。