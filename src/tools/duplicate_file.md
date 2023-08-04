# 实践|Linux 中查找和删除重复文件

如果您习惯使用下载管理器从互联网上下载各种内容，那么组织您的主目录甚至系统可能会特别困难。

通常，您可能会发现您下载了相同的 mp3、pdf 和 epub（以及各种其他文件扩展名）并将其复制到不同的目录。这可能会导致您的目录中充满各种无用的重复内容。

在本教程中，您将学习如何使用 rdfind、fdupes 和 rmlint 命令行工具以及使用名为 DupeGuru 和 FSlint 的 GUI 工具在 Linux 中查找和删除重复文件。

请注意 - 始终小心您在系统上删除的内容，因为这可能会导致不必要的数据丢失。如果您使用新工具，请首先在测试目录中尝试，在该目录中删除文件不会出现问题。



## Rdfind – 在 Linux 中查找重复文件

Rdfind 来自冗余数据查找，它是一个免费的命令行工具，用于跨多个目录或多个目录内查找重复文件。它递归地扫描目录并识别具有相同内容的文件，允许您采取适当的操作，例如删除或移动重复项。

Rdfind 使用一种算法对文件进行分类，并检测哪些重复项是原始文件，并将其余的视为重复项。

排名规则为：

- 如果在扫描输入参数时发现 A 早于 B，则 A 的排名较高。
- 如果 A 的发现深度低于 B，则 A 的排名较高。
- 如果 A 早于 B 被发现，则 A 的排名较高。

最后一条规则特别适用于在同一目录中找到两个文件时。



### Install

要在 Linux 中安装 rdfind，请根据您的 Linux 发行版使用以下命令。

```sh
$ sudo apt install rdfind         [On Debian, Ubuntu and Mint]
$ sudo yum install rdfind         [On RHEL/CentOS/Fedora and Rocky/AlmaLinux]
$ sudo emerge -a sys-apps/rdfind  [On Gentoo Linux]
$ sudo apk add rdfind             [On Alpine Linux]
$ sudo pacman -S rdfind           [On Arch Linux]
$ sudo zypper install rdfind      [On OpenSUSE]    
```

要在目录上运行 rdfind，只需键入 rdfind 和目标目录。

```sh
$ rdfind /home/user
```

![](https://s2.loli.net/2023/08/04/a4CRt6yqu85bMjr.png)



如您所见，rdfind 会将结果保存在名为 results.txt 的文件中，该文件位于运行程序的同一目录中。该文件包含 rdfind 找到的所有重复文件。如果需要，您可以查看该文件并手动删除重复的文件。

您可以做的另一件事是使用 -dryrun 选项，该选项将提供重复项列表，而无需执行任何操作：

```sh
$ rdfind -dryrun true /home/user
```

当您找到重复项时，您可以选择用硬链接替换它们。

```sh
$ rdfind -makehardlinks true /home/user
```

如果您想删除重复项，您可以运行。

```sh
$ rdfind -deleteduplicates true /home/user
```

要检查 rdfind 的其他有用选项，您可以使用 rdfind 手册。

```sh
$ man rdfind 
```



## Fdupes – 扫描 Linux 中的重复文件

Fdupes 是另一个命令行程序，可让您识别系统上的重复文件。它递归地搜索目录，比较文件大小和内容以识别重复项。

它使用以下方法来确定重复文件：

- 比较部分 md5sum 签名
- 比较完整的 md5sum 签名
- 逐字节比较验证

就像 rdfind 一样，它有类似的选项：

- 递归搜索
- 排除空文件
- 显示重复文件的大小
- 立即删除重复项
- 排除具有不同所有者的文件



### Install

要在 Linux 中安装 fdupes，请根据您的 Linux 发行版使用以下命令。

```sh
$ sudo apt install fdupes         [On Debian, Ubuntu and Mint]
$ sudo yum install fdupes         [On RHEL/CentOS/Fedora and Rocky/AlmaLinux]
$ sudo emerge -a sys-apps/fdupes  [On Gentoo Linux]
$ sudo apk add fdupes             [On Alpine Linux]
$ sudo pacman -S fdupes           [On Arch Linux]
$ sudo zypper install fdupes      [On OpenSUSE]  
```

Fdupes 语法与 rdfind 类似。只需键入命令，然后键入您要扫描的目录即可。

```sh
$ fdupes <dir>
```

要递归搜索文件，您必须指定 -r 选项，如下所示。

```sh
$ fdupes -r <dir>
```

您还可以指定多个目录并指定要递归搜索的目录。

```sh
$ fdupes <dir1> -r <dir2>
```

要让 fdupes 计算重复文件的大小，请使用 -S 选项。

```sh
$ fdupes -S <dir>
```

要收集有关找到的文件的汇总信息，请使用 -m 选项。

```sh
$ fdupes -m <dir>
```

![](https://s2.loli.net/2023/08/04/Ht93b4xe8yRrawu.png)



最后，如果您想删除所有重复项，请使用 -d 选项，如下所示。

```sh
$ fdupes -d <dir>
```

Fdupes 将询问要删除哪个找到的文件。您需要输入文件编号：

![](https://s2.loli.net/2023/08/04/Tg7mF2Rx19narMt.png)



绝对不推荐的解决方案是使用 -N 选项，这将导致仅保留第一个文件。

```sh
$ fdupes -dN <dir>
```

要获取与 fdupes 一起使用的可用选项列表，请通过运行查看帮助页面。

```sh
$ fdupes -help
```



## Rmlint – 删除重复文件

Rmlint 是一个命令行工具，用于在 Linux 系统中查找和删除重复的和类似 lint 的文件。它有助于识别具有相同内容的文件，以及各种形式的冗余或 lint，例如空文件、损坏的符号链接和孤立文件。



### Install

要在 Linux 中安装 Rmlint，请根据您的 Linux 发行版使用以下命令。

```sh
$ sudo apt install rmlint         [On Debian, Ubuntu and Mint]
$ sudo yum install rmlint         [On RHEL/CentOS/Fedora and Rocky/AlmaLinux]
$ sudo emerge -a sys-apps/rmlint  [On Gentoo Linux]
$ sudo apk add rmlint             [On Alpine Linux]
$ sudo pacman -S rmlint           [On Arch Linux]
$ sudo zypper install rmlint      [On OpenSUSE]    
```

![](https://s2.loli.net/2023/08/04/KhRBfn6j9xdPS82.png)



## dupeGuru – 在 Linux 中查找重复文件

dupeGuru 是一个开源、跨平台的工具，可用于查找 Linux 系统中的重复文件。该工具可以扫描一个或多个文件夹中的文件名或内容。它还允许您找到与您正在搜索的文件相似的文件名。

dupeGuru 有适用于 Windows、Mac 和 Linux 平台的不同版本。其快速模糊匹配算法功能可帮助您在一分钟内找到重复文件。它是可定制的，您可以提取所需的精确重复文件，并从系统中擦除不需要的文件。



### Install

要在 Linux 中安装 dupeGuru，请根据您的 Linux 发行版使用以下命令。

```sh
$ sudo apt install dupeguru         [On Debian, Ubuntu and Mint]
$ sudo yum install dupeguru         [On RHEL/CentOS/Fedora and Rocky/AlmaLinux]
$ sudo emerge -a sys-apps/dupeguru  [On Gentoo Linux]
$ sudo apk add dupeguru             [On Alpine Linux]
$ sudo pacman -S dupeguru           [On Arch Linux]
$ sudo zypper install dupeguru      [On OpenSUSE]    
```

![](https://s2.loli.net/2023/08/04/SpBDilhd9rTXNf8.png)



## FSlint – 适用于 Linux 的重复文件查找器

FSlint 是一个免费实用程序，用于查找和清理文件系统上各种形式的 lint。它还报告重复文件、空目录、临时文件、重复/冲突（二进制）名称、错误的符号链接等等。它具有命令行和 GUI 模式。

然而，值得注意的是，截至 2022 年 9 月我所知，FSlint 的最后一次更新是在 2013 年，可能不会得到积极维护或与较新的 Linux 发行版兼容。



### Install

要在 Linux 中安装 FSlint，请根据您的 Linux 发行版使用以下命令。

```sh
$ sudo apt install fslint         [On Debian, Ubuntu and Mint]
$ sudo yum install fslint         [On RHEL/CentOS/Fedora and Rocky/AlmaLinux]
$ sudo emerge -a sys-apps/fslint  [On Gentoo Linux]
$ sudo apk add fslint             [On Alpine Linux]
$ sudo pacman -S fslint           [On Arch Linux]
$ sudo zypper install fslint      [On OpenSUSE] 
```

![](https://s2.loli.net/2023/08/04/kDP8QniG7Vh2Ltx.png)



##  总结

这些是在 Linux 系统上查找重复文件的非常有用的工具，但删除此类文件时应该非常小心。

如果您不确定是否需要某个文件，最好在删除该文件之前创建该文件的备份并记住其目录。如果您有任何问题或意见，请在下面的评论部分提交。