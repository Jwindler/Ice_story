# 如何在 Linux 中设置 SSH 无密码登录

SSH（Secure SHELL）是一种开源且可信的网络协议，用于登录远程服务器以执行命令和程序。

它还用于使用安全复制 (SCP) 命令和 rsync 命令通过网络将文件从一台计算机传输到另一台计算机。

在[本文](https://www.tecmint.com/ssh-passwordless-login-using-ssh-keygen-in-5-easy-steps/ "Source")中，我们将向您展示如何在基于 RHEL 的 Linux 发行版（例如 CentOS、Fedora、Rocky Linux 和 AlmaLinux）以及基于 Debian 的发行版（例如 Ubuntu 和 Mint）上设置无密码登录，使用 ssh 密钥连接到远程Linux服务器无需输入密码。

使用带有 SSH 密钥的无密码登录将增加两个 Linux 服务器之间的信任，以便轻松同步或传输文件。



## 我的设置环境

```sh
SSH Client : 192.168.0.12 ( Fedora 36 )
SSH Remote Host : 192.168.0.11 ( CentOS 8 )
```

如果您正在处理多个 Linux 远程服务器，那么 SSH 无密码登录是自动化任务的最佳方法之一，例如使用脚本自动备份、使用 SCP 命令同步文件以及远程命令执行。

在本例中，我们将设置 SSH 无密码自动登录，从服务器 192.168.0.12 以用户 howtoing 登录到 192.168.0.11 以用户 sheena 登录。





## 1. 创建身份验证 SSH-Keygen 密钥

首先使用用户howtoing登录服务器192.168.0.12，并使用以下命令生成一对公钥。

```sh
$ ssh-keygen -t rsa
```

![](https://s2.loli.net/2023/08/16/vkh8CtNEd2IoPWG.png)



## 2. 上传 SSH 密钥

从服务器 192.168.0.12 使用 SSH，并在服务器 192.168.0.11 的 sheena 的 .ssh 目录下上传新生成的公钥（id_rsa.pub），文件名为authorized_keys。

```sh
$ ssh-copy-id sheena@192.168.0.11
```

确保对远程服务器上的 ~/.ssh 目录和 ~/.ssh/authorized_keys 文件设置正确的权限。

```sh
$ ssh sheena@192.168.0.11 "chmod 700 ~/.ssh && chmod 600 ~/.ssh/authorized_keys"
```



## 3. 禁用密码验证（可选）

为了提高安全性，您可以在远程服务器上禁用密码身份验证，仅允许 SSH 密钥身份验证。为此，请打开远程服务器上的 SSH 服务器配置文件：

```sh
$ sudo nano /etc/ssh/sshd_config
OR
$ sudo vi /etc/ssh/sshd_config
```

找到包含PasswordAuthentication 的行并将其设置为no。

```sh
PasswordAuthentication no
```

保存文件并重新启动 SSH 服务。

```sh
$ sudo systemctl restart sshd
```



## 4. 测试 SSH 无密码登录

从现在开始，您可以以 sheena 用户身份从服务器 192.168.0.12 以 howtoing 用户身份登录 192.168.0.11，无需密码。

```sh
$ ssh sheena@192.168.0.11
```

![](https://s2.loli.net/2023/08/16/tNw8gx2f79y3K1u.png)



在本文中，您学习了如何使用 ssh 密钥设置 SSH 无密码登录。我希望这个过程很简单。如果您有任何疑问，请在下面的评论部分发表。