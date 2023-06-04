# 如何修复 SSH Client_loop: send disconnect: Broken pipe Error

SSH 是 Secure Shell 的缩写，是一种远程网络协议，用于通过 TCP/IP 网络安全地连接到远程设备，例如服务器和网络设备。

它是一种加密网络协议，可提供强大的加密技术和散列法来保护网络上两个设备之间的通信。

SSH使用两种主要的认证方式：密码认证和公钥认证。使用密码验证时，用户提供远程主机的 IP 地址或 FQDN（完全限定域名）和密码进行验证。

公钥认证使用SSH 密钥对进行认证，SSH 密钥对由两个SSH 密钥组成：私钥和公钥。

私钥驻留在用户的机器上，应始终保密和安全。公钥保存在用户连接的远程主机上。在身份验证期间，比较两个密钥的身份并授予访问权限。

通过 SSH 连接到远程系统时，您可能会遇到错误 Client_loop: send disconnect: Broken pipe。

![](https://s2.loli.net/2023/06/04/7wy9fpRHgWL2UvA.png)



在[本教程](https://www.tecmint.com/client_loop-send-disconnect-broken-pipe/ "Source")中，我们将了解为什么会发生这种情况并解决错误。



## Client_loop: send disconnect: Broken pipe 错误

该错误只是一条断开连接消息，通知您已超过 SSH 连接超时。

这是一个不活动的时期，在此期间没有从客户端执行或发出任何 Linux 命令。发生这种情况时，SSH 会话将终止，从而有效地断开您与远程服务器的连接。

大多数用户通常会按“ENTER”或键盘上的某个键，以避免空闲 SSH 会话导致与主机断开连接。然而，这可能是乏味且浪费时间的。

值得庆幸的是，SSH 默认配置设置提供了一些参数，您可以配置这些参数以使 SSH 连接保持活动状态更长的时间。



## 修复 Client_loop: send disconnect: Broken pipe 错误

要解决此问题，您需要增加客户端上的 SSH 连接超时时间。为此，请修改通常位于 /etc/ssh/sshd_config 的默认 SSH 配置文件。

```sh
sudo vi /etc/ssh/sshd_config
```

请务必找到这两个参数：ClientAliveInterval 和 ClientAliveCountMax。让我们看看他们做了什么。

- ClientAliveInterval——这是一段不活动的时间，在此之后 SSH 服务器向连接到它的远程客户端发送一条活动消息。

- ClientAliveCountMax – 这是服务器尝试将活动消息从服务器发送到客户端的次数。

我们将这两个值设置如下：

```sh
ClientAliveInterval	300
ClientAliveCountMax	3
```

![](https://s2.loli.net/2023/06/04/iFLbeBj8zKoSWpI.png)



这意味着在客户端不活动的前 300 秒（5 分钟）之后，服务器将向客户端发送一条活动消息以保持 SSH 会话处于活动状态。

如果在接下来的 300 秒内（在 600 秒标记处）没有从客户端收到任何数据或响应，服务器将再次发送另一条活动消息。最后，在客户端不活动 900 秒后，SSH 连接将终止或断开。

请务必保存对文件所做的更改，然后退出。然后重新启动 SSH 守护程序。

```sh
sudo systemctl restart sshd
```

或者，您可以通过以秒（300 秒）为单位指定 ServerAliveInterval 参数来连接到您的远程客户端 Linux 系统，这意味着您的 SSH 会话处于活动状态最多 5 分钟。

```sh
ssh -o ServerAliveInterval=300 username@server_ip_address
```

![](https://s2.loli.net/2023/06/04/jwZ4OWx3SAsL78r.png)



在本教程中，我们演示了如何解决 Client_loop: send disconnect: Broken pipe 错误。如您所见，您只需在 SSH 配置文件中执行一些调整。