# 实战|如何在Linux 系统上免费托管网站

Web 服务器可用于指代硬件和软件，或者两者一起工作。出于本指南的目的，我们将重点关注软件方面，并了解如何在 Linux 机器上托管网站。

Web 服务器是一种通过 HTTP/HTTPS 协议接收并响应客户端请求的软件程序。其主要目的是显示网站内容，这些内容通常采用文本、图像和视频的形式。

Web 服务器可以提供静态或动态内容。静态内容，顾名思义，是指几乎不会改变并且必然保持不变的内容。服务器按原样将内容发送回用户的浏览器。

动态内容是经常变化或不断更新的内容。为了提供动态内容，Web 服务器还必须与数据库服务器和服务器端脚本语言一起工作。

[本指南](https://www.tecmint.com/host-website-locally/ "Source")将演示如何设置 Apache Web 服务器以在 Linux 系统上免费托管网站。



## 依赖

要按照本指南进行操作，请确保您具备以下条件。

- 可以从您的 ISP 获取专用公共 IP 地址。
- Linux 盒子，可以是您首选操作系统变体的 Linux 服务器安装。在本指南中，我们将使用 Debian 11。

您还需要安装 LAMP 服务器，它是 Linux、Apache 和 MySQL（也可以是 MariaDB）的缩写。



## 如何在 Linux 服务器上托管网站

在本节中，我们将继续讨论 Web 服务器的主要组件。



### 什么是 Apache？

Apache 是一种流行的免费开源跨平台 Web 服务器，在 Apache License 2.0 下发布。它是使用最广泛的 Web 服务器之一，占据近 32.2% 的 Web 服务器市场份额。

要检查可用的 Apache 最新版本以及您的服务器上是否安装了该版本，请运行以下命令：

```sh
apt-cache policy apache2 (On Debian-based OS)
```

从输出中，您可以看到参数 Installed: (none) 表示尚未安装。您还可以获得有关 Debian / Ubuntu 存储库提供的最新版本的信息，在本例中为 2.4.52。

![](https://s2.loli.net/2023/07/05/vUkXjC4SsmrcpAG.png)



在现代 Red Hat 发行版上，您可以使用以下 dnf 命令检查 Apache 的可用性，如下所示。

```sh
dnf search httpd
```

![](https://s2.loli.net/2023/07/05/VQz5uxY2rMyPI4o.png)



从上面的输出中，您可以看到 Apache httpd 包可供下载。如果您的系统上未安装 Apache，请使用“apt”或“dnf”包管理器来安装 Apache，如图所示。

在基于 Debian 的系统上：

```sh
$ sudo apt install apache2 -y 	 
$ sudo systemctl start apache2	 
$ sudo systemctl enable apache2	 
$ sudo systemctl status apache2
```

![](https://s2.loli.net/2023/07/05/mpbCQJy9SviMEog.png)



在基于 Red Hat 的系统上：

```sh
# dnf install httpd -y 	 
# systemctl start httpd	 
# systemctl enable httpd	 
# systemctl status httpd
```

![](https://s2.loli.net/2023/07/05/yDnQhElAYWVvupP.png)



### 什么是 MariaDB？

MariaDB 是 MySQL 的一个分支，是最流行的开源关系数据库管理系统之一。如今，它比 MySQL 更受欢迎，因为它具有更快的复制和执行查询速度以及安全性和大量的存储引擎。

要在基于 Debian 的系统上安装 MariaDB：

```sh
$ sudo apt install mariadb-server mariadb-client -y	 
$ sudo systemctl start mariadb	 
$ sudo systemctl enable mariadb	 
$ sudo systemctl status mariadb	 
```

以下输出显示 MariaDB 已安装并按预期运行。

![](https://s2.loli.net/2023/07/05/3X5Bps6JH7cuniQ.png)



要在基于 RHEL 的系统上安装 MariaDB：

```sh
# dnf install mariadb-server -y	 
# systemctl start mariadb	 
# systemctl enable mariadb	 
# systemctl status mariadb	
```

![](https://s2.loli.net/2023/07/05/tqL6hsnXaWJFOe9.png)





### 什么是 PHP？

PHP 是 PHP 超文本预处理器的递归缩写，它是一种流行的通用脚本语言，主要用于 Web 开发。

要在基于 Debian 的系统上安装 PHP：

```sh
$ sudo apt update
$ sudo apt upgrade
$ sudo apt install  ca-certificates apt-transport-https software-properties-common
$ sudo add-apt-repository ppa:ondrej/php
$ sudo apt update
$ sudo apt install php8.0 libapache2-mod-php8.0 
```

要在基于 RHEL 的系统上安装 PHP，您需要首先启用 EPEL 存储库。

```sh
$ sudo dnf install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-9.noarch.rpm  [RHEL 9]
$ sudo dnf install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm  [RHEL 8]
```

接下来，启用 Remi 存储库，它在基于 RHEL 的系统上提供最新版本的 PHP。

```sh
$ sudo dnf install -y https://rpms.remirepo.net/enterprise/remi-release-9.rpm  [RHEL 8]
$ sudo dnf install -y https://rpms.remirepo.net/enterprise/remi-release-8.rpm  [RHEL 8]
```

在系统上启用 EPEL 和 Remi 存储库后，您可以如图所示安装 PHP。

```sh
# dnf module list php
# dnf module enable php:remi-8.0 -y 
# dnf install php php-cli php-common
```

安装所有组件后，您现在可以使用 WordPress CMS 构建网站，该软件使用户可以轻松开发和管理网站，而无需了解 HTML、CSS、PHP 和 Javascript 等网页设计语言。



## WordPress 建站

为了进行演示，我们将在 Debian 11 和 RHEL 9 系统上安装 WordPress，这将提供一个示例网站，可以根据您的喜好进一步定制。

本节假设您已经安装了 LAMP 堆栈。

### 1. 安装附加 PHP 模块

要继续，请安装 WordPress 所需的其他 PHP 模块，如图所示。

要在基于 Debian 的系统上安装 PHP 模块：

```sh
$ sudo apt install php libapache2-mod-php php-pear php-cgi php-common php-mbstring php-zip php-net-socket php-gd php-mysql php-bcmath
```

要在基于 RHEL 的系统上安装 PHP 模块：

```sh
# dnf install php-gd php-soap php-intl php-mysqlnd php-pdo php-bcmath php-curl php-zip php-xmlrpc wget
```



### 2. 为 WordPress 创建数据库

WordPress 用 PHP 编写，是一个数据驱动的、免费的开源内容管理系统。数据库是 WordPress 的重要组成部分。

该数据库用于存储所有博客文章、页面、类别、评论、主题、插件以及 WordPress 配置文件。

要为 WordPress 创建数据库，请登录 MariaDB 数据库服务器：

```sh
$ sudo mysql -u root -p
```

接下来，创建数据库，如图所示

```sh
CREATE DATABASE wordpress_db;
```

接下来，创建一个数据库用户并将数据库上的所有权限分配给该用户。

```sh
GRANT ALL PRIVILEGES ON wordpress_db.* to wordpress_user@localhost identified by 'P@ssword321';
```

然后最后重新加载授权表以保存所做的更改并退出数据库。

```sh
FLUSH PRIVILEGES;
QUIT;
```



### 3. 下载WordPress

数据库就位后，继续使用 wget 命令下载最新的 WordPress tarball 文件。

```sh
$ wget https://wordpress.org/latest.tar.gz
```

下载后，使用 tar 命令解压缩压缩文件。

```sh
$ tar -xvzf latest.tar.gz
```

该命令将文件的内容提取到名为 wordpress 的文件夹中。将文件夹移动或复制到 Apache Web 服务器的文档根目录中。

```sh
$ sudo mv wordpress/ /var/www/html/
```

接下来，分配以下权限和所有权。

```sh
$ sudo chmod 755 -R /var/www/html/wordpress/
$ sudo chown -R www-data:www-data /var/www/html/wordpress/
```



### 4. 为 WordPress 创建 Apache 虚拟主机

术语虚拟主机是指在单个服务器上托管多个网站的做法。如果您打算在一台服务器上托管多个网站，则需要为每个网站创建一个虚拟主机。

在这种情况下，您需要为 WordPress 网站创建虚拟主机，如下所示。

```sh
$ sudo nano /etc/apache2/sites-available/wordpress.conf  [On Debian]
# vi /etc/httpd/conf/httpd.conf [On RHEL]
```

粘贴以下代码行来定义虚拟主机。对于 ServerName 指令，提供服务器的 IP 地址或完全限定域名，它应指向专用公共 IP 地址。

```sh
<VirtualHost *:80>
     ServerAdmin admin@your_domain.com
     DocumentRoot /var/www/html/wordpress
     ServerName 192.168.0.100

     <Directory /var/www/html/wordpress>
          Options FollowSymlinks
          AllowOverride All
          Require all granted
     </Directory>

     ErrorLog ${APACHE_LOG_DIR}/your-domain.com_error.log
     CustomLog ${APACHE_LOG_DIR}/your-domain.com_access.log combined

</VirtualHost>
```

保存更改并退出文件。

要连接到数据库，需要进行一些额外的修改。因此，导航到 wordpress 文件夹。

```sh
$ cd /var/www/html/wordpress/
```

接下来，使用 wp-config-sample.php 文件的内容更新 wp-config.php 文件。

```sh
$ cp wp-config-sample.php wp-config.php
$ sudo nano wp-config.php
```

接下来，使用数据库详细信息更新数据库名称、数据库用户名和密码指令。

接下来，在基于 Debian 的系统上启用新的 WordPress 站点，如下所示。

```sh
$ sudo ln -s /etc/apache2/sites-available/wordpress.conf /etc/apache2/sites-enabled/wordpress.conf
$ sudo a2ensite wordpress
$ sudo a2enmod rewrite
$ sudo a2dissite 000-default
```

要使更改生效，请重新启动 Apache。

```sh
$ sudo systemctl restart apache2   [On Debian]
# systemctl restart httpd  [On RHEL]
```



### 5. 在浏览器上完成 WordPress 设置

要完成设置，请浏览 Web 服务器的 IP 地址，如下所示：

```sh
http://server-ip
```

您应该会看到 WordPress 欢迎页面，如图所示。选择您的首选语言，然后单击“继续”。

![](https://s2.loli.net/2023/07/05/CX2Y4tkHjxv6JgA.png)



接下来，填写站点详细信息。

![](https://s2.loli.net/2023/07/05/Sb6qGfvPHw9RgKx.png)



然后单击“安装 WordPress”以完成 WordPress 设置。

![](https://s2.loli.net/2023/07/05/BDneKyFYOLdw2PE.png)



如果一切顺利，您将收到安装成功的确认信息。要登录，请单击“登录”按钮。

![](https://s2.loli.net/2023/07/05/oWPChyHgFOs6Mfa.png)



如您所见，这将引导您进入 WordPress 仪表板。此时，您可以尝试使用各种主题来增强示例网站的外观。

![](https://s2.loli.net/2023/07/05/oeCkr84L2MgGVhB.png)





### 6. 使用端口转发访问 WordPress

由于您是通过家里的 Linux 系统或局域网 (LAN) 自托管 Web 服务器，因此下一步是让外部用户或 LAN（局域网）之外的用户可以访问它。这就是端口转发的用武之地。

端口转发，也称为端口映射，是一种允许外部设备通过 Internet 访问专用网络内的服务器或资源的技术。整个想法是从外部访问专用网络，否则这是不可能的，因为外部设备无法与内部 IP 地址通信。

在您的设置中，您需要转发 Web 服务器正在侦听的端口（在大多数情况下，对于 HTTP 流量是端口 80，对于 HTTPS 是端口 443）以及 Web 服务器的静态专用 IP 地址。

因此，登录您的路由器并前往端口转发部分。在我们的示例中，我们使用 DLink 路由器将 Web 服务器的端口（80 和 443）和私有 IP (192.168.0.100) 端口转发到 ISP 分配的专用 IP 公共 IP。

根据您的情况，指定 Web 服务器的端口和专用 IP 并保存更改。

![](https://s2.loli.net/2023/07/05/KcHzOIUiSmghGML.png)



要保存更改，您可能需要重新启动路由器。所以，继续做吧。

正确执行端口转发后，您现在可以通过公共 IP 地址访问网络外部的 Web 服务器。



## 总结

在本指南中，我们演示了如何在 Linux 机器上使用 Apache 自行托管 Web 服务器。欢迎您对本指南提供反馈。