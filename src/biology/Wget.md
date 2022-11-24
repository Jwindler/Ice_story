# wget



## 导读

本文将介绍`wget`的基本使用方法，和一些高级用法，比如**递归下载**等。对于经常在`FTP`网页下载数据的读者来说，可以说是必备的技能之一。



## 1. 介绍

[Wget](https://www.hostinger.com/tutorials/wget-command-examples/ "Wget") 是由 `GNU `项目创建的计算机工具。您可以使用它从各种 `Web `服务器检索内容和文件。这个名字是万维网和`get`这个词的组合。它支持通过` FTP、SFTP、HTTP 和 HTTPS `下载。

`Wget` 是用可移植的 C 语言创建的，可在任何` Unix` 系统上使用。也可以在` Mac OS X、Microsoft Windows、AmigaOS `和其他流行平台上实现。



## 2. 安装

- Ubuntu 18.04

```sh
sudo apt-get install wget
```



- CentOS 7

```sh
sudo yum install wget
```



## 3. Command



### 3.1. 单个文件

```sh
# 下载单个文件到当前文件夹
wget https://example.zip  # wget url
```



### 3.2. 多个文件

需要将多个件的`url`写到一个`txt`文件中，再利用`wget`下载

```sh
# example.txt
https://example1.zip

https://example2.zip

https://example3.zip
```

- 下载上面`example.txt`文件中的文件

```sh
wget -i example.txt
```



### 3.3. 命名

- 给下载文件重命名

```sh
wget -O example.zip https://ttt.zip  

# 利用-o 选项，重命名文件为`example.zip`
```



### 3.4. 指定目录

- 将下载文件保存到指定目录

```sh
wget -P documents/archives/ https://example.zip

# 利用-p 选项，将文件保存到`documents/archives/`目录下
```



### 3.5. 限制下载速度

```sh
wget --limit-rate=500k https://example.zip

# 将下载速度最高限制为 500k
```



### 3.6. 重试尝试次数

```sh
wget -tries=100 https://example.zip

# 设置重新连接次数为100
```



### 3.7. 后台下载

当下载的文件非常大的时候，可以将下载任务放置到后台

```sh
wget -b https://example.zip
```



### 3.8. FTP下载

```sh
wget --ftp-user=YOUR_USERNAME --ftp-password=YOUR_PASSWORD ftp://example.com/something.tar
```

- `--ftp-user`：FTP用户名
- `--ftp-password`：密码



### 3.8. 断点续连

当再下载途中，链接中断时，可以使用`-c`选项，进行重新连接,继续上次下载。

```sh
wget -c https://example.zip
```



### 3.9. 检索全站

- 下载整个网站的内容

```sh
wget --mirror --convert-links --page-requisites --no-parent -P documents/websites/ https://example.com
```

| 参数                 | 作用                             |
| -------------------- | -------------------------------- |
| **–mirror**          | 递归下载                         |
| **–convert-links**   | 所有链接都将转换为正确的脱机使用 |
| **–page-requisites** | 下载将包括CSS、JS和图像          |
| **–no-parent**       | 不检索父目录                     |
| **-P**               | 指定保存目录                     |



### 3.10. 查找断开链接

- 查找网页中无法下载的连接，并输出到文件中

```sh
wget -o wget-log -r -l 5 --spider http://example.com
```

| 参数        | 作用                         |
| ----------- | ---------------------------- |
| **-o**      | 将输出收集到文件中供以后使用 |
| **-l**      | 指定递归级别                 |
| **-r**      | 递归下载                     |
| **–spider** | 将wget设置为spider模式       |

- 利用下面命令，过滤出无法下载的文件

```sh
grep -B 2 '404' wget-log | grep "http" | cut -d " " -f 4 | sort -u

# wget-log是第一步的输出结果
```



### 3.11. 下载编号文件

如果文件名是按照数字编号时，可以同时下载。

```sh
wget http://example.com/images/{1..50}.jpg
```



## 4. 常用

```sh
nohup wget -c -r -np -L -P ./ ftp://download.big.ac.cn/gsa/CRA004538 > download.log 2>&1 &
# 一定要是ftp连接，不然容易变成爬网站
# -P 表示下载到哪个目录
# -r 表示递归下载，下载指定网页某一目录下（包括子目录）的所有文件
# -np 不要追溯到父目录
# -k 表示将下载的网页里的链接修改为本地链接.（下载整个站点后脱机浏览网页，最好加上这个参数
# -p 获得所有显示网页所需的元素，如图片等
# -c 断点续传
# -nd 递归下载时不创建一层一层的目录，把所有的文件下载到当前目录
# -o 将log日志指定保存到文件（新建一个文件）
# -a, –append-output=FILE 把记录追加到FILE文件中
# -A 指定要下载的文件样式列表，多个样式用逗号分隔
# -A zip 只下载指定文件类型（zip）
# -N 不要重新下载文件除非比本地文件新
# -O test.zip 下载并以不同的文件名保存
# -nc 不要覆盖存在的文件或使用.#前缀
# -m, –mirror 等价于 -r -N -l inf -nr
# -L 递归时不进入其它主机，如wget -c -r www.xxx.org 如果网站内有一个这样的链接： www.yyy.org，不加参数-L，就会像大火烧山一样，会递归下载www.yyy.org网站
# -i 后面跟一个文件，文件内指明要下载的URL，常用于多个url下载
# -nc 不要重复下载已存在的文件 --no-clobber
```

