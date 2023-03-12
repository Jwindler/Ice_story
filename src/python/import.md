# Python 包 install+import+管理用法总结



## 1. 简介

在Python中，包是一种组织代码的方式，它可以包含多个模块，并且可以将它们组合在一起以提供更高级别的功能。一个包是一个文件夹，里面包含一个名为__init__.py的空文件，这个文件告诉Python这个文件夹是一个包。包中可以包含子包，它们也是包的形式。



## 2. 安装

在Python中，我们可以通过pip来安装包。pip是Python的包管理器，可以用来安装和管理Python包。

### 2.1. 基本用法

以下是安装包的基本步骤：

1. 确认你已经安装了Python和pip

如果你还没有安装Python和pip，你需要先安装它们。你可以在Python官网上下载Python并安装，pip已经默认随Python一起安装了，所以你只需要在终端（或命令提示符）中输入以下命令来检查pip是否已安装：

```python
pip --version
```

如果你看到了pip的版本信息，那么pip已经安装好了。如果没有，你需要在终端中运行以下命令来安装pip：

```python
python get-pip.py
```

2. 使用pip来安装包

一旦你确认了pip已经安装好了，你就可以使用pip来安装包了。在终端中运行以下命令来安装一个包：

```python
pip install package_name
```

其中，`package_name`是你想要安装的包的名称。

如果你想要安装的是特定版本的包，可以使用以下命令：

```python
pip install package_name==version_number
```

其中，`version_number`是你想要安装的包的版本号。



### 2.2. 高阶方法

1. 安装特定版本的包

如果你想要安装特定版本的包，你可以使用以下命令：

```python
pip install package_name==version_number
```

其中，`version_number`是你想要安装的包的版本号。

2. 安装本地的包

如果你有一个本地的包文件（通常是以`.tar.gz`或`.whl`结尾的文件），你可以使用以下命令来安装它：

```python
pip install path/to/package_file
```

其中，`path/to/package_file`是你本地包文件的路径。

3. 升级已安装的包

如果你想要升级已安装的包，你可以使用以下命令：

```python
pip install --upgrade package_name
```

其中，`package_name`是你想要升级的包的名称。

4. 列出已安装的包

如果你想要列出已安装的包，你可以使用以下命令：

```python
pip list
```

这将列出你当前Python环境中已安装



## 3. 导包方法

在Python中，导入包可以使用`import`语句，也可以使用`from`语句导入包中的特定模块或函数。下面是一些常用的Python导入包的方法及其相应的代码实例。

### 3.1. import

使用`import`语句可以导入整个包或模块。一般格式如下：

```python
import package_name
```

其中，`package_name`为要导入的包名。

例如，导入`numpy`包，可以使用以下代码：

```python
import numpy
```

这样，就可以使用`numpy`包中的所有函数和模块。



### 3.2. from

使用`from`语句可以导入包中的特定模块或函数。一般格式如下：

```python
from package_name import module_name
```

其中，`package_name`为要导入的包名，`module_name`为要导入的模块名。

例如，导入`numpy`包中的`array`模块，可以使用以下代码：

```py
from numpy import array
```

这样，就可以直接使用`array`函数，而不需要加上`numpy.`前缀。



### 3.3. as

使用`as`关键字可以为导入的包或模块指定别名，使得代码更加简洁易懂。一般格式如下：

```python
import package_name as alias_name
```

或者

```python
from package_name import module_name as alias_name
```

例如，为`numpy`包指定别名`np`，可以使用以下代码：

```python
import numpy as np
```

这样，就可以使用`np`代替`numpy`。



### 3.4. 通配符

使用`*`通配符可以导入包中的所有模块和函数。一般格式如下：

```python
from package_name import *
```

例如，导入`numpy`包中的所有模块和函数，可以使用以下代码：

```python
from numpy import *
```

但是，使用`*`通配符会导致代码可读性下降，容易出现命名冲突等问题，因此不建议在生产环境中使用。



## 4. 包管里

在Python中，有一些高级的方法可以用于包的管理：

1. Virtualenv：虚拟环境是一种创建隔离的Python环境的方法，可以用于不同的项目之间，避免不同项目中的依赖冲突。
2. Pip：pip是Python的包管理器，可以用于安装、卸载和管理Python包。
3. Setuptools：Setuptools是一个工具集，它可以帮助你创建和分发Python包。
4. PyPI：PyPI是Python Package Index的缩写，是Python的软件仓库，其中包含了数以万计的Python包。
5. Anaconda：Anaconda是一个开源的Python发行版，包含了数以千计的科学计算和数据分析的Python包。

总之，包是Python中组织代码的一种方式，可以用于实现更高级别的功能。高阶方法可以帮助你更好地管理和使用Python包。

