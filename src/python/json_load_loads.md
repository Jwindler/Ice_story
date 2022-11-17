# JSON load() and loads()

使用 json.load() 和 loads() 进行 Python JSON 解析

https://pynative.com/python-json-load-and-loads-to-parse-json/#:~:text=Parameter%20used%3A-,The%20json.,that%20contains%20a%20JSON%20document.

本文演示如何使用 Python 的 json.load() 和 json.loads() 方法从文件和字符串中读取 JSON 数据。使用 json.load() 和 json.loads() 方法，您可以将 JSON 格式的数据转换为 Python 类型，这个过程称为 JSON 解析。Python 内置模块 json 提供了以下两种解析 JSON 数据的方法。



要从 URL 或文件解析 JSON，请使用 json.load()。要解析包含 JSON 内容的字符串，请使用 json.loads()。

![JSON parsing](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221117220835801.png)



## 语法

我们可以使用 load 和 loads() 方法进行许多 JSON 解析操作。首先，了解它的语法和参数，然后我们逐一介绍它的用法。



### load()

```python
json.load(fp, *, cls=None, object_hook=None, parse_float=None, parse_int=None, parse_constant=None, object_pairs_hook=None, **kw)

```



### loads()

```python
json.loads(s, *, cls=None, object_hook=None, parse_float=None, parse_int=None, parse_constant=None, object_pairs_hook=None, **kw)

```



### 参数

所有参数在两种方法中具有相同的含义。

json.load() 用于从文件中读取 JSON 文档，json.loads() 用于将 JSON 字符串文档转换为 Python 字典。

- fp

用于读取文本文件、二进制文件或 JSON 文件的文件指针。

- object_hook

是可选函数，将使用任何对象文字解码的结果调用。

- object_pairs_hook

是一个可选函数，将使用任何对象文字的结果调用，该对象文字是用有序的对列表解码的。

- `parse_float`

是可选参数，但如果指定，将使用要解码的每个 JSON 浮点数和整数的字符串调用。

- `parse_int`

如果指定，它将使用要解码的每个 JSON int 的字符串调用。默认情况下，这等同于 int(num_str)。