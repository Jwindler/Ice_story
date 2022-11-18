# JSON load() and loads()

使用 `json.load()` 和 `loads()` 进行 `Python` ``JSON`` 解析



## 导读

[本文](https://pynative.com/python-json-load-and-loads-to-parse-json/ "Source")演示如何使用 `Python` 的 `json.load()` 和 `json.loads()` 方法从文件和字符串中读取 `JSON` 数据。使用 `json.load()` 和 `json.loads()` 方法，您可以将 `JSON` 格式的数据转换为 Python 类型，这个过程称为 `JSON` 解析。Python 内置模块 json 提供了以下两种解析 `JSON` 数据的方法。



要从 URL 或文件解析 `JSON`，请使用 `json.load()`。要解析包含 `JSON` 内容的字符串，请使用 `json.loads()`。

![`JSON` parsing](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221117220835801.png)



## 语法

我们可以使用 `load` 和 `loads()` 方法进行许多 `JSON` 解析操作。首先，了解它的语法和参数，然后我们逐一介绍它的用法。



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

`json.load()` 用于从文件中读取 `JSON` 文档，`json.loads()` 用于将 `JSON` 字符串文档转换为 Python 字典。

- `fp`

用于读取文本文件、二进制文件或 `JSON` 文件的文件指针。

- `object_hook`

是可选函数，将使用任何对象文字解码的结果调用。

- `object_pairs_hook`

是一个可选函数，将使用任何对象文字的结果调用，该对象文字是用有序的对列表解码的。

- `parse_float`

是可选参数，但如果指定，将使用要解码的每个 `JSON` 浮点数和整数的字符串调用。

- `parse_int`

如果指定，它将使用要解码的每个 `JSON` int 的字符串调用。默认情况下，这等同于 int(num_str)。



## json.load

`json.load()` 从文件中读取 `JSON` 数据并将其转换为字典。使用 `json.load()` 方法，我们可以从文本、`JSON` 或二进制文件中读取 `JSON` 数据。 `json.load()` 方法以 Python 字典的形式返回数据。然后我们使用这个字典来访问和操作我们的应用程序或系统中的数据。

`json.load()` 和 `json.loads()` 方法在解码时使用转换表，参考如下

- 解析转换表

| `JSON`          | Python |
| :------------ | :----- |
| object        | dict   |
| array         | list   |
| string        | str    |
| number (int)  | int    |
| number (real) | float  |
| true          | True   |
| false         | False  |
| null          | None   |



### 例子

现在，我正在读取硬盘上的“developer.json”文件。此文件包含以下 `JSON` 数据。

![developer.json](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221118153716740.png)



- 读取代码

```python
import json

print("Started Reading `JSON` file")
with open("developer.json", "r") as read_file:
    print("Converting `JSON` encoded data into Python dictionary")
    developer = json.load(read_file)

    print("Decoded `JSON` Data From File")
    for key, value in developer.items():
        print(key, ":", value)
    print("Done reading json file")
```



- 结果

```python
Started Reading `JSON` file
Converting `JSON` encoded data into Python dictionary

Decoded `JSON` Data From File
name : jane doe
salary : 9000
skills : ['Raspberry pi', 'Machine Learning', 'Web Development']
email : JaneDoe@pynative.com
projects : ['Python Data Mining', 'Python Data Science']

Done reading json file
```



### key

如果您想直接访问 `JSON` key 而不是从文件中迭代整个 `JSON`，使用以下代码

```python
import json

print("Started Reading `JSON` file")
with open("developer.json", "r") as read_file:
    print("Converting `JSON` encoded data into Python dictionary")
    developer = json.load(read_file)

    print("Decoding `JSON` Data From File")
    print("Printing `JSON` values using key")
    print(developer["name"])
    print(developer["salary"])
    print(developer["skills"])
    print(developer["email"])
    print("Done reading json file")
```



- 结果

```python
Started Reading `JSON` file
Converting `JSON` encoded data into Python dictionary

Decoding `JSON` Data From File
Printing `JSON` values using key

jane doe
9000
['Raspberry pi', 'Machine Learning', 'Web Development']
JaneDoe@pynative.com

Done reading json file
```



## json.loads

`json.loads()` 将 `JSON` 字符串转换为字典。有时我们会收到字符串格式的 `JSON` 数据。所以要在我们的应用程序中使用它，需要将 `JSON` 字符串转换为 Python 字典。使用 `json.loads()` 方法，我们可以将包含 `JSON` 文档的原生字符串、字节或字节数组实例反序列化为 Python 字典。



### 例子

- 参考如下

```python
import json

developerJsonString = """{
    "name": "jane doe",
    "salary": 9000,
    "skills": [
        "Raspberry pi",
        "Machine Learning",
        "Web Development"
    ],
    "email": "JaneDoe@pynative.com",
    "projects": [
        "Python Data Mining",
        "Python Data Science"
    ]
}
"""

print("Started converting `JSON` string document to Python dictionary")
developerDict = json.loads(developerJsonString)

print("Printing key and value")
print(developerDict["name"])
print(developerDict["salary"])
print(developerDict["skills"])
print(developerDict["email"])
print(developerDict["projects"])

print("Done converting `JSON` string document to a dictionary")
```



- 结果

```python
Started converting `JSON` string document to Python dictionary

Printing key and value
jane doe
9000
['Raspberry pi', 'Machine Learning', 'Web Development']
JaneDoe@pynative.com
['Python Data Mining', 'Python Data Science']

Done converting `JSON` string document to a dictionary
```



## 嵌套

- 解析和检索嵌套的 `JSON` 键值。

假设您有一个如下所示的 `JSON` 数据：

```python
developerInfo = """{
    "id": 23,
    "name": "jane doe",
    "salary": 9000,
    "email": "JaneDoe@pynative.com",
    "experience": {"python":5, "data Science":2},
    "projectinfo": [{"id":100, "name":"Data Mining"}]
}
"""
```



- 嵌套解析参考如下

```python
import json

print("Started reading nested `JSON` array")
developerDict = json.loads(developerInfo)

print("Project name: ", developerDict["projectinfo"][0]["name"])
print("Experience: ", developerDict["experience"]["python"])

print("Done reading nested `JSON` Array")
```



- 结果

```python
Started reading nested `JSON` array
Project name:  Data Mining
Experience:  5
Done reading nested `JSON` Array
```



## 有序字典

- 将 `JSON` 解析为 OrderedDict

正如我们上面讨论的那样，`json.load()` 方法的 object_pairs_hook 参数是一个可选函数，它将使用任何对象文字的结果调用，并使用有序的对列表进行解码。



- 参考如下

```python
import json
from collections import OrderedDict

print("Ordering keys")
OrderedData = json.loads('{"John":1, "Emma": 2, "Ault": 3, "Brian": 4}', object_pairs_hook=OrderedDict)
print("Type: ", type((OrderedData)))
print(OrderedData)
```



- 结果

```python
Ordering keys
Type:  <class 'collections.OrderedDict'>
OrderedDict([('John', 1), ('Emma', 2), ('Ault', 3), ('Brian', 4)])
```



## 类型

假设 `JSON` 文档包含许多浮点值，并且您希望将所有浮点值四舍五入到两位小数。在这种情况下，我们需要定义一个自定义函数来执行您想要的任何舍入。我们可以将这样的函数传递给 parse_float kwarg。当然 parse_int kwarg 也是如此。

- 参考如下

```python
import json

def roundFloats(salary):
    return round(float(salary), 2)

def salartToDeduct(leaveDays):
    salaryPerDay = 465
    return int(leaveDays) * salaryPerDay

print("Load float and int values from `JSON` and manipulate it")
print("Started Reading `JSON` file")
with open("developerDetails.json", "r") as read_file:
    developer = json.load(read_file, parse_float=roundFloats,
                          parse_int=salartToDeduct)
    # after parse_float
    print("Salary: ", developer["salary"])

    # after parse_int
    print("Salary to deduct: ", developer["leavedays"])
    print("Done reading a `JSON` file")
```



- 结果

```python
Load float and int values from `JSON` and manipulate it
Started Reading `JSON` file
Salary:  9250.542
<class 'float'>
Salary to deduct:  3
Done reading a `JSON` file
```

