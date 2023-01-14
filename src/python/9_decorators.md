# 9个Python 内置装饰器: 显著优化代码

装饰器是应用“Python 之禅”哲学的最佳 Python 特性。装饰器可以帮助您编写更少、更简单的代码来实现复杂的逻辑并在任何地方重用它。

更重要的是，有许多很棒的内置 Python 装饰器使我们的生活变得更加轻松，因为我们只需使用一行代码就可以为现有的函数或类添加复杂的功能。

让我们看看我精心挑选的 9 个装饰器，[本文](https://medium.com/techtofreedom/9-python-built-in-decorators-that-optimize-your-code-significantly-bc3f661e9017 "Source")将向您展示 Python 是多么优雅。



## 1. @lru_cache

使用缓存技巧加速 Python 函数的最简单方法是使用 @lru_cache 装饰器。

这个装饰器可以用来缓存一个函数的结果，这样后续调用相同参数的函数就不会再执行了。它对于计算量大或使用相同参数频繁调用的函数特别有用。

让我们看一个直观的例子：

```python
import time


def fibonacci(n):
    if n < 2:
        return n
    return fibonacci(n - 1) + fibonacci(n - 2)


start_time = time.perf_counter()
print(fibonacci(30))
end_time = time.perf_counter()
print(f"The execution time: {end_time - start_time:.8f} seconds")
# The execution time: 0.18129450 seconds
```

上面的程序使用 Python 函数计算第 N 个斐波那契数。计算fibonacci(30)的时候很耗时，很多前面的Fibonacci数在递归过程中会计算很多次。

现在，让我们使用@lru_cache 装饰器来加速它：

```py
from functools import lru_cache
import time


@lru_cache(maxsize=None)
def fibonacci(n):
    if n < 2:
        return n
    return fibonacci(n - 1) + fibonacci(n - 2)


start_time = time.perf_counter()
print(fibonacci(30))
end_time = time.perf_counter()
print(f"The execution time: {end_time - start_time:.8f} seconds")
# The execution time: 0.00002990 seconds
```

如上代码所示，使用@lru_cache装饰器后，我们可以在0.00002990秒内得到相同的结果，比之前的0.18129450秒快了很多。

@lru_cache 装饰器有一个 maxsize 参数，指定要存储在缓存中的最大结果数。当缓存已满并且需要存储新结果时，最近最少使用的结果将从缓存中逐出以为新结果腾出空间。这称为最近最少使用 (LRU) 策略。

默认情况下，maxsize 设置为 128。如果将其设置为 None，如我们的示例，LRU 功能将被禁用，并且缓存可以无限增长。



## 2. @total_ordering

`functools` 模块中的 `@total_ordering` 装饰器用于根据定义的方法为 Python 类生成缺少的比较方法。

```py
from functools import total_ordering


@total_ordering
class Student:
    def __init__(self, name, grade):
        self.name = name
        self.grade = grade

    def __eq__(self, other):
        return self.grade == other.grade

    def __lt__(self, other):
        return self.grade < other.grade


student1 = Student("Alice", 85)
student2 = Student("Bob", 75)
student3 = Student("Charlie", 85)

print(student1 < student2)  # False
print(student1 > student2)  # True
print(student1 == student3)  # True
print(student1 <= student3) # True
print(student3 >= student2) # True
```

如上面的代码所示，在 Student 类中没有定义 __ge__、__gt__ 和 __le__ 方法。但是，感谢@total_ordering 装饰器，我们在不同实例之间的比较结果都是正确的。

这个装饰器的好处是显而易见的：

- 它可以使您的代码更清晰并节省您的时间。因为你不需要写所有的比较方法。
- 一些旧类可能没有定义足够的比较方法。将 @total_ordering 装饰器添加到它以供进一步使用会更安全。



## 3. @contextmanager

Python 有一个上下文管理器机制来帮助你正确地管理资源。

```py
with open("test.txt",'w') as f:
    f.write("Yang is writing!")
```

如上面的代码所示，我们可以使用 with 语句打开一个文件，这样它会在写入后自动关闭。我们不需要显式调用 f.close() 函数来关闭文件。

有时，我们需要为一些特殊的需求定义一个自定义的上下文管理器。在这种情况下，@contextmanager 装饰器是我们的选择。

例如，下面的代码实现了一个简单的自定义上下文管理器，它可以在文件打开或关闭时打印相应的信息。

```py
from contextlib import contextmanager

@contextmanager
def file_manager(filename, mode):
    print("The file is opening...")
    file = open(filename,mode)
    yield file
    print("The file is closing...")
    file.close()

with file_manager('test.txt', 'w') as f:
    f.write('Yang is writing!')
# The file is opening...
# The file is closing...
```



## 4. @property

Getter 和 setter 是面向对象编程 (OOP) 中的重要概念。对于类的每个实例变量，getter 方法返回它的值，而 setter 方法设置或更新它的值。鉴于此，getter 和 setter 也分别称为访问器和修改器。它们用于保护您的数据不被直接和意外地访问或修改。不同的 OOP 语言有不同的机制来定义 getter 和 setter。在 Python 中，我们可以简单地使用 @property 装饰器。

```py
class Student:
    def __init__(self):
        self._score = 0

    @property
    def score(self):
        return self._score

    @score.setter
    def score(self, s):
        if 0 <= s <= 100:
            self._score = s
        else:
            raise ValueError('The score must be between 0 ~ 100!')

Yang = Student()

Yang.score=99
print(Yang.score)
# 99

Yang.score = 999
# ValueError: The score must be between 0 ~ 100!
```

如上例所示，score 变量不能设置为 999，这是一个无意义的数字。因为我们使用 @property 装饰器在 setter 函数中限制了它的可接受范围。

毫无疑问，添加这个setter可以成功避免意外的错误或结果。



## 5. @cached_property

Python 3.8 为 functool 模块引入了一个新的强大装饰器——@cached_property。它可以将一个类的方法转换为一个属性，该属性的值计算一次，然后在实例的生命周期内作为普通属性缓存。

```py
from functools import cached_property


class Circle:
    def __init__(self, radius):
        self.radius = radius

    @cached_property
    def area(self):
        return 3.14 * self.radius ** 2


circle = Circle(10)
print(circle.area)
# prints 314.0
print(circle.area)
# returns the cached result (314.0) directly
```

在上面的代码中，我们通过@cached_property 修饰了area 方法。所以没有对同一个不变实例的circle.area进行重复计算。



## 6. @classmethod

在 Python 类中，有 3 种可能的方法类型：

- 实例方法：绑定到实例的方法。他们可以访问和修改实例数据。在类的实例上调用实例方法，它可以通过 self 参数访问实例数据。
- 类方法：绑定到类的方法。他们不能修改实例数据。在类本身上调用类方法，它接收类作为第一个参数，通常命名为 cls。
- 静态方法：未绑定到实例或类的方法。

实例方法可以定义为普通的 Python 函数，只要它的第一个参数是 self。但是，要定义一个类方法，我们需要使用@classmethod 装饰器。

为了演示，以下示例定义了一个类方法，可用于通过直径获取 Circle 实例：

```py
class Circle:
    def __init__(self, radius):
        self.radius = radius

    @classmethod
    def from_diameter(cls, diameter):
        return cls(diameter / 2)

    @property
    def diameter(self):
        return self.radius * 2

    @diameter.setter
    def diameter(self, diameter):
        self.radius = diameter / 2


c = Circle.from_diameter(8)
print(c.radius)  # 4.0
print(c.diameter)  # 8.0
```



## 7. @staticmethod

如前所述，静态方法不绑定到实例或类。它们被包含在一个类中只是因为它们在逻辑上属于那个类。

静态方法通常用于执行一组相关任务（例如数学计算）的实用程序类。通过将相关函数组织到类中的静态方法中，我们的代码将变得更有条理，也更容易理解。

要定义一个静态方法，我们只需要使用@staticmethod 装饰器。让我们看一个例子：

```py
class Student:
    def __init__(self, first_name, last_name):
        self.first_name = first_name
        self.last_name = last_name
        self.nickname = None

    def set_nickname(self, name):
        self.nickname = name

    @staticmethod
    def suitable_age(age):
        return 6 <= age <= 70


print(Student.suitable_age(99)) # False
print(Student.suitable_age(27)) # True
print(Student('yang', 'zhou').suitable_age(27)) # True
```



## 8. @dataclass

@dataclass装饰器（Python 3.7引入）可以自动为一个类生成几个特殊的方法，如__init__、__repr__、__eq__、__lt__等。

因此，它可以为我们节省大量编写这些基本方法的时间。如果一个类主要用于存储数据，那么@dataclass 装饰器是最好的选择。

为了演示，下面的示例只定义了一个名为 Point 的类的两个数据字段。感谢 @dataclass 装饰器，它足以被使用：

```py
from dataclasses import dataclass

@dataclass
class Point:
    x: float
    y: float

point = Point(1.0, 2.0)
print(point)
# Point(x=1.0, y=2.0)
```



## 9. @atexit.register

来自 atexit 模块的 @register 装饰器可以让我们在 Python 解释器退出时执行一个函数。

这个装饰器对于执行最终任务非常有用，例如释放资源或只是说再见！

```py
import atexit

@atexit.register
def goodbye():
    print("Bye bye!")

print("Hello Yang!")
```

- 输出是：

```py
Hello Yang!
Bye bye!
```

如示例所示，由于使用了@register 装饰器，终端打印了“Bye bye!”即使我们没有显式调用再见函数。