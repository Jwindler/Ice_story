# Python 闭包使用详解(内含实例)



## 1. 简介

闭包（closure）是一个函数对象，它与它的环境变量（包括自由变量）的引用组合而成的实体。闭包可以保留函数定义时所在的环境变量，即使这些变量在定义时不在该函数的作用域内也可以使用。闭包可以让函数像对象一样被传递、赋值、作为参数传递，甚至在函数内部定义函数。



## 2. 使用场合

闭包可以在以下场合使用：

1. 当需要一个函数在不同的环境中具有不同的行为时；
2. 当需要一个函数在其定义之后被动态修改其行为时；
3. 当需要在一个函数内部定义其他函数时；
4. 当需要将一个函数作为参数传递给另一个函数时。



## 3. 代码实例

```python
def outer_function(x):
    def inner_function(y):
        return x + y
    return inner_function

closure_example = outer_function(10)
print(closure_example(5)) # 输出 15
```

在这个示例中，`outer_function` 定义了一个内部函数 `inner_function`，它使用 `outer_function` 的参数 `x`。`outer_function` 返回 `inner_function`，并且 `inner_function` 被分配给变量 `closure_example`。

现在， `closure_example` 变成了一个闭包，因为它保留了 `outer_function` 的参数 `x`，并且在后续调用 `closure_example` 时可以使用。在上面的示例中，我们首先调用 `outer_function(10)` 并将其返回值分配给 `closure_example`。然后，我们调用 `closure_example(5)`，它将使用先前的 `x` 值 `10`，返回 `15`。



以下是另一个闭包示例，展示了如何在函数内部定义其他函数：

```python
def outer_function():
    def inner_function():
        print("Hello from inner function!")
    return inner_function

closure_example = outer_function()
closure_example() # 输出 "Hello from inner function!"

```

在这个示例中，`outer_function` 定义了一个内部函数 `inner_function`。然后，`outer_function` 返回 `inner_function`，并且 `inner_function` 被分配给变量 `closure_example`。最后，我们调用 `closure_example` 函数，它输出了 "Hello from inner function!"。

以上是闭包的简单介绍和两个示例，闭包可以在函数式编程中非常有用，因为它们允许我们创建高度抽象的功能。



## 4. 常用场景

### 4.1 计数器

闭包可以用来实现计数器，例如：

```python
def counter():
    count = 0
    def inner():
        nonlocal count
        count += 1
        return count
    return inner

c = counter()
print(c())  # 输出 1
print(c())  # 输出 2
print(c())  # 输出 3

```

在这个例子中，`counter` 函数返回了一个闭包 `inner`，它会在每次调用时返回计数器的值并将其增加。通过调用 `counter()` 函数，我们创建了一个闭包 `c`，然后每次调用 `c()` 都会返回一个新的计数器值。



### 4.2. 保护私有变量

闭包可以用来保护私有变量，例如：

```python
def create_savings_account(initial_balance):
    balance = initial_balance
    def deposit(amount):
        nonlocal balance
        balance += amount
    def withdraw(amount):
        nonlocal balance
        if amount > balance:
            print("Insufficient funds")
        else:
            balance -= amount
    def get_balance():
        return balance
    return deposit, withdraw, get_balance

deposit, withdraw, get_balance = create_savings_account(1000)
print(get_balance())  # 输出 1000
deposit(500)
print(get_balance())  # 输出 1500
withdraw(2000)        # 输出 "Insufficient funds"
withdraw(500)
print(get_balance())  # 输出 1000

```

在这个例子中，`create_savings_account` 函数返回了三个闭包 `deposit`、`withdraw` 和 `get_balance`，它们可以对私有变量 `balance` 进行操作并返回其值。通过调用 `create_savings_account` 函数，我们创建了一个带有初始余额的账户，并可以通过 `deposit` 和 `withdraw` 函数进行存款和取款，并可以通过 `get_balance` 函数查询余额。



### 4.3. 缓存函数结果

闭包可以用来缓存函数的结果，以提高函数的性能，例如：

```python
def memoize(func):
    cache = {}
    def inner(n):
        if n not in cache:
            cache[n] = func(n)
        return cache[n]
    return inner

@memoize
def fibonacci(n):
    if n <= 1:
        return n
    return fibonacci(n-1) + fibonacci(n-2)

print(fibonacci(10))  # 输出 55

```

在这个例子中，`memoize` 函数返回了一个闭包 `inner`，它会缓存 `func` 函数的结果并返回。通过将 `@memoize` 装饰器应用到 `fibonacci` 函数上，我们可以对其进行缓存，以避免重复计算。