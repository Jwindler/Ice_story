# 如何在 Linux 中使用 Bash For 循环

在编程语言中，循环是必不可少的组件，当您想要一遍又一遍地重复代码直到满足指定条件时使用。

在 Bash 脚本中，循环扮演着几乎相同的角色，并用于自动执行重复性任务，就像在编程语言中一样。

在 Bash 脚本中，有 3 种类型的循环：for 循环、while 循环和 until 循环。这三个用于迭代值列表并执行一组给定的命令。

在[本指南](https://www.tecmint.com/bash-for-loop-linux/ "Source")中，我们将重点介绍 Linux 中的 Bash For 循环。



## 循环语法

如前所述，for 循环遍历一系列值并执行一组 Linux 命令。

For 循环采用以下语法：

```bash
for variable_name in value1 value2 value3  .. n
do
    command1
    command2
    commandn
done
```

现在让我们检查 bash for 循环的几个示例用法。



## 循环示例

在最简单的形式中，for 循环采用以下基本格式。在此示例中，变量 n 遍历一组用花括号括起来的数值，并将它们的值打印到标准输出。

```bash
for n in {1 2 3 4 5 6 7};
do
   echo $n
done
```

![](https://s2.loli.net/2023/06/15/IuJORUBYLg8f1iW.png)





## 带有范围的循环

在前面的示例中，我们明确列出了要由 for 循环迭代的值，效果很好。但是，您只能想象如果您要迭代（例如，一百个值），那将是一项多么繁琐和耗时的任务。这将迫使您键入从 1 到 100 的所有值。

要解决此问题，请指定一个范围。为此，请指定以两个句点分隔的开始和停止编号。

在此示例中，1 是第一个值，而 7 是范围中的最后一个值。

```bash
#!/bin/bash

for n in {1..7};
do
   echo $n
done
```

执行 shell 脚本后，将列出范围内的所有值，类似于我们在简单循环中的情况。

![](https://s2.loli.net/2023/06/15/Fc9q5e7oRZmUKxM.png)



此外，我们可以在范围的末尾包含一个值，该值将导致 for 循环以增量步骤迭代这些值。

以下 bash 脚本打印 1 到 7 之间的值，从第一个值开始在这些值之间增加 2 个步长。

```bash
#!/bin/bash

for n in {1..7..2};
do
   echo $n
done
```

![](https://s2.loli.net/2023/06/15/GWszubUvIoS7Btw.png)



从上面的示例中，您可以看到循环将花括号内的值递增了 2 个值。



## 数组循环

您还可以使用 for 循环轻松地遍历数组中定义的值。在以下示例中，for 循环遍历 fruits 数组中的所有值并将它们打印到标准输出。

```bash
#!/bin/bash

fruits=("blueberry" "peach" "mango" "pineapple" "papaya") 

for n in ${fruits[@]}; 
do
    echo $n
done
```

![](https://s2.loli.net/2023/06/15/jHeX2TEQJva6UxK.png)



@ 运算符访问或定位所有元素。这使得一个一个地遍历所有元素成为可能。

此外，您可以通过指定其在数组中的位置来访问单个元素。

例如，要访问“mango”元素，请将 @ 运算符替换为元素在数组中的位置（第一个元素从 0 开始，因此在这种情况下，“mango”将用 2 表示）。

这就是 for 循环的样子。

```bash
#!/bin/bash

fruits=("blueberry" "peach" "mango" "pineapple" "papaya") 

for n in ${fruits[2]}; 
do
    echo $n
done
```

![](https://s2.loli.net/2023/06/15/6zAljLUJhs3cx7q.png)



## C 风格的循环

您可以在循环内使用变量来迭代一系列元素。这就是 C 风格的 for 循环的用武之地。以下示例说明了 C 风格的 for 循环，它打印出从 1 到 7 的数值列表。

```bash
#!/bin/bash

n=7
for (( n=1 ; n<=$n ; n++ )); 
do
    echo $n
done
```

![](https://s2.loli.net/2023/06/15/2w8FgOhobRrNAzH.png)



 ##  C 风格的带有条件语句的循环

您可以在 C 风格的 for 循环中包含条件语句。在下面的示例中，我们包含了一个 if-else 语句，用于检查并打印出 1 到 7 之间的偶数和奇数。

```bash
#!/bin/bash

for (( n=1; n<=7; n++ ))
do  
    # Check if the number is even or not
    if (( $n%2==0 ))
    then
        echo "$n is even"
    else
        echo "$n is odd"
    fi  
done
```

![](https://s2.loli.net/2023/06/15/cEYP78CAq2wpthg.png)



## 使用“Continue”语句

“continue”语句是控制脚本运行方式的内置命令。除了 bash 脚本之外，它还用于 Python 和 Java 等编程语言。

continue 语句在满足特定条件时停止循环内的当前迭代，然后恢复迭代。

考虑如下所示的 for 循环。

 ```bash
 #!/bin/bash
 for n in {1..10}
 do
         if [[ $n -eq '6' ]]
         then
               echo "Target $n has been reached"
               continue
         fi
         echo $n
 done
 ```

![](https://s2.loli.net/2023/06/15/fVbM81crCgQLZAj.png)



这是代码的作用：

- 第 2 行：标记 for 循环的开始，并将变量 n 从 1 迭代到 10。
- 第 4 行：检查 n 的值，如果变量等于 6，则脚本向标准输出回显一条消息并在第 2 行的下一次迭代中重新启动循环。
- 第 9 行：仅当第 4 行的条件为假时才将值打印到屏幕。

以下是运行脚本后的预期输出。

![](https://s2.loli.net/2023/06/15/sBMG3AuYtrze8IR.png)



## 使用“break”语句

顾名思义，“break”语句会在满足条件时停止或结束迭代。

考虑下面的 For 循环。

```bash
#!/bin/bash
for n in {1..10}
do
        if [[ $n -eq '6' ]]
        then
                echo "Target $n has been reached"
                break
        fi
        echo $n
done
echo "All done"
```

![](https://s2.loli.net/2023/06/15/YD3kxBeOvCyHMrm.png)



这是代码的作用：

- 第 2 行：标记 for 循环的开始，并将变量 n 从 1 迭代到 10。
- 第 4 行：检查 n 的值，如果变量等于 6，则脚本向标准输出回显一条消息并停止迭代。
- 第 9 行：仅当第 4 行的条件为假时才将数字打印到屏幕上。

从输出中可以看出，一旦变量满足循环条件，循环就会停止。



