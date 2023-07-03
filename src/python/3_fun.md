# 这3个Python 函数你知道吗？



作为21世纪最流行的语言之一，Python当然有很多有趣的功能值得深入探索和研究。今天将介绍其中的三个，每个都从理论上和通过实际示例进行介绍。

我想要介绍这些函数的主要原因是它们可以帮助您避免编写循环。在某些情况下，循环的运行成本可能很高，除此之外，这些函数将有助于提高速度。

以下是[本文](https://towardsdatascience.com/top-3-python-functions-you-dont-know-about-probably-978f4be1e6d "Source")将涵盖的功能：

1. `map()`
2. `filter()`
3. `reduce()`

即使您以前听说过这些功能，通过更多的理论和示例来加强您的知识也没有什么坏处。

因此，事不宜迟，让我们开始吧！



## map

map() 函数接受另一个函数作为参数，以及某种数组。这个想法是将一个函数（作为参数传入的函数）应用于数组中的每个项目。

这派上用场有两个原因：

1. 你不必写一个循环
2. 它比循环更快

让我们看看它的实际效果。我将声明一个名为 num_func() 的函数，它将一个数字作为参数。该数字被平方并除以 2 并原样返回。请注意，操作是任意选择的，您可以在函数内做任何您想做的事情：

![](https://s2.loli.net/2023/06/09/U4HNCxFY9PvG6q2.png)



现在让我们声明一个数字数组，我们要在其上应用 num_func()。请注意 map() 本身将返回一个地图对象，因此您需要将其转换为列表：

![](https://s2.loli.net/2023/06/09/ruMgnEUXAalboN2.png)



似乎该过程已成功完成。这里没有什么开创性的，但尽可能避免循环是一件好事。



## filter

这是另一个可以节省您时间的不错的函数——无论是在编写还是在执行上。顾名思义，这个想法是只将满足特定条件的项目保留在数组中。

就像 map() 一样，我们可以预先声明函数，然后将它与可迭代列表一起传递给 filter()。

让我们看看这个在行动中。我已经声明了一个名为 more_than_15() 的函数，顾名思义，如果作为参数给出的项目大于 15，它将返回 True：

![](https://s2.loli.net/2023/06/09/lgiyYWQwXIu5dac.png)



接下来，我们声明一个数字数组并将它们作为第二个参数传递给 filter() 函数：

![](https://s2.loli.net/2023/06/09/wzqQDdbGvhMXxac.png)



正如预期的那样，只有三个值满足给定条件。再一次，这里没有什么开创性的，但看起来比循环好得多。



## reduce

现在 reduce() 与前两个有点不同。首先，我们必须从 functools 模块中导入它。这背后的主要思想是它将给定的函数应用于项目数组并返回单个值作为结果。

最后一部分很关键——reduce() 不会返回一个项目数组，它总是返回一个值。让我们看一张图来具体说明这个概念。

![](https://s2.loli.net/2023/06/09/PMrj6bqONTQmk2t.png)



这是在案例图不是 100% 清楚的情况下写出的逻辑：

1. 5 加到 10，结果是 15
2. 15 加 12，结果是 27
3. 27 加 18，结果是 45
4. 45 加到 25，结果是 70

70 是返回的值。从代码实现开始，让我们从 functools 模块导入 reduce 函数并声明一个返回两个数字之和的函数：

![](https://s2.loli.net/2023/06/09/OgFU6xjkRn9o1NA.png)



现在我们可以重新访问代码中的图表，并验证一切是否正常工作：

![](https://s2.loli.net/2023/06/09/LgrMDpjPz1FHfwX.png)



暂时不要进入评论部分——我完全知道还有其他方法可以对列表中的项目求和。这只是展示该功能如何工作的最简单示例。