# 如何反转ggplot2中的图例键顺序

在[本教程](https://datavizpyr.com/reverse-legend-key-order-in-ggplot2/ "Source")中，我们将学习如何反转 ggplot2 中图例键的顺序。

在 ggplot2 中，当我们在 aes() 中使用颜色或填充参数为变量着色时，我们会得到一个带有键的图例，显示哪些键匹配哪些颜色。在这里，我们将展示如何使用 guides() 参数为两种类型的图反转图例键的顺序，一种是带有由“颜色”参数制作的图例的散点图，另一种是带有颜色的条形图添加了“填充”参数。

![](https://s2.loli.net/2023/05/09/McESFH6wq1pxCUN.png)



让我们从加载 tidyverse 开始。

```R
library(tidyverse)
theme_set(theme_bw(16))
```

我们将使用 tidyverse 提供的钻石数据。

```R
diamonds %>% head()
```

![](https://s2.loli.net/2023/05/09/XFQCO3PBamjkdWz.png)



## 带彩色点的散点图

让我们在两个变量之间绘制散点图，并使用 aes() 中的颜色参数为第三个（分类）变量着色。

在这里，我们使用从钻石数据中随机抽取的 200 个数据点，使用 slice_sample() 函数制作散点图。

```R
diamonds %>% 
  slice_sample(200) %>%
  ggplot(aes(x=carat, y=price, color=cut))+
  geom_point()
ggsave("how_to_reverse_legend_key_order_legend_with_color.png")

```

这就是使用默认图例键排序的散点图的样子。

![](https://s2.loli.net/2023/05/13/ZT29WGew4Yu1lpj.png)



我们可以使用带有颜色参数的 guides() 函数来反转图例键顺序。我们使用颜色参数来反转，因为我们之前在 aes() 函数中使用颜色参数创建了图例。 reverse = TRUE 的 guide_legend() 函数实际上颠倒了 kegend 键顺序。

```R
diamonds %>% 
  slice_sample(n=200) %>%
  ggplot(aes(x=carat, y=price, color=cut))+
  geom_point()+
  guides(color = guide_legend(reverse = TRUE))
ggsave("reverse_legend_key_order_legend_with_color.png")

```

![](https://s2.loli.net/2023/05/09/sFgpO7wG9QYVndt.png)



## 带填充颜色的条形图

在第二个示例中，让我们制作一个条形图，其中填充了第二个变量指定的颜色。我们在这里使用 aes() 中的 fill 参数来添加颜色，用颜色填充条形图。

```R
diamonds %>% 
  ggplot(aes(cut, fill=clarity))+
  geom_bar()+
  scale_fill_brewer(palette="Dark2")
ggsave("how_to_reverse_legend_key_order_legend_with_fill.png")
```

![image-20230509110948230](https://s2.loli.net/2023/05/09/GrkyjUg1aHWcP4u.png)



我们可以使用 guides() 函数，但这次使用 fill 参数来反转此处的图例键顺序，因为图例是使用 aes() 中的 fill 参数创建的。

```R
diamonds %>% 
  ggplot(aes(cut, fill=clarity))+
  geom_bar()+
  scale_fill_brewer(palette="Dark2")+
  guides(fill = guide_legend(reverse = TRUE))
ggsave("reverse_legend_key_order_for_legend_with_fill.png")
```

![](https://s2.loli.net/2023/05/09/YHMuXZ5neigys7A.png)