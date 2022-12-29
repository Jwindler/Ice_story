# 5个数据相似性指标



## 1. 介绍

相似性度量是许多数据分析和机器学习任务中的重要工具，使我们能够比较和评估不同数据片段之间的相似性。有许多不同的指标可用，每个指标各有利弊，适用于不同的数据类型和任务。

[本文](https://towardsdatascience.com/5-data-similarity-metrics-f358a560855f "Source")将探讨一些最常见的相似性指标并比较它们的优缺点。通过了解这些指标的特点和局限性，我们可以选择最适合我们特定需求的指标，并确保结果的准确性和相关性。



## 2. 指标

### 2.1. 欧几里得距离

该指标计算 n 维空间中两点之间的直线距离。它常用于连续的数值数据，易于理解和实现。但是，它可能对异常值很敏感，并且没有考虑不同特征的相对重要性。

```python
from scipy.spatial import distance

# Calculate Euclidean distance between two points
point1 = [1, 2, 3]
point2 = [4, 5, 6]

# Use the euclidean function from scipy's distance module to calculate the Euclidean distance
euclidean_distance = distance.euclidean(point1, point2)
```



### 2.2. 曼哈顿距离

该指标通过考虑两点坐标在每个维度中的绝对差异并将它们相加来计算两点之间的距离。它对离群点的敏感性不如欧氏距离，但在某些情况下可能无法准确反映点与点之间的实际距离。

```python
from scipy.spatial import distance

# Calculate Manhattan distance between two points
point1 = [1, 2, 3]
point2 = [4, 5, 6]

# Use the cityblock function from scipy's distance module to calculate the Manhattan distance
manhattan_distance = distance.cityblock(point1, point2)

# Print the result
print("Manhattan Distance between the given two points: " + \
      str(manhattan_distance))
```



### 2.3. 余弦相似度

该指标通过考虑角度来计算两个向量之间的相似度。它通常用于文本数据并且可以抵抗向量大小的变化。但是，它没有考虑不同特征的相对重要性。

```python
from sklearn.metrics.pairwise import cosine_similarity

# Calculate cosine similarity between two vectors
vector1 = [1, 2, 3]
vector2 = [4, 5, 6]

# Use the cosine_similarity function from scikit-learn to calculate the similarity
cosine_sim = cosine_similarity([vector1], [vector2])[0][0]

# Print the result
print("Cosine Similarity between the given two vectors: " + \
      str(cosine_sim))Jaccard Similarity
```



### 2.4. Jaccard相似度

该指标通过考虑两个集合的交集和并集的大小来计算两个集合之间的相似性。它通常用于分类数据并且可以抵抗集合大小的变化。但是，它不考虑集合的顺序或元素的频率。

```python
def jaccard_similarity(list1, list2):
    """
    Calculates the Jaccard similarity between two lists.
    
    Parameters:
    list1 (list): The first list to compare.
    list2 (list): The second list to compare.
    
    Returns:
    float: The Jaccard similarity between the two lists.
    """
    # Convert the lists to sets for easier comparison
    s1 = set(list1)
    s2 = set(list2)
    
    # Calculate the Jaccard similarity by taking the length of the intersection of the sets
    # and dividing it by the length of the union of the sets
    return float(len(s1.intersection(s2)) / len(s1.union(s2)))

# Calculate Jaccard similarity between two sets
set1 = [1, 2, 3]
set2 = [2, 3, 4]
jaccard_sim = jaccard_similarity(set1, set2)

# Print the result
print("Jaccard Similarity between the given two sets: " + \
      str(jaccard_sim))
```



### 2.5. 皮尔逊相关系数

该指标计算两个变量之间的线性相关性。它通常用于连续的数值数据，并考虑不同特征的相对重要性。但是，它可能无法准确反映非线性关系。

```python
import numpy as np

# Calculate Pearson correlation coefficient between two variables
x = [1, 2, 3, 4]
y = [2, 3, 4, 5]

# Numpy corrcoef function to calculate the Pearson correlation coefficient and p-value
pearson_corr = np.corrcoef(x, y)[0][1]

# Print the result
print("Pearson Correlation between the given two variables: " + \
      str(pearson_corr))
```

