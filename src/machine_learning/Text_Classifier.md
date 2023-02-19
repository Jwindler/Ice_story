# 循序渐进的基础知识：文本分类器

在 Python 中构建监督机器学习文本分类器的指导指南和流程图



## 引言

构建文本分类器和理解自然语言处理 (NLP) 的世界涉及很多步骤。这些步骤必须按特定顺序执行。如果数据中的目标类别不平衡，则需要更多步骤。从头开始学习这一切可能有点雷区。网上有很多学习资源，但事实证明，要找到涵盖高层次所有内容的整体指南非常棘手。因此，我写这篇[文章](https://towardsdatascience.com/step-by-step-basics-text-classifier-e666c6bac52b "Source")的目的是希望通过 10 个简单的步骤指南为这个过程提供一些透明度。

我将首先提供一个流程图，该流程图包含所有必要的步骤和要理解的关键点，从阐明任务到部署训练有素的文本分类器。

- 首先，什么是文本分类器？

> 文本分类器是一种算法，它学习单词的存在或模式以预测某种目标或结果，通常是一个类别，例如电子邮件是否是垃圾邮件。

在这里值得一提的是，我将专注于使用监督机器学习方法构建文本分类器。另一种方法是使用深度学习方法，例如神经网络。

让我们看一下该流程图。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230219142649395.png)



## 1. 明确任务

这是任何数据科学项目中最重要的步骤之一。确保您已完全理解所问的问题。您是否有可用的相关数据来回答问题？您的方法是否符合利益相关者的期望？如果您需要利益相关者的支持，请不要构建一些难以解释的超级复杂模型。从简单开始，让每个人都和你一起踏上这段旅程。



## 2. 数据质量检查

任何项目的另一个重要步骤。您的模型只会和输入的数据一样好，因此请确保删除重复项并相应地处理缺失值。



## 3. 探索性数据分析 (EDA)

现在我们可以进行一些特定于文本数据的分析。 EDA 就是要了解数据并了解您可以从中得到什么。此步骤的关键点之一是了解目标类分布。您可以使用 pandas .value_counts() 方法或绘制条形图来可视化数据集中每个类的分布。您将能够看到哪些是多数类和少数类。

![](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20230219143028344.png)



模型在处理不平衡数据时表现不佳。该模型通常会忽略少数类，因为根本没有足够的数据来训练模型来检测它们。 las，如果您发现自己的数据集不平衡且严重偏向目标类别之一，那还不是世界末日。这其实很正常。在您的模型构建过程之前了解这一点很重要，这样您就可以在以后进行调整。

不平衡数据集的存在还应该让您考虑应该使用哪些指标来评估模型性能。在这种情况下，“准确性”（正确预测的比例）真的不是你的朋友。假设您有一个包含二元目标类的数据集，其中 80% 的数据标记为“红色”，20% 的数据标记为“蓝色”。您的模型可以简单地预测整个测试集的“红色”，并且仍然有 80% 的准确率。因此，鉴于您的模型可以简单地预测多数类别，模型的准确性可能会产生误导。

可以使用的一些更好的指标是召回率（正确预测的真阳性的比例）、精度（正确预测的阳性预测的比例）或两者的平均值，即 F1 分数。进入模型构建阶段后，请密切注意少数类的这些分数。您将希望提高这些分数。



## 4. 文本预处理

现在开始一些有趣的事情！文本数据可能包含大量对任何机器学习模型都没有用的东西（取决于任务的性质）。这个过程实际上是关于去除数据集中的“噪音”，将单词同质化并将其剥离回裸露的骨骼，以便只保留有用的单词和最终的特征。

通常，您需要删除标点符号、特殊字符、停用词（如“this”、“the”、“and”等词）并将每个词缩减为词条或词干。您可以尝试制作自己的函数，以便在清理数据之前了解数据中的内容。以下面的函数为例：

```python
#  exploring patterns in the text to assess how best to cleanse the data
pat_list = [r'\d', '-', '\+', ':', '!', '\?', '\.', '\\n'] # list of special characters/punctuation to search for in data

def punc_search(df, col, pat):
    """
    function that counts the number of narratives
    that contain a pre-defined list of special
    characters and punctuation
    """
    for p in pat:
        v = df[col].str.contains(p).sum() # total n_rows that contain the pattern
        print(f'{p} special character is present in {v} entries')

punc_search(df, 'text', pat_list)

# the output will look something like this:

"""
\d special character is present in 12846 entries
- special character is present in 3141 entries
\+ special character is present in 71 entries
: special character is present in 1874 entries
! special character is present in 117 entries
\? special character is present in 53 entries
\. special character is present in 16962 entries
\n special character is present in 7567 entries
"""
```

然后，当您对需要从数据中删除的内容有了更好的了解时，可以尝试编写一个函数，一次性为您完成所有工作：

```python
lemmatizer = WordNetLemmatizer()  # initiating lemmatiser object

def text_cleanse(df, col):
    """
    cleanses text by removing special
    characters and lemmatizing each
    word
    """
    df[col] = df[col].str.lower()  # convert text to lowercase
    df[col] = df[col].str.replace(r'-','', regex=True) # replace hyphens with '' to join hyphenated words together
    df[col] = df[col].str.replace(r'\d','', regex=True) # replace numbers with ''
    df[col] = df[col].str.replace(r'\\n','', regex=True) # replace new line symbol with ''
    df[col] = df[col].str.replace(r'\W','', regex=True)  # remove special characters
    df[col] = df[col].str.replace(r'\s+[a-zA-Z]\s+',' ', regex=True) # remove single characters
    df[col] = df.apply(lambda x: nltk.word_tokenize(x[col]), axis=1) # tokenise text ready for lemmatisation
    df[col] = df[col].apply(lambda x:[lemmatizer.lemmatize(word, 'v') for word in x]) # lemmatise words, use 'v' argument to lemmatise versbs (e.g. turns past participle of a verb to present tense)
    df[col] = df[col].apply(lambda x : " ".join(x)) # de-tokenise text ready for vectorisation
```

然后，您可以在清理后的数据上再次运行第一个函数，以检查您想要删除的所有内容是否确实已被删除。

对于那些注意到上述功能的人，不要删除任何停用词，很好地发现。您可以在矢量化过程中通过几步时间删除停用词。



## 5. 训练-测试拆分

这是有自己的子标题的，因为在开始摆弄这些功能之前执行此步骤非常重要。使用 sklearn 的 train_test_split() 函数拆分数据，然后单独保留测试数据，这样就没有数据泄漏的风险。

如果您的数据不平衡，您可以在测试训练拆分中指定一些可选参数（'shuffle' 和 'stratify'），以确保在目标类之间均匀拆分。这可以确保您的少数类不会完全出现在您的训练或测试集中。

```python
# create train and test data split
X_train, X_test, y_train, y_test = train_test_split(df['text'], # features
                                                    df['target'], # target
                                                    test_size=0.3, # 70% train 30% test
                                                    random_state=42, # ensures same split each time to allow repeatability
                                                    shuffle = True, # shuffles data prior to splitting
                                                    stratify = df['target']) # distribution of classes across train and test
```



## 6. 文本矢量化

模型无法解释文字。相反，必须使用称为矢量化的过程将单词转换为数字。矢量化有两种方法；词袋和词嵌入。 Bag of Words 方法寻找文本之间单词的精确匹配，而 Word Embedding 方法考虑单词上下文，因此可以在文本之间寻找相似的单词。

对于 Bag of Words 方法，句子被标记化，然后每个独特的单词成为一个特征。数据集中的每个唯一单词都对应一个特征，其中每个特征都有一个整数，具体取决于该单词在文本中出现的次数（字数统计向量——sklearn 的 CountVectorizer()）或一个表示重要性的加权整数文本中的单词（一个 TF-IDF 向量——sklearn 的 TfidVectorizer()）。

请务必在训练数据上训练向量化器对象，然后使用它来转换测试数据。



## 7. 模型选择

尝试一些分类模型以查看哪种分类模型最适合您的数据是个好主意。然后，您可以使用性能指标来选择最合适的模型进行优化。我通过运行一个 for 循环来做到这一点，该循环使用 cross_validate() 函数迭代每个模型。

```python
#  defining models and associated parameters
models = [RandomForestClassifier(n_estimators = 100, max_depth=5, random_state=42), 
          LinearSVC(random_state=42),
          MultinomialNB(), 
          LogisticRegression(random_state=42)]

kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=1) # With StratifiedKFold, the folds are made by preserving the percentage of samples for each class.

scoring = ['accuracy', 'f1_macro', 'recall_macro', 'precision_macro']

#  iterative loop print metrics from each model
for model in tqdm(models):
    model_name = model.__class__.__name__
    result = cross_validate(model, X_train_vector, y_train, cv=kf, scoring=scoring)
    print("%s: Mean Accuracy = %.2f%%; Mean F1-macro = %.2f%%; Mean recall-macro = %.2f%%; Mean precision-macro = %.2f%%" 
          % (model_name, 
             result['test_accuracy'].mean()*100, 
             result['test_f1_macro'].mean()*100, 
             result['test_recall_macro'].mean()*100, 
             result['test_precision_macro'].mean()*100))
```



## 8. Baseline model

在您为提高这些性能指标而调整所选模型的超参数而得意忘形之前，请停止。在开始优化之前记下模型的性能。您只能通过将模型与基线分数进行比较来知道（并证明）您的模型有所改进。如果您处于被要求介绍您的方法的位置，它可以帮助您获得利益相关者的支持和讲故事。

创建一个空的 DataFrame，然后在每次模型迭代之后，附加您选择的指标以及迭代的次数或名称，这样您就可以清楚地看到您的模型在优化尝试中的进展情况。



## 9. 模型调优——纠正不平衡数据

通常，微调模型可能涉及调整其超参数和特征工程，以提高模型的预测能力。然而，对于本节，我将重点介绍可用于减少类不平衡影响的技术。

除了为少数类收集更多数据外，还有 5 种方法（据我所知）可用于解决类不平衡问题。大多数是特征工程的一种形式，其目的是对少数类进行过采样或对多数类进行欠采样以平衡整体类分布。

让我们快速浏览一下每种方法：

### 9.1. 添加少数类惩罚

分类算法有一个参数，通常称为“class_weight”，您可以在训练模型时指定该参数。这本质上是一个惩罚函数，如果少数类别被错误分类，将给予更高的惩罚以阻止错误分类。您可以选择自动论证，也可以根据类别手动分配罚分。请务必阅读您正在使用的算法的文档。



### 9.2. 过采样少数类

随机过采样涉及从少数类中随机复制示例并将它们添加到训练数据集中以创建均匀的类分布。这种方法可能会导致过度拟合，因为没有生成新的数据点，所以一定要检查这一点。

python 库 imblearn 包含用于过采样和欠采样数据的函数。重要的是要知道任何过采样或欠采样技术仅适用于训练数据。

如果您使用交叉验证方法将数据拟合到模型中，则需要使用管道来确保仅对训练折叠进行过采样。 Pipeline() 函数可以从 imblearn 库中导入。

```python
over_pipe = Pipeline([('RandomOverSample', RandomOverSampler(random_state=42)), 
                      ('LinearSVC', LinearSVC(random_state=42))])

params = {"LinearSVC__C": [0.001, 0.01, 0.1, 1, 10, 100]}

svc_oversample_cv = GridSearchCV(over_pipe, 
                                 param_grid = params, 
                                 cv=kf, 
                                 scoring='f1_macro',
                                 return_train_score=True).fit(X_train_vector, y_train)
svc_oversample_cv.best_score_  # print f1 score
```



### 9.3. 欠采样多数类

上述方法的另一种方法是对多数类进行欠采样，而不是对多数类进行过采样。有些人可能会争辩说，如果你有数据，就不值得删除数据，但这可能是一个值得你自己尝试的选择。同样，imblearn 库具有可供使用的过采样函数。



### 9.4. 合成少数类的新实例

可以使用称为 SMOTE（合成少数过采样技术）的过程生成少数类的新实例，该过程也可以使用 imblearn 库实现。这里有一篇很棒的文章提供了一些实施 SMOTE 的示例。



### 9.5. 文本增强

可以使用现有数据的同义词生成新数据，以增加少数类的数据点数量。方法包括同义词替换和反向翻译（翻译成一种语言并返回原始语言）。

迭代地运行这些平衡处理步骤中的每一个并将分数与您的基线分数进行比较，然后您可以看到哪种方法最适合您的数据。



## 10. 部署经过训练的分类器

现在是时候将经过训练的分类器推入生产环境，并让它在未见过和未标记的数据上发挥其魔力，前提是它已经过测试。



## 总结

使用监督机器学习方法在 Python 中构建文本分类器的 10 个简单步骤。总之，我们了解到：

- 构建文本分类器所需的步骤顺序
- 检查类别分布的重要性以及了解这如何影响模型性能指标
- 文本预处理步骤
- 如何选择合适的模型并记录基线模型性能
- 解决阶级不平衡的方法

