# 序列操作神器：Seqkit



## 导读

[本文](https://bioinf.shenwei.me/seqkit/ "Ref")将介绍 `SeqKit` ：用于 `FASTA/Q` 文件操作的跨平台和超快工具包，后续提供了一些长用的示例。



## 1. 安装

- `conda` 安装

```sh
conda install -c bioconda seqkit
```

- `Mac` 安装

```sh
brew install seqkit  # 用于苹果电脑
```



## 2. 用法

### 2.1. 序列操作

```shell
seqkit seq [flags] file
```

- 参数

| 参数      | 作用                               |
| --------- | ---------------------------------- |
| -p        | 取互补序列                         |
| --dna2rna | DNA to RNA                         |
| -l        | 序列以小写字母输出                 |
| -g        | 移除组装序列中的gap                |
| -r        | 取反向序列                         |
| --rna2dna | RNA to DNA                         |
| -u        | 序列以大写字母输出                 |
| -w        | 每行指定长度数据序列（default=60） |



```sh
# 将序列转换为一行输出
seqkit seq ex.fasta -w 0 > test.fasta

# 每行输出指定碱基n
seqkit seq -w n ex.fasta

# DNA序列转换为RNA序列
seqkit seq --dna2rna ex.fasta

# 取反向互补，切每行100碱基
seqkit seq -w 100 -p -r ex.fasta > test.fasta
```



### 2.2. 格式转换

- fa2fa

```sh
# fastq 转换为 fasta
seqkit fq2fa ex1.fq -o ex2.fa

# FASTA/FASTQ 转换成 tab 格式
seqkit fx2tab ex.fa > test.fa.tab.fa
seqkit fx2tab ex.fq > test.fq.tab.fq
```



```sh
# 序列碱基含量及序列长度信息统计
seqkit fx2tab [flags]
```

- 参数

| 参数 | 作用                           |
| ---- | ------------------------------ |
| -B   | 输出碱基的含量 Ex: -B AT  -B N |
| -g   | 输出 GC 含量                   |
| -l   | 输出序列长度                   |
| -n   | 仅输出名字                     |
| -i   | 输出ID                         |
| -H   | 输出 header 行                 |



```sh
# 输出序列长度，GC含量，名字，ID
seqkit fx2tab -l -g -n -i -H ex.fasta
```



### 2.3. 序列信息统计

```sh
# 序列长度分布统计
seqkit stat [flags]
```



- 参数

| 参数 | 作用                                                    |
| ---- | ------------------------------------------------------- |
| -a   | 输出所有统计数据，包括 seq 长度的四分位数、sum_gap、N50 |



```sh
# 统计信息
seqkit stats *.f{a,q}.gz

# 结果如下图
```

![示例](https://swindler-typora.oss-cn-chengdu.aliyuncs.com/typora_imgs/image-20221112214017228.png)



### 2.4. 根据ID提取序列

```sh
seqkit grep
```

- 参数

| 参数 | 作用                                             |
| ---- | ------------------------------------------------ |
| -n   | 匹配整个序列的名字                               |
| -s   | 匹配序列                                         |
| -d   | pattern/motif 包含简并碱基                       |
| -i   | 忽略大小写                                       |
| -v   | 反向匹配                                         |
| -p   | 匹配模式，支持连续写多个模式，匹配任一模式即输出 |
| -R   | 匹配位置选择                                     |
| -r   | 使用正则表达式                                   |



```sh
# 选取有起始密码子的序列
seqkit grep -s -r -i -p ^atg ex.fa

# 根据ID提取序列
seqkit grep -f list ex.fa > new.fa

# 简并碱基使用。S 代表C or G.
seqkit grep -s -d -i -p TTSAA

# 匹配限定到某区域
seqkit grep -s -R 1:30 -i -r -p GCTGG＃
```



### 2.5. motif定位

```sh
seqkit locate [flags]
```

- 参数

| 参数 | 作用                                   |
| ---- | -------------------------------------- |
| -d   | pattern/motif contains degenerate base |
| -i   | ignore case                            |
| -P   | only search at positive strand         |
| -p   | search pattern/motif                   |
| -f   | pattern/motif file (FASTA format)      |



```sh
seqkit locate -i -d -p AUGGACUN ex.fa
```



### 2.6. 多个文件寻找相同的序列

```sh
seqkit common [flags]
```

- 参数

| 参数 | 作用                                                  |
| ---- | ----------------------------------------------------- |
| -n   | 匹配整个序列的名字，包含description部分，而不是序列id |
| -s   | match by sequence                                     |
| -i   | 忽略大小写                                            |
| -m   | use MD5 reduce memory usage                           |



```sh
# By ID (default,>后面，空格之前的名字)输出ID名字相同的。
seqkit common test1.fa test2.fa -o common.fasta

# By full name（整个序列的名字，包含description部分）。输出序列名字相同的。
seqkit common test1.fa test2.fa  -n -o common.fasta

# 输出要比较的文件中序列相同的序列
seqkit common test1.fa test2.fa -s -i -o common.fasta

# 输出要比较的文件中序列相同的序列 (for large sequences)
seqkit common test1.fa test2.fa -s -i -o common.fasta --md5
```



### 2.7. 文件切割

```sh
seqkit split [flags]
```

- 参数

| 参数 | 作用                                           |
| ---- | ---------------------------------------------- |
| -i   | split squences according to sequence ID        |
| -p   | 将一个文件分割成N 份                           |
| -s   | 将一个文件按照N 条序列一个文件进行分割         |
| -O   | 输出目录                                       |
| -2   | two-pass mode to lower memory usage(only FAST) |



```sh
# 将一个文件切割为 4 份
seqkit split ex.fa -p 4
```
