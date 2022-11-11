# Seqkit

## 序列操作

```shell
seqkit seq [flags] file
```

>-p 	#取互补序列
>
>--dna2rna 	#DNA to RNA
>
>-l	#序列以小写字母输出
>
>-g 	#移除组装序列中的gap
>
>-r 	#取反向序列
>
>--rna2dna	#RNA to DNA
>
>-u	#序列以大写字母输出
>
>-w	#每行指定长度数据序列（default=60）

```sh
#将序列转换为一行输出
seqkit seq ex.fasta -w 0 > test.fasta

#每行输出指定碱基n
seqkit seq -w n ex.fasta

#DNA序列转换为RNA序列
seqkit seq --dna2rna ex.fasta

#取反向互补，切每行100碱基
seqkit seq -w 100 -p -r ex.fasta > test.fasta
```



## 格式转换

```sh
# fastq 转换为 fasta
seqkit fq2fa ex1.fq -o ex2.fa

#FASTA/FASTQ转换成tab格式。seqkit fx2tab
seqkit fx2tab ex.fa > test.fa.tab.fa
seqkit fx2tab ex.fq > test.fq.tab.fq
```



## 序列信息统计

```sh
#序列碱基含量及序列长度信息统计
seqkit fx2tab [flags]
```

> -B	#输出碱基的含量 Ex: -B AT  -B N
>
> -g 	#print GC content
>
> -l 	#pring sequence length
>
> -n 	# only print names
>
> -i	#print ID instead of full head



```sh
seqkit fx2tab -l -g -n -i -H ex.fasta
```



```sh
#序列长度分布统计
seqkit stat [flags]

seqkit stat ex.fasta
```



## 根据ID提取序列

```sh
seqkit grep

```

>-n 	#匹配整个序列的名字
>
>-s	  #匹配序列
>
>-d 	#pattern/motif 包含简并碱基
>
>-i  	#忽略大小写
>
>-v 	#反向匹配
>
> -p 	#匹配模式，支持连续写多个模式，匹配任一模式即输出
>
>-R 	#匹配位置选择
>
>-r 	#正则使用



```sh
＃选取有起始密码子的序列
seqkit grep -s -r -i -p ^atg ex.fa

＃根据ID提取序列
seqkit grep -f list ex.fa > new.fa

＃简并碱基使用。S 代表C or G.
seqkit grep -s -d -i -p TTSAA

＃匹配限定到某区域
seqkit grep -s -R 1:30 -i -r -p GCTGG＃

```





## motif定位

```sh
seqkit locate [flags]
```

> -d, --degenerate	#pattern/motif contains degenerate base
>
> -i, --ignore-case	#ignore case
>
> -P, --only-positive-strand	#only search at positive strand
>
> -p, --pattern value	#search pattern/motif
>
> -f, --pattern-file string	#pattern/motif file (FASTA format)

```sh
seqkit locate -i -d -p AUGGACUN ex.fa

```



## 多个文件寻找相同的序列

```sh
seqkit common [flags]
```

> -n, --by-name 	#匹配整个序列的名字，包含description部分，而不是序列id
>
> -s, --by-seq 	#match by sequence
>
> -i, --ignore-case 	#ignore case
>
> -m, --md5 	#use MD5 reduce memory usage



```sh
#By ID (default,>后面，空格之前的名字)输出ID名字相同的。
seqkit common test1.fa test2.fa -o common.fasta

#By full name（整个序列的名字，包含description部分）。输出序列名字相同的。
seqkit common test1.fa test2.fa  -n -o common.fasta

#输出要比较的文件中序列相同的序列
seqkit common test1.fa test2.fa -s -i -o common.fasta

#输出要比较的文件中序列相同的序列 (for large sequences)
seqkit common test1.fa test2.fa -s -i -o common.fasta --md5

```



## 文件切割

```sh
seqkit split [flags]
```

> -i, --by-id	#split squences according to sequence ID
>
> -p, --by-part int	#将一个文件分割成N 份
>
> -s, --by-size int	#将一个文件按照N 条序列一个文件进行分割
>
> -O, --out-dir string	#output directory (default value is infile.split)
>
> -2, --two-pass	#two-pass mode to lower memory usage(only FAST)
>
>

```sh
# 将一个文件切割为 4 份
seqkit split ex.fa -p 4
```
