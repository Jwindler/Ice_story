# Python日学壹技：性能分析



## 导读

相信日常使用`Python`作为生产力的读者，一定会存在想要分析代码中每一行的运行时间与变量占用内存大小的需求，本文主要分析两个模块，用于分析每行代码的内存使用情况和运行时间情况。



## 内存使用

- [memory-profiler](https://pypi.org/project/memory-profiler/ "memory-profiler")

### 安装

```sh
pip install memory-profiler
```



### 使用方法一

1. 在需要分析的函数上，添加装饰器@profile

```python
@profile
def test1():
	c=0
	for item in xrange(100000):
 		c+=1
 	print (c)
        
```

2. 使用下面的命令运行

```python
python -m memory_profiler memory_profiler_test.py     
```



### 使用方法二

```python
 from memory_profiler import profile
 
 @profile(precision=4,stream=open('memory_profiler.log','w+'))
# @profile
def test1():
     c=0
     for item in xrange(100000):
         c+=1
     print c

# 直接运行即可    
```



### 结果

```txt
Filename: memory_profiler_test.py

Line #    Mem usage    Increment   Line Contents
================================================
     5   21.492 MiB   21.492 MiB   @profile
     6                             def test1():
     7   21.492 MiB    0.000 MiB       c=0
     8   21.492 MiB    0.000 MiB       for item in xrange(100000):
     9   21.492 MiB    0.000 MiB           c+=1
    10   21.492 MiB    0.000 MiB       print c
```

- Mem usage: 内存占用情况
- Increment: 执行该行代码后新增的内存



## 运行时间

- [line-profiler](https://github.com/rkern/line_profiler "line-profiler")

### 安装

```sh
pip install line-profiler
```



### 使用

1. 在需要分析的函数上，添加装饰器@profile

```python
@profile
def slow_function(a, b, c):
    ...
```

2. 运行

```sh
python -m line_profiler script_to_profile.py.lprof
```



### 结果

```
Pystone(1.1) time for 50000 passes = 2.48
This machine benchmarks at 20161.3 pystones/second
Wrote profile results to pystone.py.lprof
Timer unit: 1e-06 s

File: pystone.py
Function: Proc2 at line 149
Total time: 0.606656 s

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
   149                                           @profile
   150                                           def Proc2(IntParIO):
   151     50000        82003      1.6     13.5      IntLoc = IntParIO + 10
   152     50000        63162      1.3     10.4      while 1:
   153     50000        69065      1.4     11.4          if Char1Glob == 'A':
   154     50000        66354      1.3     10.9              IntLoc = IntLoc - 1
   155     50000        67263      1.3     11.1              IntParIO = IntLoc - IntGlob
   156     50000        65494      1.3     10.8              EnumLoc = Ident1
   157     50000        68001      1.4     11.2          if EnumLoc == Ident1:
   158     50000        63739      1.3     10.5              break
   159     50000        61575      1.2     10.1      return IntParIO
```

- 每列含义

> - Line #: The line number in the file.
> - Hits: The number of times that line was executed.
> - Time: The total amount of time spent executing the line in the timer's units. In the header information before the tables, you will see a line "Timer unit:" giving the conversion factor to seconds. It may be different on different systems.
> - Per Hit: The average amount of time spent executing the line once in the timer's units.
> - % Time: The percentage of time spent on that line relative to the total amount of recorded time spent in the function.
> - Line Contents: The actual source code. Note that this is always read from disk when the formatted results are viewed, *not* when the code was executed. If you have edited the file in the meantime, the lines will not match up, and the formatter may not even be able to locate the function for display.
