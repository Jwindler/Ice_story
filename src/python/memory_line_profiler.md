# performance analysis



```python
# 时间分析（需要加装饰器，不然只有总的时间）
kernprof -l example.py

python -m line_profiler example.py.lprof 


# 内存分析
from memory_profiler import profile
# @profile 装饰函数后，直接运行即可
```

