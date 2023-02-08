#!/usr/bin/env python
# coding: utf-8
# author: Guisen Chen
# datetime: 2022/09/02
# email: thecgs001@foxmail.com

import os
import re
import sys

import pyfastx

usage = f"python {os.path.basename(sys.argv[0])} <file.fq.gz> [optional, defualt: 1000]"

if len(sys.argv) != 2 and len(sys.argv) != 3:
    print(usage)
    sys.exit()

if len(sys.argv) == 2:
    sys.argv.append(1000)

fastq = pyfastx.Fastx(sys.argv[1])

sites = {
    "MseI": "TTATAA",
    "NcoI": "CCATGCATGG",
    "HindIII": "AAGCTAGCTT",
    "MboI/DpnII": "GATCGATC",
    "Arima": "GATCGATC|GA[A,T,G,C]TGATC|GA[A,T,G,C]TA[A,T,G,C]TC|GATCA[A,T,G,C]TC",
}
i = 0
counts = {"NcoI": 0, "MseI": 0, "HindIII": 0, "MboI/DpnII": 0, "Arima": 0}
while i < int(sys.argv[2]):
    i += 1
    record = next(fastq)
    seq = record[1]
    for e, s in sites.items():
        try:
            re.search(s, seq).group()
            counts[e] += 1
        except:
            pass

print("-------------------------------------------------")
candidate = sorted(counts.items(), key=lambda c: c[1], reverse=True)
if candidate[0][0] == "Arima" and candidate[1][0] == "MboI/DpnII":
    if counts["MboI/DpnII"] < counts["Arima"] - counts["MboI/DpnII"]:
        print(f"The predicted restriction enzyme is \033[35m{candidate[0][0]}\033[0m")
    else:
        print(f"The predicted restriction enzyme is \033[35m{candidate[1][0]}\033[0m")
else:
    print(f"The predicted restriction enzyme is \033[35m{candidate[0][0]}\033[0m")
print("-------------------------------------------------")
for e in counts:
    if e != "Arima":
        print(f"{e}: {counts[e]}")
    else:
        print(
            f'{e}: {counts[e]} (It contains {counts["MboI/DpnII"]} result of MboI/DpnII)'
        )
print("-------------------------------------------------")
