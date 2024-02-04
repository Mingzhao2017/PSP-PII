# -*- coding:utf-8 -*-
# @Time: 2020/12/1  11:40
# @User: Barry W
# @Author: Mingzhao Wang (wangmz2014@163.com)
# @File: physicochemical.py: 
# @Ref:
# The value of physical, chemical and physicochemical properties

import globalConstant
# DNA
# EIIP values of nucleotides
eiipDictDNA = {'A': 0.1260, 'C': 0.1340, 'G': 0.0806, 'T': 0.1335}
eiipDNA = [0.1260, 0.1340, 0.0806, 0.1335]
# EIIP values of dinucleotides
diDNA = globalConstant.di_DNA
eiipDiDNA = []
for diStr in diDNA:
    valueDiStr = 0
    for i in range(2):
        valueNucle = eiipDictDNA.get(diStr[i])
        valueDiStr = valueDiStr + valueNucle
    eiipDiDNA.append(valueDiStr)

# EIIP values of trinucleotides
triDNA = globalConstant.tri_DNA
eiipTriDNA = []
for triStr in triDNA:
    valueDiStr = 0
    for i in range(3):
        valueNucle = eiipDictDNA.get(triStr[i])
        valueDiStr = valueDiStr + valueNucle
    eiipTriDNA.append(valueDiStr)

# 15 DNA dinucleotide physicochemical properties


