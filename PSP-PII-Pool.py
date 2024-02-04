# -*- coding:utf-8 -*-
# @Time: 2021/7/7  11:02
# @User: Barry W
# @Author: Mingzhao Wang (wangmz2014@163.com)
# @File: PSP-PII-Pool.py: 
# @Ref: pool 最新版本，代码有问题。在集群上无问题，待排查！
import multiprocessing
import os
import sys
import time
import math
import globalConstant
import numpy as np

def make_PSPPII_vector(seq, lenSeq, posiData, negaData, alpha, beta):
    '''
    :param seq: The current sequence
    :param lenSeq: the lengh of sequence
    :param posiData: positive data
    :param negaData: negative data
    :param alpha: the distance of between first nucleotide and second nucleotide
    :param beta: the distance of between second nucleotide and third nucleotide
    :return: The corresponding coded feature under parameter pair (alpha, beta)
    '''
    vectorPSPPII = []
    flag = True
    for i in range(alpha + 1):
        for j in range(beta + 1):
            # alpha+beta <= (l-5)/2,否则跳出
            if i + j > math.floor((lenSeq - 5) / 2):
                break
            print('Alpha: ', i)
            print('Beta: ', j)
            seqValueFre = []
            for k in range(i+j+2, lenSeq-i-j-2):
                # 单核苷酸
                singleNucleCurrent = seq[k]
                singleNucleAlphaPrior = seq[k - i - 1]
                singleNucleAlphaBetaPrior = seq[k - i - j - 2]
                singleNucleAlphaNext = seq[k + i + 1]
                singleNucleAlphaBetaNext = seq[k + i + j + 2]
                # 二核苷酸
                diNucleAlphaForward = seq[k] + seq[k + i + 1]
                diNucleAlphaBetaForward = seq[k] + seq[k + i + j + 2]
                diNucleBetaForward = seq[k + i + 1] + seq[k + i + j + 2]
                diNucleAlphaReverse = seq[k] + seq[k - i - 1]
                diNucleAlphaBetaReverse = seq[k] + seq[k - i - j - 2]
                diNucleBetaReverse = seq[k - i - 1] + seq[k - i - j - 2]
                # 三核苷酸
                triNucleForward = seq[k] + seq[k + i + 1] + seq[k + i + j + 2]
                triNucleReverse = seq[k] + seq[k - i - 1] + seq[k - i - j - 2]
                # 计算正类数据中单核苷酸出现频率
                frePosiSingleNucleCurrent = calculate_frequency_single(k, singleNucleCurrent, posiData)
                frePosiSingleNucleAlphaPrior = calculate_frequency_single(k - i - 1, singleNucleAlphaPrior, posiData)
                frePosiSingleNucleAlphaBetaPrior = calculate_frequency_single(k - i - j - 2, singleNucleAlphaBetaPrior, posiData)
                frePosiSingleNucleAlphaNext = calculate_frequency_single(k + i + 1, singleNucleAlphaNext, posiData)
                frePosiSingleNucleAlphaBetaNext = calculate_frequency_single(k + i + j + 2, singleNucleAlphaBetaNext, posiData)
                # 计算负类数据中单核苷酸出现频率
                freNegaSingleNucleCurrent = calculate_frequency_single(k, singleNucleCurrent, negaData)
                freNegaSingleNucleAlphaPrior = calculate_frequency_single(k - i - 1, singleNucleAlphaPrior, negaData)
                freNegaSingleNucleAlphaBetaPrior = calculate_frequency_single(k - i - j - 2, singleNucleAlphaBetaPrior, negaData)
                freNegaSingleNucleAlphaNext = calculate_frequency_single(k + i + 1, singleNucleAlphaNext, negaData)
                freNegaSingleNucleAlphaBetaNext = calculate_frequency_single(k + i + j + 2, singleNucleAlphaBetaNext, negaData)
                # 计算正类数据中二核苷酸出现频率
                frePosiDiNucleAlphaForward = calculate_frequency_di(k, k + i + 1, diNucleAlphaForward, posiData)
                frePosiDiNucleAlphaBetaForward = calculate_frequency_di(k, k + i + j + 2, diNucleAlphaBetaForward, posiData)
                frePosiDiNucleBetaForward = calculate_frequency_di(k + i + 1, k + i + j + 2, diNucleBetaForward, posiData)
                frePosiDiNucleAlphaReverse = calculate_frequency_di(k, k - i - 1, diNucleAlphaReverse, posiData)
                frePosiDiNucleAlphaBetaReverse = calculate_frequency_di(k, k - i - j - 2, diNucleAlphaBetaReverse, posiData)
                frePosiDiNucleBetaReverse = calculate_frequency_di(k - i - 1, k - i - j - 2, diNucleBetaReverse, posiData)
                # 计算负类数据中二核苷酸出现频率
                freNegaDiNucleAlphaForward = calculate_frequency_di(k, k + i + 1, diNucleAlphaForward, negaData)
                freNegaDiNucleAlphaBetaForward = calculate_frequency_di(k, k + i + j + 2, diNucleAlphaBetaForward, negaData)
                freNegaDiNucleBetaForward = calculate_frequency_di(k + i + 1, k + i + j + 2, diNucleBetaForward, negaData)
                freNegaDiNucleAlphaReverse = calculate_frequency_di(k, k - i - 1, diNucleAlphaReverse, negaData)
                freNegaDiNucleAlphaBetaReverse = calculate_frequency_di(k, k - i - j - 2, diNucleAlphaBetaReverse, negaData)
                freNegaDiNucleBetaReverse = calculate_frequency_di(k - i - 1, k - i - j - 2, diNucleBetaReverse, negaData)
                # 计算正类数据中三核苷酸出现频率
                frePosiTriNucleForward = calculate_frequency_tri(k, k + i + 1, k + i + j + 2, triNucleForward, posiData)
                frePosiTriNucleReverse = calculate_frequency_tri(k, k - i - 1, k - i - j - 2, triNucleReverse, posiData)
                # 计算负类数据中三核苷酸出现频率
                freNegaTriNucleForward = calculate_frequency_tri(k, k + i + 1, k + i + j + 2, triNucleForward, negaData)
                freNegaTriNucleReverse = calculate_frequency_tri(k, k - i - 1, k - i - j - 2, triNucleReverse, negaData)
                # 计算正类数据的前向编码值
                if frePosiSingleNucleCurrent == 0 or frePosiSingleNucleAlphaNext == 0 or\
                        frePosiSingleNucleAlphaBetaNext == 0 or frePosiTriNucleForward == 0:
                    posiSubtractorForward = 0
                else:
                    posiSubtractorForward = math.log((frePosiDiNucleAlphaForward * frePosiDiNucleAlphaBetaForward *
                                                      frePosiDiNucleBetaForward) / (frePosiSingleNucleCurrent *
                                                       frePosiSingleNucleAlphaNext * frePosiSingleNucleAlphaBetaNext *
                                                       frePosiTriNucleForward), 2)
                # 计算正类数据的后向编码值
                if frePosiSingleNucleCurrent == 0 or frePosiSingleNucleAlphaPrior == 0 or\
                    frePosiSingleNucleAlphaBetaPrior == 0 or frePosiTriNucleReverse == 0:
                    posiSubtractorReverse = 0
                else:
                    posiSubtractorReverse = math.log((frePosiDiNucleAlphaReverse * frePosiDiNucleAlphaBetaReverse *
                                                      frePosiDiNucleBetaReverse) / (frePosiSingleNucleCurrent *
                                                      frePosiSingleNucleAlphaPrior * frePosiSingleNucleAlphaBetaPrior *
                                                      frePosiTriNucleReverse), 2)
                # 正类数据均值
                meanPosiFre = (posiSubtractorForward + posiSubtractorReverse) / 2
                # 计算负类数据的前向编码值
                if freNegaSingleNucleCurrent == 0 or freNegaSingleNucleAlphaNext == 0 or \
                        freNegaSingleNucleAlphaBetaNext == 0 or freNegaTriNucleForward == 0:
                    negaSubtractorForward = 0
                else:
                    negaSubtractorForward = math.log((freNegaDiNucleAlphaForward * freNegaDiNucleAlphaBetaForward *
                                                      freNegaDiNucleBetaForward) / (freNegaSingleNucleCurrent *
                                                      freNegaSingleNucleAlphaNext * freNegaSingleNucleAlphaBetaNext *
                                                      freNegaTriNucleForward), 2)
                # 计算负类数据的后向编码值
                if freNegaSingleNucleCurrent == 0 or freNegaSingleNucleAlphaPrior == 0 or\
                    freNegaSingleNucleAlphaBetaPrior == 0 or freNegaTriNucleReverse == 0:
                    negaSubtractorReverse = 0
                else:
                    negaSubtractorReverse = math.log((freNegaDiNucleAlphaReverse * freNegaDiNucleAlphaBetaReverse *
                                                      freNegaDiNucleBetaReverse) / (freNegaSingleNucleCurrent *
                                                      freNegaSingleNucleAlphaPrior * freNegaSingleNucleAlphaBetaPrior *
                                                      freNegaTriNucleReverse), 2)
                # 负类数据均值
                meanNegaFre = (negaSubtractorForward + negaSubtractorReverse) / 2
                # 做差
                valueFre = meanPosiFre - meanNegaFre
                # 当前alpha和beta值下的编码向量
                seqValueFre.append(valueFre)
            if flag:
                vectorPSPPII = seqValueFre
                flag = False
            else:
                vectorPSPPII = np.hstack((vectorPSPPII, seqValueFre))
    return vectorPSPPII

def calculate_frequency_single(first, nucleStr, data):
    '''
    :param first: the first nucleotide
    :param nucleStr: the string of nucleotide
    :param data: positive or negative data
    :return: The frequency of nucleotides at the current position in positive or negative dataset
    '''
    singleNucle = []
    for line in data:
        singleNucle.append(line[first])
    fre = float(singleNucle.count(nucleStr)) / data.size
    return fre

def calculate_frequency_di(first, second, nucleStr, data):
    '''
    :param first: the first nucleotide
    :param second:  the second nucleotide
    :param nucleStr: the string of dinucleotide
    :param data: positive or negative data
    :return: The frequency of dinucleotides at the current position in positive or negative dataset
    '''
    diNucle = []
    for line in data:
        diNucleStr = line[first] + line[second]
        diNucle.append(diNucleStr)
    fre = float(diNucle.count(nucleStr)) / data.size
    return fre

def calculate_frequency_tri(first, second, third, nucleStr, data):
    '''
    :param first: the first nucleotide
    :param second: the second nucleotide
    :param third: the third nucleotide
    :param nucleStr: the string of trinucleotide
    :param data: positive or negative data
    :return: The frequency of trinucleotides at the current position in positive or negative dataset
    '''
    triNucle = []
    for line in data:
        triNucleStr = line[first] + line[second] + line [third]
        triNucle.append(triNucleStr)
    fre = float(triNucle.count(nucleStr)) / data.size
    return fre

def read_data(filePath, nucleType):
    '''
    :param filePath: the path of data
    :param nucleType: DNA or RNA
    :return: data with list type and label
    '''
    seqList = []
    label = []
    if nucleType == 'RNA':
        startTuple = tuple(globalConstant.single_RNA)
    elif nucleType == 'DNA':
        startTuple = tuple(globalConstant.single_DNA)
    with open(filePath) as files:
        for line in files:
            if line.startswith(startTuple):
                line = line.rstrip('\n').strip()
                lenLine = len(line)
                seqList.append(line)
            else:
                seqName = line
                PosiOrNega = seqName[1]
                if PosiOrNega == 'P' or PosiOrNega == '+':
                    label.append(1)
                elif PosiOrNega == 'N' or PosiOrNega == '-':
                    label.append(2)
    return seqList, label, lenLine

def split_data_posi_nega(data, label):
    '''
    :param data: data with positive and negative
    :param label: positive and negative
    :return: positive and negative data
    '''
    posi_index = [i for i, x in enumerate(label) if x == 1]
    nega_index = [i for i, x in enumerate(label) if x == 2]
    posi_data = data[posi_index]
    nega_data = data[nega_index]
    return posi_data, nega_data

def constructPartitions(seqList, numCPU):
    seqs_num = len(seqList)
    batch_num = int(seqs_num // numCPU)
    batches = []
    for i in range(numCPU-1):
        batch = seqList[i * batch_num:(i + 1) * batch_num]
        batches.append(batch)
    batch = seqList[(numCPU-1)*batch_num:]  # 不能整除的情况
    batches.append(batch)
    return batches

def make_PSPPII_vector_by_pool(batchSeqList, lenSeq, posiData, negaData, alpha, beta):
    vectorLines = []
    for seqLine in batchSeqList:
        vectorLine = make_PSPPII_vector(seqLine, lenSeq, posiData, negaData, alpha, beta)
        vectorLines.append(vectorLine)
    return vectorLines

def save_result_by_pool(filePath, alpha, beta, methyType):
    '''
    :param filePath: the path of data
    :param alpha: the distance of between first nucleotide and second nucleotide
    :param beta: the distance of between second nucleotide and third nucleotide
    :return: None
    '''
    startTime = time.time()
    splitPath = filePath.split('/')
    nucleType = methyType.split('_')[0]
    if splitPath[2] == 'RNA_m6A':
        dataName = splitPath[4].split('.')[0]
    elif splitPath[2] == 'DNA_4mC':
        dataName = splitPath[3].split('.')[1]
    elif splitPath[2] == 'DNA_6mA':
        dataName = splitPath[3].split('.')[0]
    elif splitPath[2] == 'RNA_multiple':
        dataName = splitPath[3].split('.')[0]
    seqList, label, lenSeq = read_data(filePath, nucleType)
    npSeqList = np.array(seqList)
    posiData, negaData = split_data_posi_nega(npSeqList, label)
    if alpha == 0 and beta == 0:
        alpha = math.floor((lenSeq - 5) / 2)
        beta = math.floor((lenSeq - 5) / 2)
    else:
        if alpha + beta > math.floor((lenSeq - 5) / 2):
            error_info = 'ERROR: alpha + beta <= (lenSeq - 5) / 2'
            sys.stderr.write(error_info)
            sys.exit()
    results = []
    numCPU = int(multiprocessing.cpu_count())
    batches = constructPartitions(seqList, numCPU)
    pool = multiprocessing.Pool(numCPU)
    for batch in batches:
        temp = pool.apply_async(make_PSPPII_vector_by_pool, (batch, lenSeq, posiData, negaData, alpha, beta))
        results.append(temp)
    pool.close()
    pool.join()
    flag = True
    for vecLine in results:
        vecLineX = vecLine.get()
        if len(vecLineX) != 0:
            if flag:
                dataVectorPSPPII = vecLineX
                flag = False
            else:
                dataVectorPSPPII = np.vstack((dataVectorPSPPII, vecLineX))
    print('Done. Used time: %.2fs' % (time.time() - startTime))
    np.savetxt("./result/data_" + dataName + "_PSP_PII_alpha_" + str(alpha) + "_beta_" + str(beta) + ".csv",
               dataVectorPSPPII, delimiter=" ")
    arrayLabel = np.array(label)
    np.savetxt("./result/label_" + dataName + "_PSP_PII.csv", arrayLabel, fmt='%d', delimiter=" ")

if __name__ == '__main__':
    # methyType = 'RNA_m6A'
    # methyType = 'DNA_4mC'
    # methyType = 'DNA_6mA'
    methyType = 'RNA_multiple'
    alpha = 0  # alpha and beta sets [0,(lenSeq - 5) / 2] when alpha = 0 or beta =0
    beta = 0
    if methyType == 'RNA_m6A':
        fPath = './data/RNA_m6A/non-single-base'  # 非单分辨率
        # fPath = './data/RNA_m6A/single-base'  # 单分辨率
        for i, j, k in os.walk(fPath):
            for fName in k:
                print(fName)
                file_path = fPath + '/' + fName
                save_result_by_pool(file_path, alpha, beta, methyType)
    elif methyType == 'DNA_4mC':
        fPath = './data/DNA_4mC'
        for i, j, k in os.walk(fPath):
            for fName in k:
                print(fName)
                file_path = fPath + '/' + fName
                save_result_by_pool(file_path, alpha, beta, methyType)
    elif methyType == 'DNA_6mA':
        fPath = './data/DNA_6mA'
        for i, j, k in os.walk(fPath):
            for fName in k:
                print(fName)
                file_path = fPath + '/' + fName
                save_result_by_pool(file_path, alpha, beta, methyType)
    elif methyType == 'RNA_multiple':
        fPath = './data/RNA_multiple'
        for i, j, k in os.walk(fPath):
            for fName in k:
                print(fName)
                file_path = fPath + '/' + fName
                save_result_by_pool(file_path, alpha, beta, methyType)
    # 单个数据
    # file_path = './data/DNA_4mC/A.thaliana.fasta'
    # save_result(file_path, alpha, beta, methyType)
    print('Done all.')