# coding=gbk
"""
hmm.py
~~~~~~~~~~~~~~~~
使用马尔科夫链挖掘数据
输入：ref 行为模式库
      s   待测序列
输出：
待测序列是否为正常用户行为序列
可能的提示错误：行为模式库数量不足时可能会有异常
"""
import pymysql
import math
import sys

ITEMNUM = 1285
ITEMLEN = 4
SEQUENCELEN = 4  # 短序列长度
TV = 0.5  # 短序列匹配阈值

# viterbi算法
# 求以特定状态s开始的长度为T的极大似然状态短序列

def viterbi(A, N, s, pi):
    str = [0] * SEQUENCELEN
    str[0] = s
    delta = [[0 for i in range(N)] for i in range(SEQUENCELEN)]
    psi = [[0 for i in range(N)] for i in range(SEQUENCELEN)]
    i = 0
    while i < N:
        delta[1][i] = pi[i] + A[s][i]
        psi[1][i] = s
        i += 1

        t = 2
    while t < SEQUENCELEN:
        j = 0
        while j < N:
            minval = delta[t - 1][0] + (A[0][j])
            minvalind = 0
            i = 1
            while i < N:
                val = delta[t - 1][i] + (A[i][j])
                if val > minval:
                    minval = val;
                    minvalind = i;
                i += 1
            delta[t][j] = minval;
            psi[t][j] = minvalind;
            j += 1
        t += 1
    pprob = delta[SEQUENCELEN - 1][0]
    str[SEQUENCELEN - 1] = 0
    i = 1
    while i < N:
        if delta[SEQUENCELEN - 1][i] > pprob:
            pprob = delta[SEQUENCELEN - 1][i]
            str[SEQUENCELEN - 1] = i
        i += 1
    t = SEQUENCELEN - 2
    while t >= 1:
        str[t] = psi[t + 1][str[t + 1]]
        t -= 1
    return str;


# Testacc 函数
# 测试用
# 假设序列s与短序列长度相同，返回相似程度
def Testacc(s, sequences, A, pi):
    for index in range(len(s)):
        i = 0
        while i < len(sequences):
            if sequences[i][0] == s[index]:
                s[index] = i
                break
            i += 1
        if i == len(sequences):
            s[index] = len(sequences) - 1
    v = viterbi(A, len(sequences), s[0], pi)
    acc = 0.0
    for i in range(len(s)):
        if s[i] == v[i]:
            acc += 1
    return acc / len(s)


# Test函数
# #用Markov训练模块得到的正常特征库中短序列的长度T为滑动窗口,
# 依次从待测状态序列中截取长度为T的短序列;
# 然后与正常特征库中具有相同首状态的特征短序列相比较。
# 若两者不相同的状态数大于一定的阈值V(V是介于1～ T的整数,
# 可称之为短序列匹配阈值),则认为此次截取的短序列不正常;
# 若两者不相同的状态数小于等于阈值V,则认为此次截取的短序列正常。
def Test(s, sequences, A, pi):
    for index in range(len(s)):
        i = 0
        while i < len(sequences):
            if sequences[i][0] == s[index]:
                s[index] = i
                break
            i += 1
        if i == len(sequences):
            s[index] = len(sequences) - 1
    sup = 0.0
    for i in range(len(s) - 3):
        v = viterbi(A, len(sequences), s[0], pi)
        acc = 0.0
        for j in range(SEQUENCELEN):
            if s[i + j] == v[j]:
                acc += 1
        if (acc / SEQUENCELEN) > TV:
            sup += 1
    return sup / (len(s) - SEQUENCELEN + 1)


def readref():
    ref = []
    try:
        conn = pymysql.connect(host='localhost', user='root', passwd='root', port=3306, charset='utf8')
        cur = conn.cursor()
        cur.execute("USE genet")
        cur.execute("SELECT * FROM seq2;")
        res = cur.fetchall()
        for row in res:
            temp = [0] * ITEMLEN
            for index in range(len(row) - 1):
                temp[index] = row[index + 1]
            ref += [temp]
        cur.close()
        conn.commit()
        conn.close()
        return ref
    except  Exception:
        print("error")


# preprocess函数：
# 将系统行为按出现的概率从小到大进行排序
# 计算使∑P(i)≥ 0.9成立的最小正整数N
# 将编号为1～ N -1的行为作为一个状态,将剩余的其他行为(即罕用行为)合并为马氏链的一个状态
# 可以放在数据库中，用View完成该操作
def preprocess(items, sequences):
    itemmap = {}
    total = 0
    for item in items:
        for word in item:
            total += 1
            if word not in itemmap:
                itemmap[word] = 1.0
            else:
                itemmap[word] += 1.0
    itemmap = sorted(itemmap.items(), lambda x, y: cmp(x[1], y[1]), reverse=True)
    totalprecent = 0.0
    major = 0
    while totalprecent <= 0.9 and major < len(itemmap):
        sequences += [(itemmap[major][0], itemmap[major][1] / total)]
        totalprecent += itemmap[major][1] / total
        major += 1
    if major < len(itemmap):
        minorpercent = 0
        minor = len(itemmap) - 1
        while minor >= major:
            k = 0
            while k < ITEMNUM:
                m = 0
                while m < ITEMLEN:
                    if items[k][m] == itemmap[minor][0]:
                        items[k][m] = -1
                    m += 1
                k += 1
            minor -= 1
            minorpercent += itemmap[minor][1] / total
        sequences += [(-1, minorpercent)]


# B:辅助序列 为了获得状态转移矩阵A
def getB(sequences):
    B = [0] * len(sequences)
    i = 0
    while i < len(sequences):
        k = 0
        while k < ITEMNUM:
            m = 0
            while m < ITEMLEN - 1:
                if sequences[i][0] == ref[k][m]:
                    B[i] += 1
                m += 1
            k += 1
        i += 1
    return B


def getA(sequences):
    A = [[0 for i in range(len(sequences))] for i in range(len(sequences))]
    B = getB(sequences)
    i = 0
    while i < len(sequences):
        j = 0
        while j < len(sequences):
            k = 0
            while k < ITEMNUM:
                m = 0
                while m < ITEMLEN - 1:
                    if sequences[i][0] == ref[k][m]:
                        if sequences[j][0] == ref[k][m + 1]:
                            A[i][j] += 1
                    m += 1
                k += 1
            j += 1
        i += 1
    i = 0
    while i < len(sequences):
        j = 0
        while j < len(sequences):
            if A[i][j] != 0 and B[i] != 0:
                k = (float)(A[i][j]) / B[i]
                A[i][j] = math.log(k)
            else:
                A[i][j] = -sys.maxint
            j += 1
        i += 1
    return A


# getpi:
# 获得初始状态pi
def getpi(sequences):
    pi = []
    for seq in sequences:
        pi += [math.log(seq[1])]
    return pi


if __name__ == '__main__':
    ref = readref()
    sequences = []
    preprocess(ref, sequences)
    A = getA(sequences)
    pi=getpi(sequences)
    v = viterbi(A, len(sequences), 0, pi)
    for v1 in v:
        print sequences[v1][0]
    s = [4, 5, 26, 1]
    print Test(s, sequences, A, pi)
