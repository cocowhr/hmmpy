# coding=utf-8
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
import numpy
import os
import time
ITEMNUM = 30000
ITEMLEN = 5
TV = 0.5  # 短序列匹配阈值

# viterbi算法
# 求以特定状态s开始的长度为T的极大似然状态短序列
def getConnection():
    conn = pymysql.connect(host='localhost', db='supervision', user='root', passwd='root', port=3306,
                           charset='utf8')  # 之后可以放到配置文件中读取
    return conn


def viterbi(A, N, s, pi, len):
    str = [0] * len
    str[0] = s
    delta = [[0 for i in range(N)] for i in range(len)]
    psi = [[0 for i in range(N)] for i in range(len)]
    for i in range(0,N):
        delta[1][i] = pi[i] + A[s][i]
        psi[1][i] = s
    for t in range(2,len):
        for j in range(0,N):
            minval = delta[t - 1][0] + (A[0][j])
            minvalind = 0
            for i in range(1,N):
                val = delta[t - 1][i] + (A[i][j])
                if val > minval:
                    minval = val;
                    minvalind = i;
            delta[t][j] = minval;
            psi[t][j] = minvalind;
    pprob = delta[len - 1][0]
    str[len - 1] = 0
    for i in range(1,N):
        if delta[len - 1][i] > pprob:
            pprob = delta[len - 1][i]
            str[len - 1] = i
    t = len - 2
    while t >= 1:
        str[t] = psi[t + 1][str[t + 1]]
        t -= 1
    return str;


# Testacc 函数
# 测试用
# 假设序列s与短序列长度相同，返回相似程度
def Testacc(s, sequences, A, pi, len):
    for index in range(len(s)):
        i = 0
        while i < len(sequences):
            if sequences[i][0] == s[index]:
                s[index] = i
                break
            i += 1
        if i == len(sequences):
            s[index] = len(sequences) - 1
    v = viterbi(A, len(sequences), s[0], pi, len)
    acc = 0.0
    for i in range(len(s)):
        if s[i] == v[i]:
            acc += 1
    return acc / len(s)


# Test函数
# 以1到T为滑动窗口,
# 依次从待测状态序列中截取从滑动窗口开始到结尾的短序列;
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
    for i in range(len(s) - 1):
        slen = len(s) - i
        v = viterbi(A, len(sequences), s[i], pi, slen)
        acc = 0.0
        for j in range(1, slen):
            if s[i + j] == v[j]:
                acc += 1
        if (acc / (slen - 1)) >= TV:
            sup += 1
    return sup / (len(s) - 1)


def readref(target):
    ref = []
    conn = getConnection()
    cur = conn.cursor()
    sql = "SELECT * FROM "
    sql += target
    cur.execute(sql)
    res = cur.fetchall()
    global ITEMNUM
    ITEMNUM = len(res)
    global ITEMLEN
    ITEMLEN = len(res[0]) - 1
    for row in res:
        temp = [0] * ITEMLEN
        for index in range(len(row) - 1):
            temp[index] = row[index + 1]
        ref += [temp]
    cur.close()
    conn.commit()
    conn.close()
    return ref


# preprocess函数：
# 将系统行为按出现的概率从小到大进行排序
# 计算使∑P(i)≥ 0.9成立的最小正整数N
# 将编号为1～ N -1的行为作为一个状态,将剩余的其他行为(即罕用行为)合并为马氏链的一个状态
# 可以放在数据库中，用View完成该操作
def preprocess(items, sequences, reverse):
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
    while major < len(itemmap):
        sequences += [(itemmap[major][0], itemmap[major][1] / total)]
        totalprecent += itemmap[major][1] / total
        major += 1
    minorpercent = 1.0 / sys.maxint
    sequences += [(-1, minorpercent)]
    ###加入稀有状态###
    # while totalprecent <= 0.99 and major < len(itemmap):
    #     sequences += [(itemmap[major][0], itemmap[major][1] / total)]
    #     totalprecent += itemmap[major][1] / total
    #     major += 1
    # minorpercent = 1.0 / sys.maxint
    # if major < len(itemmap):
    #     minor = len(itemmap) - 1
    #     while minor >= major:
    #         # k = 0
    #         # while k < ITEMNUM:
    #         #     m = 0
    #         #     while m < ITEMLEN:
    #         #         if items[k][m] == itemmap[minor][0]:
    #         #             items[k][m] = -1
    #         #         m += 1
    #         #     k += 1
    #         minor -= 1
    #         minorpercent += itemmap[minor][1] / total
    # sequences += [(-1, minorpercent)]
    ### ###
    for i in range(len(sequences)):
        reverse[sequences[i][0]] = i


# B:辅助序列 为了获得状态转移矩阵A
def getB(sequences, ref):
    B = [0] * len(sequences)
    for k in range(ITEMNUM):
        for m in range(ITEMLEN - 1):
            for i in range(len(sequences) - 1):  # 去掉最后的-1
                if sequences[i][0] == ref[k][m]:
                    B[i] += 1
                    break
                if i == len(sequences) - 2:  # 为稀有情况
                    B[len(sequences) - 1] += 1
    return B


def getA(sequences, ref, target):
    if os.path.exists(target + 'A.csv'):
        A = numpy.loadtxt(open(target + "A.csv", "rb"), delimiter=",", skiprows=0)
    else:
        A = [[0 for i in range(len(sequences))] for i in range(len(sequences))]
        B = getB(sequences, ref)
        for k in range(ITEMNUM):
            for m in range(ITEMLEN - 1):
                for i in range(len(sequences)):
                    for j in range(len(sequences) - 1):  # 去掉最后的-1
                        if i == len(sequences) - 1:  # 为稀有情况
                            if sequences[j][0] == ref[k][m + 1]:
                                A[len(sequences) - 1][j] += 1
                                break
                            elif j == len(sequences) - 2:  # 为稀有情况
                                A[len(sequences) - 1][len(sequences)-1] += 1
                                break
                        elif sequences[i][0] == ref[k][m]:
                            if sequences[j][0] == ref[k][m + 1]:
                                A[i][j] += 1
                                break
                            elif j == len(sequences) - 2:  # 为稀有情况
                                A[i][len(sequences) - 1] += 1
                                break


        # i = 0
        # while i < len(sequences):
        #     j = 0
        #     while j < len(sequences):
        #         k = 0
        #         while k < ITEMNUM:
        #             m = 0
        #             while m < ITEMLEN - 1:
        #                 if sequences[i][0] == ref[k][m]:
        #                     if sequences[j][0] == ref[k][m + 1]:
        #                         A[i][j] += 1
        #                 m += 1
        #             k += 1
        #         j += 1
        #     i += 1
        for i in range(len(sequences)):
            for j in range(len(sequences)):
                if A[i][j] != 0 and B[i] != 0:
                    k = (float)(A[i][j]) / B[i]
                    A[i][j] = math.log(k)
                else:
                    A[i][j] = -sys.maxint
    return A


def saveA(A, target):
    numpy.savetxt(target + 'A.csv', A, delimiter=',')


# getpi:
# 获得初始状态pi
def getpi(sequences):
    pi = []
    for seq in sequences:
        pi += [math.log(seq[1])]
    return pi


def getviterbibyname(target, s):
    conn = getConnection()
    cur = conn.cursor()
    ref = readref(target)
    sequences = []
    reverse = {}
    preprocess(ref, sequences, reverse)
    ###处理s###
    sql = "SELECT `id` FROM "
    sql += target + "_codes"
    sql += " where `name`= '%s'" % (s)
    cur.execute(sql)
    res = cur.fetchall()
    s = reverse[res[0][0]]
    ######
    # print sequences
    A = getA(sequences, ref, target)
    saveA(A, target)
    # for a in A:
    #     print a
    pi = getpi(sequences)
    v = viterbi(A, len(sequences), s, pi, 5)
    for i in range(len(v)):
        if sequences[v[i]][0] != -1:
            sql = "SELECT `name` FROM "
            sql += target + "_codes"
            sql += " where `id`= '%s'" % (sequences[v[i]][0])
            cur.execute(sql)
            res = cur.fetchall()
            v[i] = res[0][0]
        else:
            v[i] = unicode('稀有状态', 'utf-8')
    return v


def getviterbi(target, s):
    conn = getConnection()
    cur = conn.cursor()
    ref = readref(target)
    sequences = []
    reverse = {}
    preprocess(ref, sequences, reverse)
    v = []
    A = getA(sequences, ref, target)
    saveA(A, target)
    pi = getpi(sequences)
    for si in s:
        si = reverse[si]
        vi = viterbi(A, len(sequences), si, pi, 5)
        for i in range(len(vi)):
            if sequences[vi[i]][0] != -1:
                sql = "SELECT `name` FROM "
                sql += target + "_codes"
                sql += " where `id`= '%s'" % (sequences[vi[i]][0])
                cur.execute(sql)
                res = cur.fetchall()
                vi[i] = res[0][0]
            else:
                vi[i] = unicode('稀有状态', 'utf-8')
        v += [vi]
    return v


def getsimilaritybynumber(target, s):
    ref = readref(target)
    sequences = []
    reverse = {}
    preprocess(ref, sequences, reverse)
    # print sequences
    A = getA(sequences, ref, target)
    saveA(A, target)
    # for a in A:
    #     print a
    pi = getpi(sequences)
    return Test(s, sequences, A, pi)


def getsimilaritybyname(target, s):
    conn = getConnection()
    cur = conn.cursor()
    for i in range(len(s)):
        sql = "SELECT `id` FROM "
        sql += target + "_codes"
        sql += " where `name`= '%s'" % (s[i])
        cur.execute(sql)
        res = cur.fetchall()
        s[i] = res[0][0]
    cur.close()
    conn.commit()
    conn.close()
    return getsimilaritybynumber(target, s)


if __name__ == '__main__':
    target = "warning_information"
    # v = getviterbi(target, [1, 2, 3, 4])
    # for v1 in v:
    #     for v11 in v1:
    #         print v11.encode('gbk')
    #     print "_______"
    #target="data_process_fileinfo_type"

    start = time.clock()
    ref = readref(target)
    sequences = []
    reverse = {}
    preprocess(ref, sequences, reverse)
    v = []
    A = getA(sequences, ref, target)
    end = time.clock()
    print str(end - start) + "s"
    v = getviterbi(target, [1,2,3,4])
    for v1 in v:
        for v11 in v1:
            print v11.encode('gbk')
        print "_______"