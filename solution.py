from fractions import Fraction
from fractions import gcd
from copy import deepcopy

def findterminalidx(m):
        flag = 0
        for i in range(len(m)):
                for j in range(len(m[0])):
                        if m[i][j] != 0:
                                flag = -1
                if flag != -1:
                        flag = i
                        break
                else:
                        flag = 0
        return flag

def computeprobability(m):
        sum = 0
        for i in range(len(m)):
                sum = Fraction(0)
                for j in range(len(m[0])):
                        sum += m[i][j]
                for j in range(len(m[0])):
                                if sum.numerator != 0:
                                        m[i][j] /= sum
        return m

def createidentitymatrix(n):
        m = [[] for i in range(n)]
        for i in range(n):
                m[i] += []
                for j in range(n):
                        m[i] += [Fraction(0)]
                m[i][i] = Fraction(1)
        return m


def IminusQ(I, Q):
        for i in range(len(I)):
                for j in range(len(I[0])):
                        I[i][j] = I[i][j] - Q[i][j]
        return I

def multiplyrow(m, row, k):
    for i in range(len(m[row])):
        m[row][i] *= k
    return m

def addmultipleofrow(m, sourcerow, k, targetrow):
    m1 = multiplyrow(deepcopy(m), sourcerow, k)
    for i in range(len(m[targetrow])):
        if (sourcerow != targetrow):
                m[targetrow][i] += m1[sourcerow][i]
    return m

def invertmatrix(m):
        invertedm =  createidentitymatrix(len(m))
        pass
        for col in range(len(m)):
                diagonalRow = col
                assert(m[diagonalRow][col] != 0)
                k = Fraction(1,m[diagonalRow][col])
                m = multiplyrow(m, diagonalRow, k)
                invertedm = multiplyrow(invertedm, diagonalRow, k)
                sourceRow = diagonalRow
                for targetRow in range(len(m)):
                   if (sourceRow != targetRow):
                        k = -m[targetRow][col]
                        m = addmultipleofrow(m, sourceRow, k, targetRow)
                        invertedm = addmultipleofrow(invertedm, sourceRow,
                                                         k, targetRow)
        return invertedm

def multiplymatrix(m1, m2):
    retVal = []
    rows = len(m1)
    cols = len(m2[0])

    for row in range(rows):
        retVal += [[0]*cols]

    for row in range(rows):
        for col in range(cols):
            num = Fraction(0)
            for i in range(len(m1[0])):
                num += m1[row][i] * m2[i][col]
            retVal[row][col] = num
    return retVal

def lcm(numbers):
        def lcm(a, b):
                return (a * b) / gcd(a, b)
        return reduce(lcm, numbers, 1)

def needswapping(idx1,idx2, m):
        sum1 = sum2 = 0
        for i in range(len(m[0])):
                sum1 += m[idx1][i]
                sum2 += m[idx2][i]
        if sum1 == 0 and sum2 != 0:
                return True
        else:
                return False

def swaprows(idx1, idx2, m):
        for i in range(len(m)):
                temp = m[i][idx2]
                m[i][idx2] = m[i][idx1]
                m[i][idx1] = temp
        temp = m[idx2]
        m[idx2] = m[idx1]
        m[idx1] = temp

def rearrange(m):
        while(True):
                cnt = 0 
                idx1 = 0
                for i in range(len(m)):
                        if needswapping(idx1,i,m):
                                swaprows(idx1,i,m) 
                                cnt += 1
                        idx1 = i
                if cnt == 0:
                        break
        return m

def answer(m):
        m = rearrange(m)
        terminalidx = findterminalidx(m)
        if terminalidx == 0:
                retVal = [0 for x in range(len(m))]
                retVal[0] = 1
                retVal += [1]
                return retVal
        m1 = [[] for x in range(len(m))]
        for i in range(len(m)):
                m1[i] += []
                for j in range(len(m[0])):
                        m1[i] += [Fraction(m[i][j])]
        m1 = computeprobability(m1)
        R = m1[:terminalidx]
        for i in range(len(R)):
                R[i] = R[i][terminalidx:]
        Q = m1[:terminalidx]
        for i in range(len(Q)):
                Q[i] = Q[i][:terminalidx]
        I =  createidentitymatrix(terminalidx)
        I_Q = IminusQ(I, Q)
        I_Q_inverse = invertmatrix(I_Q)
        m = multiplymatrix(I_Q_inverse, R)
        denominators = [x.denominator for x in m[0]]
        numerators = [x.numerator for x in m[0]]
        lcm_denominator =  lcm(denominators)
        retVal = []
        for n, d in zip(numerators, denominators):
                retVal += [n * (lcm_denominator/d)]
        retVal += [lcm_denominator]

        return retVal
        
