import numpy as np

def esSimetrica(A):
    N, M = A.shape
    if N != M:
        return False
    for i in range(N):
        for j in range(i, M):
            if A[i,j] != A[j,i]:
                return False
    return True

def norma2(A):
    if esSimetrica(A):
        l1, v1 = metPot(A)
    else:
        l1, v1 = metPot(A.T @ A)
        l1 = np.sqrt(l1)
    n2 = np.abs(l1)
    return n2


