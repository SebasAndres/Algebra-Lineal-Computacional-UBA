import numpy as np

def solveTriangularInferior(L: np.ndarray, b: np.ndarray) -> np.ndarray:
    """ Encontrar x: L @ x = b """
    N, _ = L.shape
    x = np.zeros(N)
    for i in range(N):
        factor = sum([L[i,j]*x[j] for j in range(0,i)])
        x[i] = (b[i] - factor) / L[i,i]
    return x


def solveTriangularSuperior(U: np.ndarray, b: np.ndarray) -> np.ndarray:
    """ Encontrar x: U @ x = b """
    N, _ = U.shape
    x = np.zeros(N)
    for i in range(N-1,-1,-1):
        factor = sum([U[i,j]*x[j] for j in range(N-1,i,-1)])
        x[i] = (b[i] - factor) / U[i,i]
    return x

def LU(A):
    ...

def solveSystemWithLU(A,b):
    L, U = LU(A)
    y = solveTriangularInferior(L,b)
    x = solveTriangularSuperior(U,y)
    return x