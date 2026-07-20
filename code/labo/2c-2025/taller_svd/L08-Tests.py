import numpy as np
from scipy.linalg import schur

def metpot2k(A, tol, K):
    # Inicializo la cantidad de operaciones en 0
    k = 0 

    # Genero un vector aleatorio de tamaño n
    n, _ = A.shape
    v = np.random.rand(n)

    # Aplico dos veces A 
    w = A@A@v / np.linalg.norm(A@A@v)

    # Actualizo el error inicial
    e = np.dot(w, v)

    # Itero hasta converger o una cantidad máxima de veces
    while (
        np.abs(e-1) > tol and \
        k < K
    ):
        # Actualizo v 
        v = w 

        # Calculo w aplicando dos veces A
        w = A@A@v / np.linalg.norm(A@A@v)

        # Actualizo el error 
        e = np.dot(v,w) 

        # Actualizo el contador de "operaciones"
        k = k + 1

    # Calculo el autovalor correspondiente al autovector "más grande"
    l = w.T @ A @ w
    error = e - 1 # ????

    # Devuelvo autovector, autovalor, cantidad operaciones, errorFinal
    return w, l, (k, error)

def canonico(i, n):
    e = np.zeros(n)
    e[i-1] = 1
    return e

def diagRH(A, tol, K):

    MayAvect, MayAval, _ = metpot2k(A,tol,K)

    n, _ = A.shape

    e1v1 = canonico(1, n) - MayAvect
    H1 = np.eye(n) - 2 * np.outer(e1v1, e1v1) / np.linalg.norm(e1v1)**2

    if n == 2: 
        S = H1
        D = H1 @ A @ H1.T

    else:
        B = H1 @ A @ H1.T
    
        S_ ,D_ = diagRH(B[1:n, 1:n], tol, K)

        D = np.zeros((n,n))
        D[0][0] = MayAval
        D[1:, 1:] = D_

        S1 = np.zeros((n, n))
        S1[0,0] = 1
        S1[1:, 1:] = S_

        S = H1 @ S1

    return S, D


# -------------
# Tests L08
def svd_reducida(A,k="max",tol=1e-15):
    """
    A la matriz de interes (de m x n)
    k el numero de valores singulares (y vectores) a retener.
    tol la tolerancia para considerar un valor singular igual a cero
    Retorna hatU (matriz de m x k), hatSig (vector de k valores singulares) y hatV (matriz de n x k)
    """

    if k=='max':
        m, k = A.shape

    # 1. Calculo los autovalores y autovectores de A*A
    AtA = A.T@A
    S, D = diagRH(AtA, tol, K=1000000)  # A = SDS*
    print("StS: ", S.T @ S)

    D_ = np.zeros((n,n))
    for i in range(n):
        D_[i,i] = np.diag(D)[i]

    D = D_

    # 2. Defino a Sigma en base a D y a V en base a S, quedandonos con las primeras k cols
    hatSig = np.sqrt(D[:k, :k]) 
    hatV = S[:, :k]

    print("D: ", D)
    print("S: ", S)

    # A = U sigma V*
    hatU = np.zeros((m,k))
    for i in range(k):
        hatU[:, i] = AtA @ hatV[:, i] / np.linalg.norm(AtA @ hatV[:, i])


    print("A: ", A)
    print("pred: ", hatU @ hatSig @ hatV.T)

    return hatU, hatSig, hatV

    
# Matrices al azar
def genera_matriz_para_test(m,n=2,tam_nucleo=0):
    if tam_nucleo == 0:
        A = np.random.random((m,n))
    else:
        A = np.random.random((m,tam_nucleo))
        A = np.hstack([A,A])
    return(A)

def test_svd_reducida_mn(A,tol=1e-15):
    m,n = A.shape
    hU,hS,hV = svd_reducida(A,tol=tol)
    nU,nS,nVT = np.linalg.svd(A)
    r = len(hS)+1
    assert np.all(np.abs(np.abs(np.diag(hU.T @ nU))-1)<10**r*tol), 'Revisar calculo de hat U en ' + str((m,n))
    assert np.all(np.abs(np.abs(np.diag(nVT @ hV))-1)<10**r*tol), 'Revisar calculo de hat V en ' + str((m,n))
    assert len(hS) == len(nS[np.abs(nS)>tol]), 'Hay cantidades distintas de valores singulares en ' + str((m,n))
    assert np.all(np.abs(hS-nS[np.abs(nS)>tol])<10**r*tol), 'Hay diferencias en los valores singulares en ' + str((m,n))

for m in [2,5,10,20]:
    for n in [2,5,10,20]:
        for _ in range(10):
            A = genera_matriz_para_test(m,n)
            test_svd_reducida_mn(A)


# Matrices con nucleo

m = 12
for tam_nucleo in [2,4,6]:
    for _ in range(10):
        A = genera_matriz_para_test(m,tam_nucleo=tam_nucleo)
        test_svd_reducida_mn(A)

# Tamaños de las reducidas
A = np.random.random((8,6))
for k in [1,3,5]:
    hU,hS,hV = svd_reducida(A,k=k)
    assert hU.shape[0] == A.shape[0], 'Dimensiones de hU incorrectas (caso a)'
    assert hV.shape[0] == A.shape[1], 'Dimensiones de hV incorrectas(caso a)'
    assert hU.shape[1] == k, 'Dimensiones de hU incorrectas (caso a)'
    assert hV.shape[1] == k, 'Dimensiones de hV incorrectas(caso a)'
    assert len(hS) == k, 'Tamaño de hS incorrecto'
