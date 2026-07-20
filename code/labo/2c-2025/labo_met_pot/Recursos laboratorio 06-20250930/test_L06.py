# Test L06-metpot2k, Aval

import numpy as np

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

    e1v1 = canonico(1, n)-MayAvect
    H1 = np.eye(n) - 2 * np.outer(e1v1, e1v1)/(np.linalg.norm(e1v1, ord=2))

    if n==1: 
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

#### TESTEOS
# Tests metpot2k

S = np.vstack([
    np.array([2,1,0])/np.sqrt(5),
    np.array([-1,2,5])/np.sqrt(30),
    np.array([1,-2,1])/np.sqrt(6)
              ]).T

# Pedimos que pase el 95% de los casos
exitos = 0
for i in range(100):
    D = np.diag(np.random.random(3)+1)*100
    A = S@D@S.T
    v,l,_ = metpot2k(A,1e-15,1e5)
    if np.abs(l - np.max(D))< 1e-8:
        exitos += 1
assert exitos > 95


#Test con HH
exitos = 0
for i in range(100):
    v = np.random.rand(9)
    #v = np.abs(v)
    #v = (-1) * v
    ixv = np.argsort(-np.abs(v))
    D = np.diag(v[ixv])
    I = np.eye(9)
    H = I - 2*np.outer(v.T, v)/(np.linalg.norm(v)**2)   #matriz de HouseHolder

    A = H@D@H.T
    v,l,_ = metpot2k(A, 1e-15, 1e5)
    #max_eigen = abs(D[0][0])
    if abs(l - D[0,0]) < 1e-8:         
        exitos +=1
assert exitos > 95


# Tests diagRH
D = np.diag([1,0.5,0.25])
S = np.vstack([
    np.array([1,-1,1])/np.sqrt(3),
    np.array([1,1,0])/np.sqrt(2),
    np.array([1,-1,-2])/np.sqrt(6)
              ]).T

A = S@D@S.T
SRH,DRH = diagRH(A,tol=1e-15,K=1e5)
assert np.allclose(D,DRH)
assert np.allclose(np.abs(S.T@SRH),np.eye(A.shape[0]),atol=1e-7)



# # Pedimos que pase el 95% de los casos
# exitos = 0
# for i in range(100):
#     A = np.random.random((5,5))
#     A = 0.5*(A+A.T)
#     S,D = diagRH(A,tol=1e-15,K=1e5)
#     ARH = S@D@S.T
#     e = normaExacta(ARH-A,p='inf')
#     if e < 1e-5: 
#         exitos += 1
# assert exitos >= 95



