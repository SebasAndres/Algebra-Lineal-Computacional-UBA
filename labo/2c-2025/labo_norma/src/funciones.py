from typing import (
    List, Set, Union
)
import numpy as np

Vector = Union[List, Set, np.ndarray]
Matriz = Union[List[Vector], Set[Vector], np.ndarray]

def norma(x: Vector, p: int) -> int:
    """
    Calcula la norma p de un vector x dado.
    """
    if isinstance(p, str):
        if p == 'inf':
            return abs(max(x, key=lambda c: abs(c)))
        raise Exception("Norma invalida")

    res = 0
    for cord in x:
        res += abs(cord)**p
    res = res**(1/p)
    return res


def normaliza(x:Vector, p: int) -> Vector:
    """
    Normaliza un vector x dado (norma p)
    """
    _norm:int = norma(x,p)
    y:Vector = []
    for c in x:
        y.append(c/_norm)
    return y    


def generarVectorRandom(size:int):
    """
    Genera un vector aleatorio de tamaño size
    """
    return np.random.rand(size)
    

def generarVectorRandomNormalizado(size:int, p:int):
    """
    Genera un vector random normalizado
    """
    x = generarVectorRandom(size=size)
    return x/norma(x, p)


def normaMatMC(A:Matriz, q:int, p:int, Np:int):
    """
    Método Monte Carlo:
    Estima ||A||q,p usando Np vectores de x de R^n
    generados al azar

    A: matriz original
    q, p: normas 
    Np: cantidad de vectores generados al azar

    Retorna la norma y el vector que la maximiza
    """

    vectoresRandoms = [
        A@generarVectorRandomNormalizado(size=q, p=p)
        for _ in range(Np)
    ]
    return max(
        [
           (r,norma(x=r,p=q))
           for r in vectoresRandoms
        ],
        key=lambda t: t[1]
    )

def norma1Exacta(A:Matriz):
    """
    Calcula la norma 1
    - Es la suma más grande de las columnas en modulo
    """

    if isinstance(A, np.ndarray):
        r,c = A.shape
    else:
        r = len(A)
        c = len(A[0])

    return max([
        sum(
            [abs(A[i][j]) for i in range(r)]
        )
        for j in range(c)
    ])

def normaInfExacta(A:Matriz):
    """
    Calcula la norma inf
    - Es la suma más grande de las filas en modulo
    """
    if isinstance(A, np.ndarray):
        return norma1Exacta(A.T)
    else:
        r = len(A)
        c = len(A[0])
        return max([
            sum(
                [abs(A[i][j]) for j in range(r)]
            )
            for i in range(c)
        ])
    
def condMc(A: Matriz, p:int, _ITERS=500000):
    """
    Estima el nro de condicion (k) usando la norma inducida p
    """
    norm_A = normaMatMC(
        A, 
        p=p,
        q=p, 
        Np=_ITERS
        )[1]
    norm_A_inv = normaMatMC(
        np.linalg.inv(A), 
        p=p, 
        q=p,
        Np=_ITERS
        )[1]
    estimated_k = norm_A * norm_A_inv
    return estimated_k

def variaPerc(b:Vector, perc:int):
    """ 
    Varia un vector un porcentaje dado
    """
    numeroAleatorioEntre = lambda desde, hasta: \
        np.random.uniform(desde, hasta) 
    return [
        c * numeroAleatorioEntre(1-perc/100,1+perc/100)
        for c in b
    ]

