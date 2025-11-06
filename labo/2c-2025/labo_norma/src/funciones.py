from typing import (
    List, Set, Union, Tuple
)
import numpy as np

Vector = Union[List, Set, np.ndarray]
Matriz = Union[List[Vector], Set[Vector], np.ndarray]

def norma(x: Vector, p: int) -> float:
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


def generarVectorRandom(size:int) -> Vector:
    """
    Genera un vector aleatorio de tamaño size
    """
    return np.random.rand(size)
    

def generarVectorRandomNormalizado(size:int, p:int) -> Vector:
    """
    Genera un vector random normalizado
    """
    x = generarVectorRandom(size=size)
    return x/norma(x, p)


def normaMatMC(A:Matriz, q:int, p:int, Np:int) -> Tuple[float, Vector]:
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

def norma1Exacta(A:Matriz) -> float:
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

def normaInfExacta(A:Matriz) -> float:
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
    
def condMc(A: Matriz, p:int, _ITERS=500000) -> float:
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

def variaPerc(b:Vector, perc:int) -> Vector:
    """ 
    Varia un vector un porcentaje `perc` dado
    """
    numeroAleatorioEntre = lambda desde, hasta: \
        np.random.uniform(desde, hasta) 
    return [
        c * numeroAleatorioEntre(1-perc/100,1+perc/100)
        for c in b
    ]

def generarVectorCentrado(
    x0: Vector
) -> Vector:
    """
    Genera un vector centrado a x0
    """
    return 


def normaMatLJ(
    A:Matriz,
    q:int,
    p:int,
    Np:int,
    maxiter=1000,
    variacion=2,
    rate=0.95
) -> float:
    """
    Calcula la norma matricial usando el método
    Luus-Jaakola
    """

    # Generamos Np vectores unitarios random
    vectoresRandoms = [
        generarVectorRandomNormalizado(size=q, p=p)
        for _ in range(Np)
    ]

    # Definimos como x0 al vector que maximiza ||Ax|| 
    x0, norma0 = max(
        [
            (x, norma(A@x, p)) for x in vectoresRandoms
        ],
        lambda y: y[1]
    )

    for _ in range(maxiter):
        # Generamos Np muestras centradas en x0 con valores 
        # uniformemente distribuidos entre x0 + variacion 
        # y x0-variacion. (np.random.uniform)

        muestrasCentradas = [
            generarVectorCentrado(x0) 
            for _ in range(Np)
        ]        

        muestrasCentradasNormalizadas = [
            normaliza(v,p) 
            for v in muestrasCentradas
        ]

        x1, norma1 = max(
            [
                (x, norma(A@x, p)) for x in muestrasCentradasNormalizadas
            ],
            lambda y: y[1]
        )

        if norma1 > norma0:
            x0 = x1 
            norma0 = norma1

        continue 


    return 0
