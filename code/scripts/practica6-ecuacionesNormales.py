import numpy as np

def resolverSistema(A, b):
    """ Resuelve Ax = b"""
    return np.linalg.solve(A, b)


def resolverCuadradosMinimosEN(X, Y, g):
    """
    Resuelve cuadrados mínimos usando ecuaciones normales

    Devuelve los parámetros a_0, a_1, ..., a_g que minimizan
    ||f(X) - Y||_2  

    Con f definida como:
        
    Si x_i in R 
        --> f(x) = a0 + a1 x + a2 x^2 + ... + ag x^g
    Si x_i in R^f (f > 1) 
        --> f(x) = a0 + a1 x1 + ... + ag xg

    Parámetros:
    - X: np.array
        N Datos de entrada
    - Y: np.array
        N Salidas esperadas
    - g: int    
        Grado del polinomio a ajustar
    """

    N, f = X.shape

    A = np.zeros((N,g+1))
    A[:, g] = np.ones(N)


    if f == 1:
        for i, x in enumerate(X):
            A[i, :] = np.array([x**j for j in range(g,-1,-1)]).T
    else:
        raise NotImplementedError()

    coef = resolverSistema(
        A.T @ A,
        A.T @ Y
    ) 

    return coef


def nahiveErrorCuadratico(testX, testY, foo):
    predY = foo(testX)
    return np.linalg.norm(predY-testY)


if __name__ == '__main__':

    X = np.array([
        [x] for x in range(100)
    ])
    Y = np.array([
        np.polyval([1,0,1], x) for x in range(100)
    ])

    for i in range(5):
        print("Grado #", i)

        coef = resolverCuadradosMinimosEN(X, Y, g=i)
        foo = lambda X: np.polyval(coef, X)
        print("Coef:", coef.round())

        mse = nahiveErrorCuadratico(X, Y, foo)
        print("MSE: ", mse)
        print()