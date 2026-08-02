import numpy as np

def resolverSistema(V, b):
    U, s, Vt = np.linalg.svd(V, full_matrices=False)
    x = Vt.T @ ((U.T @ b) / s) 
    return x

def solve(X, Y, F):
    N = X.shape[0]
    M = len(F)

    V = np.zeros((N,M))
    for m, f in enumerate(F):
        V[:, m] = f(X)

    coef = resolverSistema(V, Y)
    model = lambda X: sum(
        coef[i] * f(X) for i, f in enumerate(F)
    )

    return model, coef

if __name__ == '__main__':

    X = np.array([
        x for x in range(100)
    ])

    Y = np.array([
        x**2 for x in X
    ])

    F = [
        lambda X: np.pow(X, 2),
        lambda X: X,
        lambda X: X+1
    ]

    model, coefs = solve(X, Y, F)

    pred = model(X)
    error = sum(
        [
            (y-pred[i])**2 for i, y in enumerate(Y) 
        ]
    )

    print("coefs: ", coefs)
    print("MSE: ", error)