# ej 7 interpolación

import matplotlib.pyplot as plt
import numpy as np

from mse import (
    resolverCuadradosMinimosEN,
)

if __name__ == '__main__':   

    for i in [3, 5, 10, 15]:
        print("Grado #", i)

        X = np.array([
            [-1 + (2/i) * d] for d in range(i+1)
        ])
        Y = np.array([
            1/(1+25*x[0]**2) for x in X
        ])

        coef = resolverCuadradosMinimosEN(X, Y, g=i)
        print("Coef:", coef.round())

        foo = lambda X: np.polyval(coef, X)

        testX = np.linspace(-1, 1, 1000).reshape(-1, 1)
        testY = np.array([1/(1+25*x[0]**2) for x in testX])
        predY = foo(testX)

        plt.plot(testX, testY, color='blue')
        plt.plot(testX, predY, color='red')

        error_inf = np.max(np.abs(testY - predY.flatten()))
        print("error inf: ", error_inf)
        plt.show()