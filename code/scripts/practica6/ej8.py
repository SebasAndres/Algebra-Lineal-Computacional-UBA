import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser


from mse import resolverCuadradosMinimosEN

if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument(
        "-p",
        default=5
    )
    args = parser.parse_args()

    # Datos del ejercicio 5
    X = np.array(
        [
            [1950],
            [1960],
            [1970],
            [1980],
            [1990],
            [2000]
        ]
    )
    Y_og = np.array([
        17, 
        20.5,
        23.9,
        27.9,
        32.6,
        36.9
    ])

    # Linealizo
    Y = np.log(Y_og)

    p = args.p
    if p != 5:
        print(f"No es interpolación, grado={p}!=5")
    else:
        print(f"Interpolando con grado={p}")

    coef = resolverCuadradosMinimosEN(X,Y,g=5)

    testX = np.linspace(1950, 2010)
    predY = np.polyval(coef, testX)

    plt.scatter(X, Y, color='blue')
    plt.plot(testX, predY, color='blue')
    plt.title("Poblacion en Argentina")
    plt.show()