import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':

    X = np.array([60,85,100,150,250])
    Y = np.array([2.3, 4, 5, 9, 19.5])

    lgX = np.log(X)
    lgY = np.log(Y) 

    coef = np.polyfit(lgX, lgY, deg=1)

    b, a = coef

    model = lambda X: np.exp(a) * np.pow(X, b)
    predY = model(X)

    print("Prediccion para 93kg: ", model(93))

    plt.scatter(X, Y, color='blue')
    plt.scatter(X, predY, color='red')
    plt.show()

