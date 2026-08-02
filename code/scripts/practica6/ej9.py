
import numpy as np


if __name__ == '__main__':

    X = np.array([0,1,2,4,6])
    Y = np.array([200,195,180,120,25])

    coefs = np.polyfit(x=np.exp2(X), y=Y, deg=1)
    g = -2 * coefs[0]

    print("g: ", g)