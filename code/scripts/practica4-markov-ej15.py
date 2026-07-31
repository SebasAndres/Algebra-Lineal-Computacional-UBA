import numpy as np


def crearMatrizMarkovDelCaso(n,p):
    P = np.zeros((n,n)) 
    for k in range(n):
        if k == 0:
            P[:,k] = np.array([1]+[0 for _ in range(n-1)])
        elif k == n-1:
            P[:,k] = np.array([0 for _ in range(n-1)]+[1])
        else:
            P[:, k] = np.array([0 for _ in range(n)])
            P[k-1,k] = 1-p
            P[k+1, k] = p
    return P


def generateRandomV0(n):
    v0 = np.zeros(n)
    v0[np.random.randint(1,n-1)] = 1
    return v0

def main():

    N = 20
    p = 0.5
    NUM_SIMS = 1000

    P = crearMatrizMarkovDelCaso(n=N,p=p)
    
    v0 = generateRandomV0(N)    

    # Simulación
    x_SIM = np.linalg.matrix_power(P, NUM_SIMS) @ v0
    autovalores, autovectores = np.linalg.eig(P)

    breakpoint()

if __name__ == '__main__':
    main()