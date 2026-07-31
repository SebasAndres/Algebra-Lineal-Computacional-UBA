import numpy as np


def metpot_phi(A, epochs=int(1e4)):
    N, _ = A.shape

    c = np.random.random(N)          
    phi = lambda vec: c @ vec

    x_k = np.random.random(N)
    for _ in range(epochs):
        Ax = A @ x_k
        _lambda = phi(Ax) / phi(x_k)
        x_k = Ax / np.linalg.norm(Ax)   

    return _lambda

def main():
    N = 100
    A = np.random.random((N,N))
    l1 = metpot_phi(A)
    eigvals, eigvecs = np.linalg.eig(A)
    print("Lambda 1 estimado: ", l1)
    print("Lambda 1 numpy: ", eigvals[0])

if __name__ == '__main__':
    main()