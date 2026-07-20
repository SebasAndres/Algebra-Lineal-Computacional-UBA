import numpy as np
import scipy.linalg as la

M=100
N=100

A = np.random.random((M,N))

# SVD
U_svd, S, Vs = np.linalg.svd(A)

# Cholesky
# L_chol = np.linalg.cholesky(A)

# LU
P, L_lu, U_lu= la.lu(A)

# Schur
T, U_schur = la.schur(A, output='real')

breakpoint()
