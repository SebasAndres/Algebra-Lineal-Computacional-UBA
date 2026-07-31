# No validado

import numpy as np

def qr_gram_schmidt(A):
    """ Calcular QR con GS """    
    _, N = A.shape
    
    current_basis = [A[:,0]]         
    for j in range(1,N-1):
        col_j = A[:,j]
        proj_j = col_j
        for vec in current_basis:
            proj_j -= np.dot(vec,col_j)/np.dot(vec, vec) @ vec
        current_basis.append(proj_j)

    Q = np.array(current_basis)
    R = Q.T @ A
    
    return Q, R


def qr_householder(A):
    """ Calcular QR con HH """
    _, N = A.shape
    Q = np.eye(N)
    R = A
    for j in range(N):        
        # tienen dim R^{N-j}
        x_j = A[j:,j]
        y_j = np.zeros(N-j)

        v_j = x_j - y_j 
        H_j = np.eye(N-j) - 2/np.dot(v_j, v_j) * np.outer(v_j, v_j)

        if j > 0:
            tmp = np.zeros((N,N))
            tmp[j:,j:] = H_j
            H_j = tmp

        Q = Q @ H_j
        R = Q @ R   
    return Q, R