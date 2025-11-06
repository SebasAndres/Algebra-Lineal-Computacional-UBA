#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Eliminacion Gausianna
"""
import numpy as np

def elim_gaussiana(A):
    cant_op = 0
    m=A.shape[0]
    n=A.shape[1]
    Ac = A.copy() # donde guardo los A^{(k)}
    
    if m!=n:
        print('Matriz no cuadrada')
        return
    
    ## desde aqui -- CODIGO A COMPLETAR
    for j in range(m-1): # j representa la "k-esima iteracion"
        pivot = Ac[j,j]
        if pivot != 0:
            # Para cada valor debajo del pivot
            # en esa columna
            for i in range(j+1, m):
                # Calculo el coeficiente tal que 
                # A[i,j] - m_i * A[j,j] = 0
                m_i = Ac[i,j]/pivot

                # Resto para cada posicion en esa fila el 
                # valor m_i
                for k in range(j+1,m):
                    Ac[i,k] = Ac[i,k] - m_i*Ac[j,k]

                # [Para optimizar espacio]:
                # Guardo ese coeficiente (invertido)
                # en la posicion A[i,j] que quedaría en 0
                Ac[i,j] = m_i

                # Aumento la cantidad de operaciones
                cant_op += 2 * (m - j - 1) + 1

    # Obtenemos la matriz L como la matriz triangular 
    # inferior con 1s en la diagonal
    L = np.eye(m)
    for i in range(1, m):
        for j in range(i):
            L[i, j] = Ac[i, j]

    # Obtenemos la matriz U como la matriz triangular superior
    U = np.zeros((m, n))
    for i in range(m):
        for j in range(i, n):
            U[i, j] = Ac[i, j]

    ## hasta aqui, calculando L, U y la cantidad de operaciones sobre 
    ## la matriz Ac            
    
    return L, U, cant_op


def main():
    n = 7
    B = np.eye(n) - np.tril(np.ones((n,n)),-1) 
    B[:n,n-1] = 1
    print('Matriz B \n', B)

    L,U,cant_oper = elim_gaussiana(B)
    
    print('Matriz L \n', L)
    print('Matriz U \n', U)
    print('Cantidad de operaciones: ', cant_oper)
    print('B=LU? ' , 'Si!' if np.allclose(np.linalg.norm(B - L@U, 1), 0) else 'No!')
    print('Norma infinito de U: ', np.max(np.sum(np.abs(U), axis=1)) )

if __name__ == "__main__":
    main()
    
    
