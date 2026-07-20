Sea $A\in\mathbb{K}^{n \times n}$

~~~python
def calculateLu(A:Matriz) -> Matri:
    """
    Estamos guardando en A tanto L como U
    separados en la diagonal.

    Complejidad: O(n^3)
    """

    n, _ = A.shape
    for k in range(1,n-1):
        pivot = A[k,k]
        if pivot != 0:           
            for i in range(k+1, n):
                m_i = A[i,k]/pivot
                for j in range(k+1, n):
                    A[i,j] = A[i,j] - m_i*A[k,j]


def resolverSistemaTriangularInferior(L, b):
    for i in range(1,n):
        for j in range(1,i):            

            suma = 0
            for k in range(1, i-1):
                suma += L[i,j]*x[j]

            x_i = (b[i] - S) / L[i,i]

~~~