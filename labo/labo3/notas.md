# Notas Labo 3 

Para encontrar la matriz asociada a una transformacion lineal, es util ver cómo la TL afecta a los vectores canónicos.

Ej: 
    Si T : R² -> R², ver T(0,1) y T(1,0)
    Sup. T(1,0) = (0.5, 0) y T(0,1) = (0, 0.5)
    --> T = [[0.5, 0], [0, 0.5]]
    --> T^{-1} = [[2, 0], [0, 2]]
    Y despues ver que T.T^{-1} = I


@ = np.dot = multiplicacion matricial

* != @
* : Es el producto "lugar a lugar"

Hay un margen de error en el calculo de operaciones con matrices, entonces cómo valido si A @ B == C?
|-> np.allclose(A @ B, C)

np.eye(n) es la identidad I de dimension n