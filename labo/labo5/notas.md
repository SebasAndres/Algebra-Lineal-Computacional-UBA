# Labo 5

Calcular factoriazacion $LU$ es $O(n³)$

Resolver simplemente $Ax=b$ con eliminacion gaussiana es $O(n²)$

Pero conviene LU porque evita invertir una matriz.

Resolver $LUx = b$ es numericamente mas estable que $x = A^{⁻1}b$

Calcular $LU$ requiere menos flops (operaciones de punto flotante) que invertir

Si $A$ es epsarsa, $L$, $U$ tambien lo son pero $A^{-1}$ no necesariamente

Nunca invertir una matrizz!!!!
> Argumentos en : [gregoryundersen.com/blog/2020/12/09/matrix-inversion](https://gregorygundersen.com/blog/2020/12/09/matrix-inversion/)

La formula tradicional del determinante es O(n!)

Pero si la matriz es LU es O(n):
$det(A) = det(L)*det(U) = 1 * det(U) = 1 * \prod_k U_{kk} = \prod_k U_{kk}$

## Tips haciendo la factorizacion LU:

> Cuando el pivote está en la fila $i$, hacer en esa iteracion:
---> $f_j \leftarrow f_j + \sigma f_i, \text{ con } \sigma \in K \text{t.q. } f_j + \sigma f_i = 0 \land j>i$
---> $f_i \leftarrow f_i / [f_i]$


Donde cada transformacion / iteracion la podemos notar con una matriz que se escribe como:
$M_j = I + m * e_1^t$

