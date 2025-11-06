
# Primer parcial

## Bases, sistemas y proyecciones

#### Hallar un sistema de generadores para un subespacio dado un sistema de ecuaciones:
Unificamos las ecuaciones y reescribimos un X generico de la solucion a partir de generadores.
Ejemplo: 
    Si nos quedan $(x1, x2, x3, x4) = (a, -a+2b, 0, c)$ entonces vale que
    $< (1, -1, 0, 0), (0, 2, 0, 0), (0, 0, 0, 1) >$

#### Probar que existe una única TL dadas $T(v)$
Queremos ver que T queda unívocamente definida para todo el dominio,
Basta con mostrar que las V_i dadas son LI y por lo tanto generan una base del dominio de la transformacion.

Si son Li, QED. 
Encontrar la $V_j$ que es combinacion lineal y ver si cumple que 
$V_j = \sum_{i} \alpha_i V_i$
$T(V_j) = \sum_{i} \alpha_i T(V_i)$

#### Hallar la base de la imagen de una transformación
Tengo que definir la transformación para los generadores del dominio y validar si esos vectores ($W_i = T(V_i)$) son LI entre si.
Si lo son, entonces generan todas las soluciones de la imagen y son generadores de la imagen. $Im(T) = < W_1, ..., W_k >$
Sino, entonces saco el que es LD y vuelvo a validar que el resto lo sea repitiendo el proceso hasta que todos sean LI.

#### Intersección de subespacios
Es juntar los requerimientos que cada subespacio.

#### Subespacio suma
...

#### Probar que una matríz es inversible
$A \text{ es inversible} \iff \exists v \in \R^n, v \neq 0 \text{ tal que } Av = 0$

#### Probar que una matriz es definida positiva
$$A \text{ es definida positiva } \iff A=LU \text{ con } U_{ii} > 0 \text{ para } \forall i \in \{1 .. n\}$$

### Proyecciones
Dado un subespacio $S \subset \mathbb{R}^n$ se define una proyección ortogonal sobre $S$ a 
$$P_S: \mathbb{R}^n \rightarrow \mathbb{R}^n$$

Una matríz de proyección cumple que $P² = P$
$e \in Ker(A^T) \iff e \perp c(A)$

#### Ortogonalización
$v . w = 0 \rightarrow v, w \text{ son ortogonales}$
Decimos que $Q$ es matríz ortogonal cuando $QQ^t = I$ ($Q \in \mathbb{R}^{NxN}$) y $Q^t = Q^-1$.

Cuando nos dan una matríz $A$ de columnas $v_j$, si queremos ortogonalizarla aplicamos (por ejemplo) el algoritmo de Gram-Schmidt.

Algoritmo:
1. Tenemos dos vectores $a$ y $b$.
2. Definimos $A = a$ y $B = b - p_A(B)$ con $p_A(b) = \frac{A^Tb}{A^TA}A$ y nos queda $A \perp B$ (lo vemos porque $A^T B = 0$)
3. Tomamos $A' = \frac{A}{|| A ||}, B' = \frac{B}{||B||}$

## Normas, LU, Cholesky y QR.

#### Propiedades de normas matriciales
$$ || A || = \sup_{v \neq 0} \{ \frac{||Mv||}{||v||} \} < 1$$

$$ cond_*(A) = || A ||_* * || A^{-1} ||_* $$

$$ cond_*(A) \geq \sup_{H \text{ singular}} \{ \frac{|| A || } {|| A - H || }\}$$

$$ cond_*(A) \leq \inf_{H \text{ singular}} \{ \frac{|| A - H|| } {|| A || }\}$$

### LU
Operar / Triangular la matriz usando como pivote la diagonal.

### LU con permutación
Armar una matríz de permutación $P$ de modo que la diagonal no tenga 0s y luego aplicar LU.

### Descomposición Cholesky ($A=\hat{L}\hat{L}^t$)
Pasos:
i. Calculo $A=LU$
ii. Defino D de la siguiente forma:
    $$ D_{ii} = U_{ii} \land D_{ij} = 0 (\forall i \neq j)$$
iii. Defino $D_1$ de tal forma que $D_1 * D_1 = D$, o sea:
    $$ [D_1]_{ii} = \sqrt{D_{ii}} $$ 
iv. Defino $\hat{L} = L*D_1$
v. Luego $\hat{L}^t$ es simplemente tomar la transpuesta de la calculada en el paso anterior y tengo $A=\hat{L}\hat{L}^t$

### Hallar la descomposición QR (con Householder)
Podemos calcularla usando Gram-Schmidt o Householder.

Tenemos $A$ y tomamos $v_i$ como la iésima columna de A `v_i = (A[i:][])`, 
y $w_i = \begin{bmatrix}
        \| v_i \| \\
        0 \\
        \vdots \\
        0
      \end{bmatrix}$

Definimos $u = \frac{v-w}{||v-w||}$
Definimos $H_{u_i} = I - 2 u_i u_i^t$
Vemos $A' = H_{u_i}A = 
    \begin{bmatrix} 
    1 & 0 &  0 \\
    0 & H_{21} & H_{22} \\
    0 & H_{23} & H_{24} 
    \end{bmatrix}$

Y repito el paso hasta tener $A'$ diagonal superior (la llamo $R$),
Defino $Q = H_r * ... * H_2 * H_1$

Luego $A = H_r * ... * H_2 * H_1 * R = QR$


## Autovalores

### Definicion
Son los $v$ tal que $\exists \lambda \in K \text{ tal que } Av = \lambda v$

### Propiedades
> $$ A \text { no es inversible } \iff \lambda=0 \text{ es autovalor}$$
> $$ A\in K^{N\times N} \text { es diagonalizable } \iff \text{ Existen N autovectores LI}$$
> $$ A\in K^{N\times N} \text { es diagonalizable } \iff \text{ Para todo } \lambda_i \text{ autovalor de A, vale que } mg_A(\lambda_i) = ma_A(\lambda_i)$$
> $$A, B \in K^{N\times N} \text{ son semejantes } \iff \exists c \in K^{N \times N} \text{ tal que } A = CBC^{-1}$$ 
> $$ A\in K^{N\times N} \text { es diagonalizable } \iff A = PDP^{-1}$$
> $$ A\in K^{N\times N} \text { es diagonalizable } \iff A^m = PD^mP^{-1}$$
