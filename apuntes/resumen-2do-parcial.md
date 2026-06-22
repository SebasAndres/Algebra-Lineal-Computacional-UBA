# Temas segundo parcial

## Índice

1. [📊 Existencia de Descomposiciones Matriciales](#1)
2. [¿Cómo diagonalizo una matriz?](#2)
   2.1 [Existencia](#2-1)
   2.2 [Procedimiento](#2-2)
   2.3 [Observación](#2-3)
3. [Descomposición Shur](#3)
   3.1 [Procedimiento](#3-1)
4. [Descomposición SVD](#4)
   4.1 [Observaciones](#4-1)
   4.2 [Procedimiento](#4-2)
   4.3 [Propiedades para $A=U\Sigma V^*$](#4-3)
5. [Procesos de Markov](#5)
   5.1 [Proceso de Markov](#5-1)
   5.2 [Matriz de transición](#5-2)
   5.3 [Existencia de estado límite/equilibrio](#5-3)
   5.4 [Convergencia de estado límite/equilibrio](#5-4)
   5.5 [Cadenas reducibles o irreducibles](#5-5)
6. [Cuadrados mínimos](#6)
7. [Métodos Iterativos. Convergencia](#7)
   7.1 [Método Jacobi](#7-1)
   7.2 [Método Gauss-Seidel](#7-2)
   7.3 [Método SOR](#7-3)
   7.4 [Método gradiente](#7-4)
   7.5 [Propiedades](#7-5)
   7.6 [Radio espectral](#7-6)
   7.7 [Más propiedades](#7-7)

---

## 1. 📊 Existencia de Descomposiciones Matriciales

| Descomposición | Forma | Tipo de Matriz Requerida | ¿Existe para **Toda Matriz** $\boldsymbol{A}$? | Notas Clave |
| :--- | :--- | :--- | :--- | :--- |
| **SVD** (Valores Singulares) | $\boldsymbol{A} = \boldsymbol{U} \boldsymbol{\Sigma} \boldsymbol{V}^*$ | Ninguna (Cualquier matriz $m \times n$). | **SÍ** ✅ | Es la factorización más general y universal. |
| **Schur** (Triangularización) | $\boldsymbol{A} = \boldsymbol{U} \boldsymbol{T} \boldsymbol{U}^*$ | Cuadrada ($n \times n$). | **SÍ** (Si es cuadrada) ✅ | $\boldsymbol{T}$ es triangular superior y $\boldsymbol{U}$ es unitaria. |
| **QR** (Ortogonal/Triangular) | $\boldsymbol{A} = \boldsymbol{Q} \boldsymbol{R}$ | Ninguna (Cualquier matriz $m \times n$ con $m \ge n$). | **SÍ** ✅ | $\boldsymbol{Q}$ es ortogonal/unitaria, $\boldsymbol{R}$ es triangular superior. |
| **LU** (Lower/Upper) | $\boldsymbol{A} = \boldsymbol{L} \boldsymbol{U}$ | Cuadrada ($n \times n$). | **NO** ❌ | Solo existe si todos los **menores principales** son distintos de cero. Siempre existe con permutación: $\boldsymbol{P}\boldsymbol{A} = \boldsymbol{L}\boldsymbol{U}$. |
| **Cholesky** | $\boldsymbol{A} = \boldsymbol{L} \boldsymbol{L}^*$ (o $\boldsymbol{L}\boldsymbol{L}^T$) | Cuadrada ($n \times n$), **Hermitiana** (Simétrica), y **Definida Positiva**. | **NO** ❌ | Es la descomposición más restrictiva. Es única si existe. |

## 2. ¿Cómo diagonalizo una matríz?

### 2.1 Existencia
$A \text{ diagonalizable} \\
\iff \text{Los vectores columna de A forman una base} \\
\iff \text{Existen n autovectores LI (vectores columna de A)} \\
\iff \text{Para todo } \lambda_i \text{ autovalor de A, vale que } mg_A(\lambda_i) = ma_A(\lambda_i) \\
\iff \text{Es semejante a una matriz diagonal} \\
\iff A = PDP^{-1} \\
\iff A^m = PD^mP^{-1}$

### 2.2 Procedimiento:
1. Determinar si A es diagonalizable.
2. Hallar autovalores $\lambda_1, ..., \lambda_n$ de A. 
3. Definir $D = \begin{pmatrix} \lambda_1 & 0 & ... & 0 & 0 \\
0 & \lambda_2 & . & 0 & 0 \\
. & . & . & . & . \\
. & . & . & . & \lambda_n \\
\end{pmatrix}$
4. Hallar autovectores $v_1, ..., v_n$ de A. 
5. Definir $Q = (v_1| ... | v_n)$.
6. Escribir $A = QDQ^*$

### 2.3 Observación: 
En toda diagonalizacion de una matriz, la diagonal tiene los autovalores de esa matriz.

## 3. Descomposición Shur

- Toda matriz $A$ es unitariamente semejante a una matríz triangular superior ($\exists U \text{ unitaria y } T \text{ triangular: } A=UTU^*$).

- Cualquier matríz cuadrada $A$ puede ser escrita como $A = QUQ^*$, con $Q$ matríz unitaria ($Q* = Q^{-1}$) y $D$ diagonal.

### 3.1 Procedimiento:
1. Encontrar un autovector
2. Completar una b.o.n para $Q_1$ con el autovector
3. Calcular $Q_1^* A Q_1$
4. Repetir los pasos anteriores para la submatriz del resultado, pisando los valores en las submatrices de $Q$ y $A$.

## 4. Descomposición SVD

La descomposición en valores singulares de una matríz $A \in \mathbb{C} ^{m\times n}$ es un producto de la forma
$$A = U \Sigma V^*$$

Con:
- $U \in \mathbb{C}^{m \times m}$: Las columnas $u_1, ..., u_m$ vienen dadas por la relación $Av_j = \sigma_j u_j$ con $j=1...n$.
- $V \in \mathbb{C}^{n \times n}$: Las columnas son los autovectores (de $A^*A$). 
- $\Sigma \in \mathbb{C}^{m \times n}$: Diagonal real y no negativa.

Pueden pasar dos casos: 
- $m>n$: En tal situacion se completan las columnas de $U$ para tener una b.o.n en $\mathbb{C}^m$.
- $m<n$: Hay varios $v_j$ asociados a un autovalor 0. Si $\sigma_j=0$ para algun $j\leq \min(n,m)$ entonces se puede elegir $u_j$ completando la ortonormalidad de las columnas de $U$.

### 4.1 Observaciones:
- $A^*A=(U\Sigma V^*)^*U\Sigma V^* = V \Sigma^* U^* U \Sigma V^* = V \Sigma^* \Sigma V^*$
- $AA^*=U\Sigma V^*(U\Sigma V^*)^* = U \Sigma V^* V \Sigma^* U^* = U \Sigma \Sigma^* U^*$
- Como $\Sigma$ es diagonal: 
    - **Si es cuadrada** vale que $\Sigma^* \Sigma = \Sigma^2$.
    - $\Sigma^* \Sigma$ y $\Sigma \Sigma^*$ son semejantes a $\Sigma^2$.
- Calcular SVD para $A^*$ sale de $A^* = (U \Sigma V^*)^* = V \Sigma^t U^*$

### 4.2 Procedimiento:
Sea $A\in\mathbb{C}^{m\times n}$:
1. Calcular $A^*A$ (o $AA^*$)
2. Calcular los autovalores y autovectores de $A^*A$
3. Formar las matrices $U, \Sigma, V^T$:
    - $V$ son los autovectores normalizados de $A^*A$.
    - $\Sigma$ es la matriz diagonal $\mathbb{C}^{m\times n}$ con los $\sigma_i = \sqrt{\lambda_i}$ (de $A^*A$ o $AA^*$).
    - $U$ la calculamos con la relación $Av_i = \sigma_i u_i$ para cada $v_i$ $(i = 1...n)$ y luego completando una base ortonormal de para $u_n ... u_m$. 

### 4.3 Propiedades para $A=U\Sigma V^*$
Sea también $r = rk(A) = \# \{ \sigma_i \in \mathbb{K}_{\neq 0}: \text{ valores singulares no nulos de } A\}$

- $||A||_2 = \sigma_1(A)$

- $\text{A es matríz invertible (rango completo y cuadrada)} \rightarrow \text{cond}_2(\boldsymbol{A}) = \frac{\sigma_{\max}}{\sigma_{\min}} = \frac{\sigma_1}{\sigma_n}$

- En términos de SVD ($A = U\Sigma V^T$), una matriz es ortogonal si y solo si todos sus **valores singulares son iguales a 1** ($\Sigma = I$).

- $\max_{\|x\|_2=1} \|Ax\|_2 = \|A\|_2 = \sigma_1$

    $\text{Demo: } \\
    \max_{\|x\|_2=1} \|Ax\|_2 =  \\
    \max_{\|x\|_2=1} \|U\Sigma V^* x\|_2 = \\
    \max_{\|x\|_2=1} \|\Sigma V^* x\|_2 = \\
    \max_{\|x\|_2=1} \|\Sigma y \|_2 =    \\
    \sum_{i=1}^n{\sigma_i^2 y_i^2} \\
    \leq \sigma_1 * \sum_{i=1}^n{y_i^2}
    = \sigma_1 * \|y\|_2 \\
    = \sigma_1 
    $


- $\min_{\|x\|_2=1} \|Ax\|_2 = \sigma_n$

    $\text{Demo: } \\
    \max_{\|x\|_2=1} \|Ax\|_2 =  \\
    \max_{\|x\|_2=1} \|U\Sigma V^* x\|_2 = \\
    \max_{\|x\|_2=1} \|\Sigma V^* x\|_2 = \\
    \max_{\|x\|_2=1} \|\Sigma y \|_2 =    \\
    \sum_{i=1}^n{\sigma_i^2 y_i^2} \\
    \geq \sigma_n * \sum_{i=1}^n{y_i^2}
    = \sigma_n * \|y\|_2 \\
    = \sigma_n 
    $


- $Im(A) = < u_1, ..., u_r >$
- $Im(A^t) = < v_1, ..., v_r >$
- $Nu(A) = < v_{r+1}, ..., v_{n} >$
- $Nu(A^t) = < u_{r+1}, ..., u_{m} > = Nu(A^tA)$

| Tipo de Matriz | Condición de Forma Cuadrática | Condición de Autovalores ($\lambda_i$) | Consecuencia en $\text{det}(A)$ |
| :--- | :--- | :--- | :--- |
| **Definida Positiva (DP)** | $\boldsymbol{x}^* \boldsymbol{A} \boldsymbol{x} > 0$ para todo $\boldsymbol{x} \ne \boldsymbol{0}$ | **Todos $\boldsymbol{\lambda_i > 0}$** | $\text{det}(\boldsymbol{A}) > 0$ |
| **Semidefinida Positiva (SDP)** | $\boldsymbol{x}^* \boldsymbol{A} \boldsymbol{x} \ge 0$ para todo $\boldsymbol{x}$ | **Todos $\boldsymbol{\lambda_i \ge 0}$** | $\text{det}(\boldsymbol{A}) \ge 0$ |


## 5. Procesos de Markov.

#### 5.1 Proceso de Markov
$v^{(k+1)} = Av^{(k)}$

#### 5.2 Matríz de transición
Una matriz de transición $A$ cumple las siguientes propiedades:
- $A_{ij} \geq 0$
- $\sum_i A_ij = 1$, $\forall j$
- 1 es autovalor de $A$.
- Si $\lambda$ es autovalor entonces $|\lambda| < 1$.
- $mg(\lambda) = ma(\lambda)$
- No toda matríz de Markov es diagonalizable, pero si lo es $A = PDP^{-1} \rightarrow v^{(k+1)} = A^k v^{(0)} = PD^{k}P^{-1}v^{(0)}$

#### 5.3 Existencia de estado límite/equilibrio
Un vector $v$ se dice estado de equilibrio si $Av=v$. Es un autovector asociado para $\lambda=1$. Toda matriz de Markov tiene estado de equilibrio: $Av^* = v^*$

#### 5.4 Convergencia de estado límite/equilibrio
El método ($v^{(k+1)}=Av^{(k)}$) converge para cualquier $v^{(0)} \iff \\ \exists ! \lambda_i : \text{ autovalor de A tal que } (\lambda_i = 1 \land (\forall \lambda_j: \text{ autovalor de A: } \lambda_j \neq \lambda_i \rightarrow |\lambda_j|<1)) \\ \iff \exists A^{\infty} \\ \iff \text{La cadena de Markov asociada es irreducible}$ 
 

##### 5.5 Cadenas reducibles o irreducibles
- Se dice que una **cadena de Markov es irreducible** si $dim(Ker(A-I))=dim(E_A(\lambda=1)) = 1$.
- Si la **cadena de Markov es reducible** ($dim(E_A(\lambda=1))$) entonces la convergencia depende de $v^{(0)}$. 

    Si $v^{(0)}$ es ortogonal a los autovectores asociados a los autovalores de módulo 1, entonces $v^{(k)}$ sí puede converger. 
    Si $v^{(0)}$ tiene componente en dirección de esos autovectores, entonces $v^{(k)}$ no converge (por ejemplo, oscila).

## 6. Cuadrados mínimos

Queremos aproximar una solución de $Ax=b$.

$$||Ax - b|| \rightarrow 0$$

Para resolverlo usamos ecuaciones normales
$$A^tAx= A^tb$$
$$x=(A^tA)^{-1}A^tb$$

Props:
- Si existen soluciones para $Ax=b \rightarrow z=A^tb$ es solucion. Si hay infinitas soluciones, $z=A^tb$ es la solucion de norma 2 mínima.
- $Ax=b$ tiene solucion $\iff AA^tb = b \iff b\in Col(A)$
- Las ecuaciones normales no están bien condicionadas por eso las soluciones usan QR o SVD ($\text{cond}(A^tA) = \text{cond}_2(A)²$ y $\text{cond}_2(A) = \frac{\sigma_1}{\sigma_n}$).
- $A \text{ tiene columnas LI} \iff \text{ Cuadrados mínimos tiene única solución}$.

Desarrollo de la fórmula: 
- Aparece minimizando una función de error $L(y, \hat{y}) = (y-\hat{y})^2$ o también $E = \sum_{i=1}^n (y_i-\hat{y_i})^2$.
- Luego, buscamos...
$$\argmin_{a,b}{E(a,b)}=\argmin_{a,b}{\sum_{i=1}^n (y_i-\hat{y_i})^2}=\argmin_{a,b}{\sum_{i=1}^n (y_i-(a*x_i+b))^2}$$
- O también, matricialmente..
$$\argmin_{a,b}{
    \begin{Vmatrix}
    \begin{pmatrix}
    y_0 \\
    ... \\
    y_n
    \end{pmatrix}
    - 
    \begin{pmatrix}
    a*x_0+b \\
    ... \\
    a*x_n+b
    \end{pmatrix}
    \end{Vmatrix}^2_2
} = 
\argmin_{a,b}{
    \begin{Vmatrix}
    \begin{pmatrix}
    y_0 \\
    ... \\
    y_n
    \end{pmatrix}
    - 
    \begin{pmatrix}
    x_0 & 1 \\
    ... \\
    x_n & 1
    \end{pmatrix}

    \begin{pmatrix}
    a \\
    b 
    \end{pmatrix}
    \end{Vmatrix}_2^2
} =
\argmin_{a,b}{
    \begin{Vmatrix}
    Y - A \begin{pmatrix} a \\ b \end{pmatrix} 
    \end{Vmatrix}_2^2
}
$$

Si llamamos $\hat{X}$ a los parámetros $\begin{pmatrix} a \\ b \end{pmatrix}$ entonces la solución apróximada queda..
$$\argmin_{\hat{X}}{
    \begin{Vmatrix}
    A\hat{X}-b
    \end{Vmatrix}_2^2}$$


## 7. Métodos Iterativos. Convergencia.

$$A = D + L + U$$

Definimos un método iterativo de la siguiente forma...
$$X = TX + c$$

Teo: $x^{n} \text{ converge } \iff \rho(T) < 1$

Con: $\rho(T) = \max{ |\lambda| : \lambda \text{ autovalor de T} }$

A veces se escribe $B=T$, como la matríz de iteración.

#### 7.1 Método Jacobi

$$
Ax = (D+L+U)x = Dx+(L+U)x = b\\
Dx = -(L+U)x + b \\
x = -D^{-1}(L+U)x + D^{-1}b \\
x^{(n+1)} = -D^{-1}(L+U)x^{(n)} + D^{-1}b
$$

- Puede paralelizarse ...

~~~py
def metodoJacobi(A: np.ndarray, b: np.ndarray, epochs: int) -> np.ndarray:
	"""
    Estima la solución x usando
    x^{(n+1)} = -D^{-1}(L+U)x^{(n)} + D^{-1}b

    Parámetros
    ----------
    - A: np.ndarray
        Matriz del sistema
    - b: np.ndarray
        Solucion a partir del vector x
    - epochs: int
        Cantidad de iteraciones del método

    Retorna
    -------
    - x: np.ndarray
        Solucion estimada de x
	"""
	m, n = A.shape

	D, L, U = splitMatrix(A)
    D_inv = alc.inv(D)
	Bj = -D_inv@(L+U)
	c = D_inv@b
	
	x = np.random.random((n, 1))
	for i in range(epochs):
		x = Bj @ x + c	
	return x
~~~


#### 7.2 Método Gauss-Seidel
$$
Ax = (D+L+U)x = (D+L)x+Ux = b \\
(D+L)x = - Ux + b \\
x = -(D+L)^{-1}Ux + (D+L)^{-1}b \\
x^{(n+1)} = -(D+L)^{-1}Ux^{(n)} + (D+L)^{-1}b
$$
- No puede implementarse de manera paralelizada.

~~~py
def metodoGaussSeidel(A: np.ndarray, b: np.ndarray, epochs: int) -> np.ndarray:
	"""
    Estima la solución x usando
    x^{(n+1)} = -(D+L)^{-1}Ux^{(n)} + (D+L)^{-1}b

    Parámetros
    ----------
    - A: np.ndarray
        Matriz del sistema
    - b: np.ndarray
        Solucion a partir del vector x
    - epochs: int
        Cantidad de iteraciones del método

    Retorna
    -------
    - x: np.ndarray
        Solucion estimada de x
	"""
	m, n = A.shape

	D, L, U = splitMatrix(A)
    Aux_inv = alc.inv(D+L)
	Bgs = -Aux_inv@U
	c = Aux_inv@b
	
	x = np.random.random((n, 1))
	for i in range(epochs):
		x = Bgs @ x + c	
	return x
~~~

#### 7.3 Método SOR
$s$

#### 7.4 Método gradiente
$s$

#### 7.5 Propiedades
- $T = -M^{⁻1}N \land \lambda \text{ autovalor de T} \iff det(\lambda M + N) = 0$ 
- El método converge en $n$ pasos si $x^{n+1} = x^{n} \iff T^{n+1} = 0$.
- $A$ matriz cuadrada y tridiagonal, entonces $\rho(B_{GS}) = \rho(B_j)^2$
- Para toda norma subordinada $||.||$ vale que: $\rho(B) = lim_{n \rightarrow \inf} || B^n ||^{1/n}$

#### 7.6 Radio espectral
El radio espectral $\rho(T)$ de la matriz de iteración $T$ determina si un método iterativo converge y qué tan rápido lo hace. 
- Si $\rho(B)<1$: el método converge .
- Si $\rho(B)=1$: puede no converger o converger muy lentamente.
- Si $\rho(B)>1$ : el método diverge .
- Cuanto menor sea $\rho(B)$, más rápida será la convergencia .

### 7.7 Más propiedades

- $A$ es hermitiana $\iff$ $A=A^*$. 

- Las matrices hermitianas tienen autovalores reales y son diagonalizables ortogonalmente con $A=UDU^*$

- Las matrices hermitanas son un subconjunto de las normales.

- $A$ hermitiana -> todo sus autovectores son ortogonales entre si.

- $A$ simétrica real $\rightarrow$ $A$ hermitana. 

- $A$ es normal $\iff$ $A^*A = AA^*$

- <u>Teorema Espectral</u>: Si $A$ es normal, entonces se puede diagonalizar en una b.o.n como $A=QDQ^*$. 

- $A^tA$ es simétrica

- $\text{max } ||Ax||_2 = \text{max}_i |\sigma_i|$

- $\text{min } ||Ax||_2 = \text{min}_i |\sigma_i|$

- $||UA|| = ||A||, \text{ U matriz unitaria}$

- Si $M$ es real y simétrica y $v,w$ son autovectores de autovalores diferentes, entonces $v, w$ son ortogonales.

- Pseudoinversa: $A^{+} = (A^tA)^{-1}A^t$


