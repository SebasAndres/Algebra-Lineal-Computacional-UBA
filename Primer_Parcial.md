# Primer Parcial

## Propiedades matriciales

- $A \text{ es inversible} \iff [Ax = 0 \iff x=0] \iff det(A) \neq 0$

- $S^⊥ =^{def}$?

- Cualquier matriz $A\in\mathbb{K}^{n\times n}$ puede escribirse como $A=\sum_i \sum_j a_{ij} E_{ij}$, siendo $E_{ij}$ la matríz definida por 
$$E_{ij} = \begin{cases} 
    1 & \text{si i=j} \\
    0 & \text{cc.}
\end{cases}$$


## Proyecciones

Dado un subespacio $S \subset \mathbb{R}^n$ tal que $S = \text{col}(A)$ se define una proyección ortogonal sobre $S$ a $P_S: \mathbb{R}^n \rightarrow \mathbb{R}^n$, tal que:

- $P_S(s) = s \text{ } \forall s\in S $
- $P_S(x) \in S \text{ } \forall x \in S$
- $x−P_S(x) \in S^⊥$

Se define de la siguiente forma:  
$$P_s = A(A^tA)^{-1}A^t$$

Donde a nivel operativo se busca una $\hat A = Q$, formada con los vectores columna de $A$ ortonormalizados para que $\hat A^t \hat A = Q^t Q = I$ entonces:

$$P_S = QQ^t$$

O también, si tengo una b.o.n de $S$ de la forma $\{ v_1, ..., v_r \}$, puedo obtener el proyector ortogonal a S como:

$$P_S(x) = P * x = (\sum_{i=1}^r v_i v_i^t) * x$$ 

### Ortonormalizar una matriz

#### Proceso (Gram-Schmidt)
Para ortonormalizar una matríz se puede usar Gram-Schmidt.

<u> Algoritmo: </u>

Sea $a_1,..., a_n$ las columnas de la matríz origen $A$.

1. Fijamos $v_1 := a_1$.  

3. Para cada $i \in [1..n]$:

- $v_{i+1} = a_{i+1} - \sum_{j=1}^i p_{v_j}(a_{i+1})$

4. Para cada $v_i$ calculado tomamos $v_i' = \frac{v_i}{||v_i||}$

Donde:
-  $p_a(b) = \frac{a^Tb}{a^Ta}a$

#### Propiedades de matrices ortonormales
- $Q^tQ = I$
- $Im(A) = Im(Q) \text{ pero la transformación lineal es distinta, es decir, generalmente } Ax \neq Qx$
- $Q$ es más estable y tiene mejor numero de condicion.

## Normas vectoriales
Una norma de un $\mathbb{K}$-espacio vectorial es una función $||.|| : V \rightarrow \mathbb{R}_{\geq 0}$ que cumple las siguientes propiedades:
1. $||av|| = |a| ||v|| \text{ para } a\in \mathbb{K} \land v \in V$
2. Si $||v|| = 0$, entonces $v=0$.
3. $||u+v|| \leq ||u|| + ||v|| \text{ para todo } u,v \in V$

#### Normas vectoriales (en $\mathbb{K}^n$) comunes:
- Norma-1: $||v|| = |v_1| + ... + |v_n|$
- Norma-infinito: $||v||_{\infty} = \text{máx} \{ |v_1|, ..., |v_n| \}$
- Norma-p: $||v||_p = (|v_1|^p + ... + |v_n|^p)^{1/p}$

##### Equivalencia entre normas:
Sean $||.||$ y $||.||_*$ dos normas en mismo $\mathbb{K}$-espacio vectorial $V$. Son equivalentes si $\exists c,C > 0$ tales que para todo $x\in V$:
$$c||x||_* \leq ||x|| \leq C ||x||_*$$

##### Convergencia de un vector a una norma
Sucede si $||v_n - v|| \rightarrow 0 \text{ cuando } n\rightarrow \infty$.

## Normas matriciales
Dada $A \in \mathbb{K}^{n \times m}$ y un par de normas vectoriales $||.||_n$ y $||.||_m$ en $\mathbb{K} \text{ y } \mathbb{K}^m$. 

$$
\|A\|_{n,m} = \max_{\substack{x \in \mathbb{K}^m \\ x \ne 0}} \frac{\|Ax\|_n}{\|x\|_m}
= \max_{\substack{x \in \mathbb{K}^m, ||x||_m = 1}} ||Ax||_n
$$

#### Normas matriciales comunes
- Norma-infinito: $||A||_{\infty} = \max_{\substack{i \leq i \leq n}} \{ \sum_{j=1}^n |a_{ij}|\}$
- Norma-1: $||A||_1 = \max_{\substack{i \leq j \leq n}} \{ \sum_{j=1}^n |a_{ij}|\}$

#### Propiedades:
- $||Ax|| \leq ||A|| ||x||$
- $||AB|| \leq ||A|| ||B||$
- $||A|| = \sup_{v \neq 0} \{ \frac{||Mv||}{||v||} \} < 1$
- $Q$ unitaria, $||Qx|| = ||x||$
- $cond_*(A) = ||A||_* ||A^{-1}||_* $
- $cond_*(A) \geq \sup_{H \text{ singular}} \{ \frac{||A|| } {||A - H ||}\}$
- $cond_*(A) \leq \inf_{H \text{ singular}} \{ \frac{||A-H|| } {||A|| }\}$

## Descomposición $A=LU$
$$A = LU$$
- $L$: triangular inferior (con la diagonal principal).
- $U$: triangular superior.

Para matrices cuadradas.
Existe sii durante la eliminación gaussiana, todos los pivotes deben ser $\neq 0$ (sin intercambiar filas), esto sucede si $det(A) \neq 0$.

<u>Algoritmo:</u>
Triangular la matríz original con operaciones elementales entre filas hasta que quede triangular superior.
La matriz triangular superior resultante es $U$, a $L$ la construimos agregando en cada columna en el valor $L_{ij}$ el inverso en signo del multiplicador usado para triangular la celda $ij$.

### Descomposición $PA=LU$
Es equivalente al algoritmo de descomposición $LU$, pero se aplica previamente una permutación expresada con el producto con una matríz $P$.
Esta descomposición existe para cualquier matríz cuadrada.

## Descomposición Cholesky $A=\hat{L}\hat{L}^t$

<u>Algoritmo: </u>
1. Calculo $A=LU$
2. Defino D de la siguiente forma:
    $$ D_{ij} = \begin{cases}
    L_{ii} & \text{si i=j} \\
    0 & \text{cc.}
    \end{cases}$$
3. Defino $D_1$ de tal forma que $D_1 * D_1 = D$, o sea:
    $$ (D_1){ii} = \sqrt{D_{ii}} $$ 
4. Defino $\hat{L} = L*D_1$
5. Luego $\hat{L}^t$ es simplemente tomar la transpuesta de la calculada en el paso anterior y tengo $A=\hat{L}\hat{L}^t$

## Descomposición $A=QR$

Podemos calcularla usando Gram-Schmidt o Householder.

<u>Algoritmo:</u>
1. Tenemos $A$ y tomamos $v_i$ como la iésima columna de A `v_i = (A[i:][])` y $w_i = \begin{bmatrix}
        \| v_i \| \\
        0 \\
        \vdots \\
        0
      \end{bmatrix}$

2. Definimos $u = \frac{v-w}{||v-w||}$
3. Definimos $H_{u_i} = I - 2 u_i u_i^t$
4. Vemos $A' = H_{u_i}A = 
    \begin{bmatrix} 
    1 & 0 &  0 \\
    0 & H_{21} & H_{22} \\
    0 & H_{23} & H_{24} 
    \end{bmatrix}$

5. Repito el paso hasta tener $A'$ diagonal superior (la llamo $R$),
6. Defino $Q = H_r * ... * H_2 * H_1$
7. Luego $A = H_r * ... * H_2 * H_1 * R = QR$
