# Primer Parcial

## Preeliminar

### Sistemas
- [1] **Compatible determinado**: Tiene solución única. Todos los elementos de la diagonal de la matriz triangulada son distintos a 0.
- [2] **Compatible indeterminado**: Tiene infinitas soluciones. Si al menos un elemento de la diagonal en la matriz triangulada es 0 y el termino independiente de esa fila es 0 (en los homogeneos valdría siempre).
- [3] **Indeterminado**: No tiene solución. Si al menos un elemento de la diagonal en la matriz triangulada es 0 y el termino independiente de esa fila NO es 0.

### Propiedades
$$dim(S) = n - rk(A), \text{ con } S = \{x: Ax= 0\}$$
> Teo: $\mathbb{V}$ espacio de dim finita, $f: \mathbb{V} \rightarrow \text{W}$ una tl.
> Sup $ker(f)$ tiene a $B$ como base y $B' = B \cup C$ una completación de una base de $\mathbb{V}$ con $C\cap B=\empty$, entonces $f(C)$ es una base de $\mathbb{W}$

### Linealmente independiente
Veo que un cjt de vectores $\{ v_0, ..., v_i \}$ es LI viendo que el sistema homogeneo correspondiente a agregarlos como vectores columna es compatible determinado. 

### Propiedades matriciales
- $A \text{ es inversible} \iff [Ax = 0 \iff x=0] \iff det(A) \neq 0$
- $S^⊥ =^{def}$?
- Cualquier matriz $A\in\mathbb{K}^{n\times n}$ puede escribirse como $A=\sum_i \sum_j a_{ij} E_{ij}$, siendo $E_{ij}$ la matríz definida por 
$$E_{ij} = \begin{cases} 
    1 & \text{si i=j} \\
    0 & \text{cc.}
\end{cases}$$
- $Q$ matriz cuadrada $\rightarrow Im(Q) = Nu(I-Q)$

### Subespacios
Un subespacio debe cumplir que:
- [1] $0 \in S$
- [2] $v\in S, w\in S \rightarrow v+w \in S$
- [3] $\forall k\in \mathbb{K}, s \in S, kS \in S$

Más propiedades:
- $S$ y $T$ dos subespacios están en suma directa si $S\cap T = \{ 0 \}$.
- $S$ y $T$ con generadores $\{ s_0, ..., s_i \}$ y $\{ t_0, ..., t_j \}$ entonces $S+T= \langle s_0, ..., s_i, t_0, ..., t_j \rangle$
- $S = \langle s_0, ..., s_i \rangle, S \subset T \iff \forall s_j, s_j \in T$ 
- $S \cup T = \{ w : w\in S \lor w\in T\}$. En general no es un subespacio.

### Generadores y base
- Para obtener una base a partir de generadores puedo poner los vectores en filas, triangular la matríz y ver que, las filas que se vuelven 0 podemos sacarlas, y si se vuelven CL sacar 1.
- Para extender generados para formar una base puedo alinear los vectores como filas y completar las filas para que la matriz triangulada quede LI.

### Cambio de base
- **Coordenadas**: Sea $V = \langle v_0, ..., v_m \rangle$. Se llama coordenadas de $v = \sum_i^m \alpha_i v_i \in V$ en base $B=\{v_1, ..., v_m\}$ al vector $[w]_B = (\alpha_0, ..., \alpha_m)$. 
- **Matriz cambio de base**: Si quiero escribir a $w$ en otra base $B'=\{  z_1, ..., z_n \}$ y tengo las coordenadas en la base $B$. Vale que:
  - En base canónica: $w=[v_1 | ... | v_m]\begin{pmatrix} \alpha_1 \\ ... \\ \alpha_n \end{pmatrix} = C_{BE}\begin{pmatrix} \alpha_1 \\ ... \\ \alpha_n \end{pmatrix}$
  - $w = [z_1 | ... | z_n]\begin{pmatrix} \beta_1 \\ ... \\ \beta_n \end{pmatrix} = C_{B'E}\begin{pmatrix} \beta_1 \\ ... \\ \beta_n \end{pmatrix}$

Luego, juntando:
$$w=C_{BE}[w]_B = C_{B'E}[w]_{B'}$$
Entonces
$[w]_{B'} = C_{B'E}^{-1}C_{BE}[w]_B$

Donde $C_{BB'}=C_{B'E}^{-1}C_{BE}$ es la "matriz de cambio de base".

> Observación:
> $C_{AB}$ es la nomenclatura para la matriz de cambio de base que recibe en $A$ y devuelve en $B$.

### Transposición (propiedades)
- $(A^t)^t = A$
- $(\alpha A)^t = \alpha A^t$
- $(AB)^t = B^tA^t$

### Traza (propiedades)
- $tr(A+B) = tr(A) + tr(B)$
- $tr(\alpha A)=\alpha A$
- $tr(A^t) = tr(A)$
- $tr(AB) = tr(BA)$

### Determinante (propiedades)
- Si $A$ es triangular superior, $det(A) = \prod_{i=0}^n A_{ii}$
- $det(A) = det(A^t)$
- $det(kA) = k^n det(A)$
- $det(AB) = det(A) det(B)$
- Si $A$ es invertible, $det(A^{-1}) = \frac{1}{det(A)}$
- $A$ inversible $\iff det(A) \neq 0$

### Matríz adjunta (transpuesta y conjugada) 
- $A^* = \hat{A^t}$

## Transformaciones Lineales
Sean $\mathbb{V}, \mathbb{W}$ dos $\mathbb{K}$ espacios vectoriales. Una funcion $f: \mathbb{V} \rightarrow \mathbb{W}$ se dice que es una transformacion si para todo $u,v \in \mathbb{V}$ y $\lambda \in \mathbb{K}$, vale que:
$$f(\lambda u + v) = \lambda f(u) + f(v)$$

> Teorema: Si $B$ es una base de $\mathbb{V}$ y $f: \mathbb{V} \rightarrow \mathbb{W}$, entonces $f(B)$ contiene una base de $im(f)$.

> Monomorfismo:

> Isomorfismo:


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
$ P_S(x) = P * x = (\sum_{i=1}^r v_i v_i^t) * x $

- $Q \text{ matriz de proyeccion } \rightarrow Q^2 = Q$
- $Q \text{ matriz de proyeccion ortogonal } \rightarrow Q^2 = Q \land Q^t = Q$

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
- $Q$ es más estable y tiene mejor numerVo de condicion.

## Normas vectoriales
Una norma de un $\mathbb{K}$-espacio vectorial es una función $||.|| : V \rightarrow \mathbb{R}_{\geq 0}$ que cumple las siguientes propiedades:
1. $||av|| = |a| ||v|| \text{ para } a\in \mathbb{K} \land v \in V$
2. Si $||v|| = 0$, entonces $v=0$.
3. $||u+v|| \leq ||u|| + ||v|| \text{ para todo } u,v \in V$

#### Normas vectoriales (en $\mathbb{K}^n$) comunes:
- Norma-1: $||v|| = |v_1| + ... + |v_n|$
- Norma-infinito: $||v||_{\infty} = \text{máx} \{ |v_1|, ..., |v_n| \}$
- Norma-p: $||v||_p = (|v_1|^p + ... + |v_n|^p)^{1/p}$

#### Desigualdad Cauchy-Schwartz:
$$|x^*y| \leq ||x|| * ||y||, \text{ para } x, y \in \mathbb{K}^n$$

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
- $$e \in Ker(A^T) \iff e \perp c(A)$$
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
    $(D_1)_{ii} = \sqrt{D_{ii}}$
4. Defino $\hat{L} = L*D_1$
5. Luego $\hat{L}^t$ es simplemente tomar la transpuesta de la calculada en el paso anterior y tengo $A=\hat{L}\hat{L}^t$

## Descomposición $A=QR$

Podemos calcularla usando Gram-Schmidt o Householder.

<u>Algoritmo con Householder:</u>
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

<u>Algoritmo con Gram-Schmidt:</u>
1. Calculo una b.o.n. a partir de $A$ y obtengo $Q$ unitaria de la forma:
    $$Q=[ q_1 | ... | q_n]$$
    Con $q_i = \frac{v_i}{||v_i||}, \text{ siendo } v_i \text{ vectores ortogonales}$.

2. Defino la matríz triangular superior $R$ como $R=Q^t A$
    $$R = \begin{pmatrix}
    \mathbf{q}_1^T \mathbf{a}_1 & \mathbf{q}_1^T \mathbf{a}_2 & \mathbf{q}_1^T \mathbf{a}_3 & \cdots \\
    0 & \mathbf{q}_2^T \mathbf{a}_2 & \mathbf{q}_2^T \mathbf{a}_3 & \cdots \\
    0 & 0 & \mathbf{q}_3^T \mathbf{a}_3 & \cdots \\
    \vdots & \vdots & \vdots & \ddots
    \end{pmatrix}$$
    Donde los elementos en la diagonal son $r_{kk} = \|\mathbf{v}_k\|$.

## Teoremas y propiedades
- $$ A \text{ simetrica} \rightarrow det(A) = \prod_i \lambda_i$$

### Condición
- $$ cond_k(A) = || A ||_k * || A^{-1} ||_k $$
- $$ cond_k(A) \geq \sup_{H \text{ singular}} \{ \frac{|| A ||_k } {|| A - H ||_k }\}$$
- $$ cond_k(A) \leq \inf_{H \text{ singular}} \{ \frac{|| A - H||_k } {|| A ||_k }\}$$
- $$\frac{||\tilde{x}-x||}{||x||} \leq cond(A) \frac{||\tilde{b}-b||}{||b||}$$

### Definida positiva
- $$A \text{ es definida positiva } \iff A=LU \text{ con } U_{ii} > 0 \text{ para } \forall i \in \{1 .. n\}$$

### Inversible
- $$A \text{ inversible} \iff det(A) = 0 \iff \exists v \neq 0 | Av = 0$$
- $$ A \text { no es inversible } \iff \lambda=0 \text{ es autovalor}$$

### Diagonalizable
- $$A\in K^{N\times N} \text{ invertible } \iff \text{Los vectores columna de A forman una base} $$
- $$ A\in K^{N\times N} \text { es diagonalizable } \iff \text{ Existen N autovectores LI (vectores columna de A)}$$
- $$ A\in K^{N\times N} \text { es diagonalizable } \iff \text{ Para todo } \lambda_i \text{ autovalor de A, vale que } mg_A(\lambda_i) = ma_A(\lambda_i)$$
- $$ A\in K^{N\times N} \text { es diagonalizable } \iff A = PDP^{-1}, D \text{ diagonal }$$
- $$ A\in K^{N\times N} \text { es diagonalizable } \iff A^m = PD^mP^{-1}$$
- $$ A \text{ es diagonalizable si es semejante a una matriz diagonal}$$

### Multiplicidades
- $$ mg_a(\lambda) = dim(E_\lambda)$$
- $$ ma_a(\lambda) = \text{ "multiplicidad de lambda como raíz en } x(\lambda) \text{" }$$
- $$1 \leq mg_a(\lambda) \leq ma_a(\lambda)$$

### Matrices semejantes
- $A, B \in K^{N\times N} \text{ son semejantes } \iff \exists c \in K^{N \times N} \text{ tal que } A = CBC^{-1}$ 
- $$ \text{A y B son semejantes si } A= CBC^{-1}$$
- $$P(x) = [ P ] x = \sum (x * e_i) * e_i, \text{ con } e_i \text{ los vectores columnas de la matriz ortogonal de P}$$
