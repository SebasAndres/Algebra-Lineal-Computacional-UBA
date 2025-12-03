# Primer Parcial

## Índice

1. [Preeliminar](#1-preeliminar)
    1. [Propiedades de la dimensión](#11-propiedades-de-la-dimensión)
    2. [Linealmente independiente](#12-linealmente-independiente)
    3. [Propiedades matriciales](#13-propiedades-matriciales)
    4. [Subespacios](#14-subespacios)
    5. [Generadores y base](#15-generadores-y-base)
    6. [Cambio de base](#16-cambio-de-base)
    7. [Traza (propiedades)](#17-traza-propiedades)
    8. [Determinante (propiedades)](#18-determinante-propiedades)
2. [Transformaciones Lineales](#2-transformaciones-lineales)
3. [Proyecciones](#3-proyecciones)
    1. [Proyecciones en General](#31-proyecciones-en-general)
    2. [Proyección Ortogonal](#32-proyección-ortogonal)
        1. [Propiedades Geométricas](#321-propiedades-geométricas)
        2. [Teorema Fundamental de la Proyección (Descomposición Ortogonal)](#322-teorema-fundamental-de-la-proyección-descomposición-ortogonal)
        3. [Propiedades de la Matriz de Proyección](#323-propiedades-de-la-matriz-de-proyección)
        4. [Fórmulas de Construcción](#324-fórmulas-de-construcción)
        5. [Proyector a partir de una BON](#325-proyector-a-partir-de-una-bon)
        6. [Proyector Complementario](#326-proyector-complementario)
        7. [Relación de Subespacios](#327-relación-de-subespacios)
    3. [Ortonormalizar una matriz](#33-ortonormalizar-una-matriz)
4. [Normas](#4-normas)
    1. [Normas vectoriales](#41-normas-vectoriales)
    2. [Normas matriciales](#42-normas-matriciales)
    3. [Condición y propiedades](#43-condicion-y-propiedades)
5. [Descomposición $A=LU$](#5-descomposición-alu)
    1. [Existencia](#51-existencia)
    2. [Algoritmo](#52-algoritmo)
    3. [Descomposición $PA=LU$](#53-descomposición)
6. [Descomposición Cholesky $A=\hat{L}\hat{L}^T$](#6-descomposición-cholesky-a-hatlhatlt)
    1. [Existencia](#61-existencia)
    2. [Algoritmo](#62-algoritmo)
7. [Descomposición $A=QR$](#7-descomposición-aqr)
    1. [Existencia](#71-existencia)
    2. [Algoritmo HouseHolder](#71-algoritmo-householder)
    3. [Algoritmo GramSchmidt](#73-algoritmo-con-gram-schmidt)
8. [Teoremas y propiedades](#8-teoremas-y-propiedades)
    1. [Definida positiva](#81-definida-positiva)
    2. [Inversible](#82-inversible)
    3. [Diagonalizable](#83-diagonalizable)
    4. [Multiplicidades](#84-multiplicidades)
    5. [Matrices semejantes](#85-matrices-semejantes)

---

## 1. Preeliminar

### 1.1. Propiedades de la dimensión

- **Teorema de la dimensión**: $\dim(\operatorname{Im}(A)) + \dim(\ker(A)) = n$, donde $n$ es el número de columnas de $A$.
- $dim(A+B) = dim(A) + dim(B) - dim(A\cap B)$.
- $rk(A) = dim(Im(A)) = $"cantidad de ecuaciones li del sistema".
- $dim(ker(A)) = n - rk(A)$.
- $ker(A) = \{ x: Ax = 0 \}$.

- **Teorema**: $\mathbb{V}$ espacio de dim finita, $f: \mathbb{V} \rightarrow \text{W}$ una tl. Sup. $ker(f)$ tiene a $B$ como base y $B' = B \cup C$ una completación de una base de $\mathbb{V}$ con $C\cap B=\empty$, entonces $f(C)$ es una base de $\mathbb{W}$

### 1.2. Linealmente independiente
$\{ v_0, ..., v_i \}$ son vectores LI $\iff \text{el sistema homogeneo } [v_0 | ... | v_i]$ es compatible determinado.

### 1.3. Propiedades matriciales
- Cualquier matriz $A\in\mathbb{K}^{n\times n}$ puede escribirse como $A=\sum_i \sum_j a_{ij} E_{ij}$, siendo $E_{ij}$ la matríz definida por 
$$E_{ij} = \begin{cases} 
    1 & \text{si i=j} \\
    0 & \text{cc.}
\end{cases}$$

### 1.4. Subespacios
Un subespacio debe cumplir que:
- [1] $0 \in S$
- [2] $v\in S, w\in S \implies v+w \in S$
- [3] $\forall k\in \mathbb{K}, s \in S, kS \in S$

Más propiedades:
- $S$ y $T$ con generadores $\{ s_0, ..., s_i \}$ y $\{ t_0, ..., t_j \}$ entonces $S+T= \langle s_0, ..., s_i, t_0, ..., t_j \rangle$
- **Suma Directa:** $L+M$ es suma directa ($L \oplus M$) si y solo si $L \cap M = \{0\}$
- **Union de subespacios:** $S \cup T = \{ w : w\in S \lor w\in T\}$. En general no es un subespacio.
- $S^\perp = \{ v \mid \langle v, s \rangle = 0, \ \forall s \in S \}$ (subespacio ortogonal a $S$)
- $S = \langle s_0, ..., s_i \rangle, S \subset T \iff \forall s_j \in T$ 

### 1.5. Generadores y base
- Para obtener una base a partir de generadores puedo poner los vectores en filas, triangular la matríz y ver que, las filas que se vuelven 0 podemos sacarlas, y si se vuelven CL sacar 1.
- Para extender generados para formar una base puedo alinear los vectores como filas y completar las filas para que la matriz triangulada quede LI.

### 1.6. Cambio de base
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

### 1.7. Traza (propiedades)
- $tr(A+B) = tr(A) + tr(B)$
- $tr(\alpha A)=\alpha tr(A)$
- $tr(A^T) = tr(A)$
- $tr(AB) = tr(BA)$

### 1.8. Determinante (propiedades)
- Solo vale para matrices cuadradas.
- Si $A$ es triangular superior, $det(A) = \prod_{i=0}^n A_{ii}$
- $det(A) = det(A^T)$
- $det(kA) = k^n det(A)$, donde $n$ es el tamaño de la matriz cuadrada $A$ (es decir, $A \in \mathbb{K}^{n \times n}$)
- $det(AB) = det(A) det(B)$
- Si $A$ es invertible $\implies det(A^{-1}) = \frac{1}{det(A)}$
- $A$ inversible $\iff det(A) \neq 0$

## 2. Transformaciones Lineales
Sean $\mathbb{V}, \mathbb{W}$ dos $\mathbb{K}$ espacios vectoriales. Una funcion $f: \mathbb{V} \rightarrow \mathbb{W}$ se dice que es una transformacion si para todo $u,v \in \mathbb{V}$ y $\lambda \in \mathbb{K}$, vale que:
$$f(\lambda u + v) = \lambda f(u) + f(v)$$

> Teorema: Si $B$ es una base de $\mathbb{V}$ y $f: \mathbb{V} \rightarrow \mathbb{W}$, entonces $f(B)$ contiene una base de $Im(f)$.

## 3. Proyecciones

Este es el resumen completo de propiedades y teoremas clave sobre proyecciones que necesitas para el examen, completando y formalizando la sección que proporcionaste:

### 3.1. Proyecciones en General
Una matriz $P$ se llama **Proyector** si, al aplicarse dos veces, el resultado es el mismo que aplicarse una sola vez.

* **Idempotencia:** $P$ es una matriz de proyección $\iff P^2 = P$.
* **Autovalores:** Los autovalores de cualquier matriz de proyección $P$ solo pueden ser **0** o **1**.
* **Relación Núcleo/Imagen (Clave):** Para cualquier proyector $P$:
    $$\operatorname{Im}(P) = \operatorname{Nu}(I - P)$$
    (Los vectores en la imagen son los que quedan fijos: $Px=x$).

***

### 3.2. Proyección Ortogonal

Dado un subespacio $S \subset \mathbb{R}^n$, la **Proyección Ortogonal** sobre $S$, $P_S$, es la transformación lineal que asocia a cada vector $x$ su vector más cercano en $S$.

#### 3.2.1. Propiedades Geométricas
* **Vectores en S (Invarianza):** $P_S(s) = s \text{ } \forall s\in S$
* **Imagen en S:** $P_S(x) \in S \text{ } \forall x \in \mathbb{R}^n$
* **Vector Error:** El vector error $e = x − P_S(x)$ siempre es ortogonal a $S$:
    $$x−P_S(x) \in S^\perp \quad \text{o, equivalentemente, } P_S(x) \text{ es el vector de } S \text{ más cercano a } x$$

#### 3.2.2. Teorema Fundamental de la Proyección (Descomposición Ortogonal)
Todo vector $x \in \mathbb{R}^n$ puede descomponerse de forma **única** en la suma de dos componentes ortogonales:
$$x = P_S(x) + P_{S^\perp}(x)$$
Donde $P_S(x) \in S$ y $P_{S^\perp}(x) \in S^\perp$.

#### 3.2.3. Propiedades de la Matriz de Proyección
La matriz $P_S$ es una matriz de proyección si y solo si cumple:
* **Idempotencia:** $P_S^2 = P_S$
* **Simetría:** $P_S^T = P_S$ (Esta es la condición extra que la hace *ortogonal*).

#### 3.2.4. Fórmulas de Construcción

Sea $S = \text{col}(A)$, donde las columnas de $A$ forman una base para $S$.

| Caso | Matriz $A$ (Base de $S$) | Fórmula del Proyector $P_S$ |
| :--- | :--- | :--- |
| **Caso General** | Columnas de $A$ **L.I.** (no necesariamente ortonormales) | $$P_S = A(A^TA)^{-1}A^T$$ |
| **Caso Especial** | Columnas de $Q$ **ortonormales** ($Q^TQ = I$) | $$P_S = Q Q^T$$ |

**Nota Operativa:** En el examen de ALC, si $A$ no tiene columnas ortogonales, usualmente se le aplica **Gram-Schmidt** a las columnas de $A$ para obtener la matriz $Q$ de la base ortonormal, y así usar la fórmula simplificada $P_S = QQ^T$.

#### 3.2.5. Proyector a partir de una BON
Si tienes una **Base Ortonormal (B.O.N.)** de $S$ dada por $\{v_1, ..., v_r\}$, el proyector ortogonal sobre $S$ es la suma de proyectores de rango uno:
$$ P_S = \sum_{i=1}^r v_i v_i^T$$
Y su acción sobre un vector $x$ es:
$$ P_S(x) = P_S x = \sum_{i=1}^r (\underbrace{v_i^T x}_{\text{proy. sobre } v_i}) v_i$$

#### 3.2.6. Proyector Complementario
El proyector ortogonal sobre el complemento ortogonal $S^\perp$ se define como:
$$P_{S^\perp} = I - P_S$$

#### 3.2.7. Relación de Subespacios (Fundamental)
* $\operatorname{Im}(P_S) = S$
* $\operatorname{Nu}(P_S) = S^\perp$

### 3.3. Ortonormalizar una matriz

#### 3.3.1. Proceso (Gram-Schmidt)
Para ortonormalizar una matríz se puede usar Gram-Schmidt.

<u> Algoritmo: </u>

Sea $a_1,..., a_n$ las columnas de la matríz origen $A$.

1. Fijamos $v_1 := a_1$.  

3. Para cada $i \in [1..n]$:

- $v_{i+1} = a_{i+1} - \sum_{j=1}^i p_{v_j}(a_{i+1})$

4. Para cada $v_i$ calculado tomamos $v_i' = \frac{v_i}{||v_i||}$

Donde:
-  $p_a(b) = \frac{a^Tb}{a^Ta}a$

#### 3.3.2. Propiedades de matrices ortonormales
- $Q^TQ = I$
- $Im(A) = Im(Q) \text{ pero la transformación lineal es distinta, es decir, generalmente } Ax \neq Qx$
- $Q$ es más estable y tiene mejor numero de condicion.

## 4. Normas

### 4.1. Normas vectoriales
Una norma de un $\mathbb{K}$-espacio vectorial es una función $||.|| : V \rightarrow \mathbb{R}_{\geq 0}$ que cumple las siguientes propiedades:
1. $||av|| = |a| ||v|| \text{ para } a\in \mathbb{K} \land v \in V$
2. Si $||v|| = 0$, entonces $v=0$.
3. $||u+v|| \leq ||u|| + ||v|| \text{ para todo } u,v \in V$

#### Ejemplos en $\mathbb{K}^n$:
- Norma-1: $||v|| = |v_1| + ... + |v_n|$
- Norma-infinito: $||v||_{\infty} = \text{máx} \{ |v_1|, ..., |v_n| \}$
- Norma-p: $||v||_p = (|v_1|^p + ... + |v_n|^p)^{1/p}$

#### Desigualdad Cauchy-Schwartz:
$$|x^*y| \leq ||x|| * ||y||, \text{ para } x, y \in \mathbb{K}^n$$

#### Equivalencia entre normas:
Sean $||.||$ y $||.||_*$ dos normas en mismo $\mathbb{K}$-espacio vectorial $V$. Son equivalentes si $\exists c,C > 0$ tales que para todo $x\in V$:
$$c||x||_* \leq ||x|| \leq C ||x||_*$$

#### Convergencia de un vector a una norma
Sucede si $||v_n - v|| \rightarrow 0 \text{ cuando } n\rightarrow \infty$.

---

### 4.2. Normas matriciales
Dada $A \in \mathbb{K}^{n \times m}$ y un par de normas vectoriales $||.||_n$ y $||.||_m$ en $\mathbb{K}^n$ y $\mathbb{K}^m$, la norma inducida (o submultiplicativa) de matrices es:

$$
\|A\|_{n,m} = \max_{\substack{x \in \mathbb{K}^m \\ x \ne 0}} \frac{\|Ax\|_n}{\|x\|_m}
= \max_{\substack{x \in \mathbb{K}^m, ||x||_m = 1}} ||Ax||_n
$$

#### Ejemplos de normas matriciales:
- Norma-infinito: $||A||_{\infty} = \max_{1 \leq i \leq n} \sum_{j=1}^m |a_{ij}|$ (máx suma filas)
- Norma-1: $||A||_1 = \max_{1 \leq j \leq m} \sum_{i=1}^n |a_{ij}|$ (máx suma columnas)

#### Propiedades:
- $||Ax|| \leq ||A|| \cdot ||x||$
- $||AB|| \leq ||A|| \cdot ||B||$
- $||A|| = \sup_{v \neq 0} \left\{ \frac{||Av||}{||v||} \right\}$
- $Q$ unitaria, $||Qx|| = ||x||$
- $e \in \ker(A^T) \iff e \perp \operatorname{col}(A)$

---

### 4.3. Condición y propiedades

El **número de condición** de una matriz $A$ respecto a una norma vectorial es:

- $$ cond_*(A) = ||A||_* \cdot ||A^{-1}||_* $$

Propiedades importantes del número de condición:

- $$cond_*(A) \geq \sup_{H \text{ singular}} \left\{ \frac{||A||_* } {||A - H||_*}\right\}$$
- $$cond_*(A) \leq \inf_{H \text{ singular}} \left\{ \frac{||A-H||_* } {||A||_* }\right\}$$
- $$\frac{||\tilde{x}-x||}{||x||} \leq cond_*(A) \cdot \frac{||\tilde{b}-b||}{||b||}$$

Un número de condición grande implica que el sistema es mal condicionado y sensible a errores de redondeo.

---

## 5. Descomposición $A=LU$
$$A = LU$$
- $L$: triangular inferior (con la diagonal principal).
- $U$: triangular superior.

### 5.1 Existencia
Para matrices cuadradas sii durante la eliminación gaussiana, todos los pivotes deben ser $\neq 0$ (sin intercambiar filas), esto sucede si $det(A) \neq 0$.

### 5.2 Algoritmo
1. Triangular la matríz original con operaciones elementales entre filas hasta que quede triangular superior.
2. La matriz triangular superior resultante es $U$, a $L$ la construimos agregando en cada columna en el valor $L_{ij}$ el inverso en signo del multiplicador usado para triangular la celda $ij$.

### 5.3 Descomposición $PA=LU$
Esta descomposición existe para cualquier matríz cuadrada.
Es equivalente al algoritmo de descomposición $LU$, pero se aplica previamente una permutación expresada con el producto con una matríz $P$.

## 6. Descomposición Cholesky $A=\hat{L}\hat{L}^T$
### 6.1 Existencia
La **descomposición de Cholesky** $A = \hat{L}\hat{L}^T$ existe **si y sólo si $A$ es simétrica y definida positiva**.

- $A$ **simétrica**: $A^T = A$
- $A$ **definida positiva**: $x^T A x > 0$ para todo $x \neq 0$

Es decir, para cualquier matriz cuadrada real $A$ tal que $a_{ij} = a_{ji}$ y $x^T A x > 0$ para todo vector $x\neq 0$, se puede escribir $A$ como el producto de una triangular inferior $\hat{L}$ y su transpuesta.
  
- La factorización es única cuando $A$ es definida positiva y simétrica.

### 6.2 Algoritmo
1. Calculo $A=LU$
2. Defino D de la siguiente forma:
    $$ D_{ij} = \begin{cases}
    L_{ii} & \text{si i=j} \\
    0 & \text{cc.}
    \end{cases}$$
3. Defino $D_1$ de tal forma que $D_1 * D_1 = D$, o sea:
    $(D_1)_{ii} = \sqrt{D_{ii}}$
4. Defino $\hat{L} = L*D_1$
5. Luego $\hat{L}^T$ es simplemente tomar la transpuesta de la calculada en el paso anterior y tengo $A=\hat{L}\hat{L}^T$

## 7. Descomposición $A=QR$

### 7.1 Existencia
La descomposición $A=QR$ existe para cualquier matriz $A\in \mathbb{R}^{m \times n}$ con $m \geq n$ y rango completo de columnas (columnas linealmente independientes). Si $A$ no tiene rango completo, aún se puede realizar la descomposición pero algunas filas de $R$ serán nulas.
- **Para matrices cuadradas ($m = n$):** Siempre se puede hacer $A = QR$.
- **Para matrices rectangulares ($m > n$):** $Q$ será $m \times m$ ortogonal/unitaria y $R$ será $m \times n$ (puede tener ceros al final).

Podemos calcularla usando Gram-Schmidt o Householder.

### 7.1 Algoritmo HouseHolder

1. Tenemos $A$ y tomamos $v_i$ como la iésima columna de A `v_i = (A[i:][])` y $w_i = \begin{bmatrix}
        \| v_i \| \\
        0 \\
        \vdots \\
        0
      \end{bmatrix}$

2. Definimos $u = \frac{v-w}{||v-w||}$
3. Definimos $H_{u_i} = I - 2 u_i u_i^T$
4. Vemos $A' = H_{u_i}A = 
    \begin{bmatrix} 
    1 & 0 &  0 \\
    0 & H_{21} & H_{22} \\
    0 & H_{23} & H_{24} 
    \end{bmatrix}$

5. Repito el paso hasta tener $A'$ diagonal superior (la llamo $R$),
6. Defino $Q = H_r * ... * H_2 * H_1$
7. Luego $A = H_r * ... * H_2 * H_1 * R = QR$

### 7.3 Algoritmo con Gram-Schmidt
1. Calculo una b.o.n. a partir de $A$ y obtengo $Q$ unitaria de la forma:
    $$Q=[ q_1 | ... | q_n]$$
    Con $q_i = \frac{v_i}{||v_i||}, \text{ siendo } v_i \text{ vectores ortogonales}$.

2. Defino la matríz triangular superior $R$ como $R=Q^T A$
    $$R = \begin{pmatrix}
    \mathbf{q}_1^T \mathbf{a}_1 & \mathbf{q}_1^T \mathbf{a}_2 & \mathbf{q}_1^T \mathbf{a}_3 & \cdots \\
    0 & \mathbf{q}_2^T \mathbf{a}_2 & \mathbf{q}_2^T \mathbf{a}_3 & \cdots \\
    0 & 0 & \mathbf{q}_3^T \mathbf{a}_3 & \cdots \\
    \vdots & \vdots & \vdots & \ddots
    \end{pmatrix}$$
    Donde los elementos en la diagonal son $r_{kk} = \|\mathbf{v}_k\|$.

## 8. Teoremas y propiedades
- $ A \text{ simetrica} \rightarrow det(A) = \prod_i \lambda_i$

### 8.1. Definida positiva
- $A \text{ es definida positiva } \iff A=LU \text{ con } U_{ii} > 0 \text{ para } \forall i \in \{1 .. n\}$

### 8.2. Inversible
- $A \text{ inversible} \iff det(A) = 0 \iff \exists v \neq 0 | Av = 0$
- $ A \text { no es inversible } \iff \lambda=0 \text{ es autovalor}$
- $A\in K^{N\times N} \text{ invertible } \iff \text{Los vectores columna de A forman una base}$

### 8.3. Diagonalizable
- $ A\in K^{N\times N} \text { es diagonalizable } \iff \text{ Existen N autovalores distintos }$
- $ A\in K^{N\times N} \text { es diagonalizable } \iff \text{ Existen N autovectores LI (vectores columna de A)}$
- $ A\in K^{N\times N} \text { es diagonalizable } \iff \text{ Para todo } \lambda_i \text{ autovalor de A, vale que } mg_A(\lambda_i) = ma_A(\lambda_i)$
- $ A\in K^{N\times N} \text { es diagonalizable } \iff A = PDP^{-1}, D \text{ diagonal }$
- $ A\in K^{N\times N} \text { es diagonalizable } \iff A^m = PD^mP^{-1}$
- $ A \text{ es diagonalizable si es semejante a una matriz diagonal}$

### 8.4. Multiplicidades
- $ mg_a(\lambda) = dim(E_\lambda)$
- $ ma_a(\lambda) = \text{ "multiplicidad de lambda como raíz en } x(\lambda) \text{" }$
- $1 \leq mg_a(\lambda) \leq ma_a(\lambda)$

### 8.5. Matrices semejantes
- $A, B \in K^{N\times N} \text{ son semejantes } \iff \exists c \in K^{N \times N} \text{ tal que } A = CBC^{-1}$ 
- $ \text{A y B son semejantes si } A= CBC^{-1}$
- $P(x) = [ P ] x = \sum (x * e_i) * e_i, \text{ con } e_i \text{ los vectores columnas de la matriz ortogonal de P} $
