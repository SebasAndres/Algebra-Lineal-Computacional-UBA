# Primer Parcial

## Índice

1. [Preeliminar](#1-preeliminar)
    1. [Propiedades de la dimensión](#11-propiedades-de-la-dimensión)
    2. [Inversible (y vectores linealmente independientes)](#12-inversible-y-vectores-li-inealmente-independientes)
    3. [Propiedades matriciales](#13-propiedades-matriciales)
    4. [Subespacios](#14-subespacios)
    5. [Generadores y base](#15-generadores-y-base)
2. [Transformaciones Lineales y Cambio de Base](#2-transformaciones-lineales-y-cambio-de-base)
    1. [Definición y Propiedades Básicas](#21-definición-y-propiedades-básicas)
    2. [Núcleo e Imagen](#22-núcleo-e-imagen)
    3. [Representación Matricial de una T.L.](#23-representación-matricial-de-una-tl)
    4. [Matriz de Cambio de Base (Vectores)](#24-matriz-de-cambio-de-base-vectores)
    5. [Cambio de Base en Transformaciones Lineales](#25-cambio-de-base-en-transformaciones-lineales)
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
        1. [Proceso (Gram-Schmidt)](#331-proceso-gram-schmidt)
        2. [Propiedades de matrices ortonormales](#332-propiedades-de-matrices-ortonormales)
4. [Normas](#4-normas)
    1. [Normas vectoriales](#41-normas-vectoriales)
    2. [Normas matriciales](#42-normas-matriciales)
    3. [Condición y propiedades](#43-condición-y-propiedades)
5. [Descomposición $A=LU$](#5-descomposición-alu)
    1. [Existencia](#51-existencia)
    2. [Algoritmo](#52-algoritmo)
    3. [Descomposición $PA=LU$](#53-descomposición-palu)
6. [Descomposición Cholesky $A=\hat{L}\hat{L}^T$](#6-descomposición-cholesky-ahatlhatlt)
    1. [Existencia](#61-existencia)
    2. [Algoritmo](#62-algoritmo)
7. [Descomposición $A=QR$](#7-descomposición-aqr)
    1. [Existencia](#71-existencia)
    2. [Algoritmo Householder](#72-algoritmo-householder)
    3. [Algoritmo con Gram-Schmidt](#73-algoritmo-con-gram-schmidt)
8. [Teoremas y propiedades](#8-teoremas-y-propiedades)
    1. [Definida positiva](#81-definida-positiva)
    2. [Diagonalizable](#83-diagonalizable)
    3. [Multiplicidades](#84-multiplicidades)
    4. [Matrices semejantes](#85-matrices-semejantes)

---

## 1. Preeliminar

### 1.1. Propiedades de la dimensión

- **Teorema de la dimensión**: $\dim(\operatorname{Im}(A)) + \dim(\ker(A)) = n$, donde $n$ es el número de columnas de $A$.
- $dim(A+B) = dim(A) + dim(B) - dim(A\cap B)$.
- $rk(A) = dim(Im(A)) = $"cantidad de ecuaciones li del sistema".
- $dim(ker(A)) = n - rk(A)$.
- $ker(A) = \{ x: Ax = 0 \}$.

- **Teorema**: $\mathbb{V}$ espacio de dim finita, $f: \mathbb{V} \rightarrow \text{W}$ una tl. Sup. $ker(f)$ tiene a $B$ como base y $B' = B \cup C$ una completación de una base de $\mathbb{V}$ con $C\cap B=\empty$, entonces $f(C)$ es una base de $\mathbb{W}$

### 1.2. Inversible (y vectores LI inealmente independientes)

**A singular**
- $A \text { no es inversible } \\
\iff \lambda=0 \text{ es autovalor} \\
\iff det(A)=0 \\
\iff dim(Nu(A)) > 0 \\
\iff \text {Los vectores columna de A son LD}$

**A inversible**
- $A \text{ inversible} \\
\iff det(A) \neq 0 \\
\iff \text{Los vectores columna de A son LI} \\
$

**Caso particular**
- $A\in K^{N\times N} \text{ invertible } \iff \text{Los vectores columna de A forman una base}$
- $A \text{ estrictamente diagonal dominante (la diafonal es mayor en mulo que el resto de valores) } \implies \text{ A invertible}$

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
- **Suma de subespacios:** Sean $S$ y $T$ subespacios con generadores $\{ s_0, ..., s_i \}$ y $\{ t_0, ..., t_j \}$. Entonces, la suma $S+T$ es el subespacio generado por la unión de sus generadores:
  $$S+T = \langle s_0, ..., s_i, t_0, ..., t_j \rangle$$

- **Intersección de subespacios:** Dados $S$ y $T$ generados por los conjuntos $\{ s_0, ..., s_i \}$ y $\{ t_0, ..., t_j \}$ respectivamente, la intersección $S \cap T$ es el conjunto de vectores que pueden escribirse tanto como combinación lineal de los generadores de $S$ como de los de $T$:
  $$S \cap T = \{ v : v \in S \text{ y } v \in T \}$$
  Para encontrar el conjunto, se resuelve el sistema $S\alpha = T\beta$ con variables libres $\alpha$ y $\beta$, es decir, se buscan los vectores que pueden escribirse en ambos sistemas de generadores simultáneamente. Alternativamente, $S\cap T$ es el subespacio generado por todos los vectores comunes a ambos subespacios.

  **Forma matricial:** Si representamos las matrices $S$ (con columnas los generadores de $S$) y $T$ (con columnas los generadores de $T$), la intersección $S \cap T$ puede obtenerse encontrando todas las soluciones del sistema homogéneo $S\alpha - T\beta = 0$, es decir:
  $$
  \begin{pmatrix} S & -T \end{pmatrix} \begin{pmatrix} \alpha \\ \beta \end{pmatrix} = 0
  $$
  Luego, los vectores de la forma $S\alpha$ (o $T\beta$) que satisfacen esta ecuación constituyen la intersección $S \cap T$. Así, se puede obtener una base del subespacio intersección de manera algorítmica usando métodos matriciales.

- **Suma Directa:** La suma $L+M$ es directa ($L \oplus M$) si, y sólo si, $L \cap M = \{0\}$.

- **Unión de subespacios:** $S \cup T = \{ w : w \in S \lor w \in T \}$, aunque en general la unión de subespacios no es un subespacio salvo casos triviales.

- **Subespacio contenido:** Si $S = \langle s_0, ..., s_i \rangle$, entonces $S \subset T$ si y sólo si $\forall s_j\, :\, s_j \in T$.

- Se verifica que $A \cap B = \{0\}$ comprobando que no es posible expresar los generadores de un subespacio como combinaciones lineales de los generadores del otro.

- **Subespacios ortogonales:** Dos subespacios $S$ y $T$ son ortogonales si todo generador de $S$ es perpendicular a todo generador de $T$ (y viceversa).

- El subespacio ortogonal de $S$ se denota:
  $$S^\perp = \{ v \mid \langle v, s \rangle = 0,\, \forall s \in S \}$$

### 1.5. Generadores y base
- Para obtener una base a partir de generadores puedo poner los vectores en filas, triangular la matríz y ver que, las filas que se vuelven 0 podemos sacarlas, y si se vuelven CL sacar 1.
- Para extender generados para formar una base puedo alinear los vectores como filas y completar las filas para que la matriz triangulada quede LI.

## 2. Transformaciones Lineales y Cambio de Base

Esta sección unifica el concepto abstracto de función lineal con su implementación práctica a través de matrices y coordenadas.

### 2.1. Definición y Propiedades Básicas

Sean $\mathbb{V}$ y $\mathbb{W}$ dos $\mathbb{K}$-espacios vectoriales. Una función $f: \mathbb{V} \rightarrow \mathbb{W}$ es una **Transformación Lineal** si conserva la estructura de espacio vectorial, es decir:

$$f(\alpha u + v) = \alpha f(u) + f(v), \quad \forall u, v \in \mathbb{V}, \forall \alpha \in \mathbb{K}$$

#### Propiedades inmediatas:

1.  **El cero va al cero:** $f(0_\mathbb{V}) = 0_\mathbb{W}$.

2.  **Inversos aditivos:** $f(-v) = -f(v)$.

3.  **Preservación de combinaciones lineales:** $f(\sum \alpha_i v_i) = \sum \alpha_i f(v_i)$.

### 2.2. Núcleo e Imagen

Son los dos subespacios fundamentales asociados a una T.L.:

* **Núcleo (Kernel):** Es el conjunto de vectores del dominio que se transforman en el cero.

    $$\ker(f) = \{ v \in \mathbb{V} : f(v) = 0_\mathbb{W} \}$$

    * **Propiedad:** $f$ es inyectiva (monomorfismo) $\iff \ker(f) = \{0\}$.

* **Imagen:** Es el conjunto de vectores del codominio que son "alcanzados" por la función.

    $$\operatorname{Im}(f) = \{ w \in \mathbb{W} : \exists v \in \mathbb{V}, f(v) = w \}$$

    * **Propiedad:** $f$ es sobreyectiva (epimorfismo) $\iff \operatorname{Im}(f) = \mathbb{W}$.

#### Teorema de la Dimensión (Teorema del Rango-Nulidad)

Si $\mathbb{V}$ es de dimensión finita:

$$\dim(\mathbb{V}) = \dim(\ker(f)) + \dim(\operatorname{Im}(f))$$

---

### 2.3. Representación Matricial de una T.L.

Toda transformación lineal entre espacios de dimensión finita puede representarse mediante una matriz. La forma de esta matriz depende de las bases elegidas.

Sean $B = \{v_1, \dots, v_n\}$ una base de $\mathbb{V}$ y $B' = \{w_1, \dots, w_m\}$ una base de $\mathbb{W}$.

La **matriz asociada a $f$ en las bases $B$ y $B'$**, denotada como $[f]_{BB'}$ (o $M(f)_{BB'}$), se construye colocando en sus columnas las coordenadas de las imágenes de los vectores de la base $B$ escritas en la base $B'$.

$$[f]_{BB'} = \begin{pmatrix} | & & | \\ [f(v_1)]_{B'} & \cdots & [f(v_n)]_{B'} \\ | & & | \end{pmatrix}$$

#### Relación fundamental:

Para transformar un vector $v$, multiplicamos sus coordenadas por la matriz:

$$[f(v)]_{BB'} = [f]_{BB'} \cdot [v]_B$$

---

### 2.4. Matriz de Cambio de Base (Vectores)

Si queremos cambiar las coordenadas de un vector de una base a otra dentro del mismo espacio $\mathbb{V}$:

Sean $B = \{v_1, \dots, v_n\}$ y $B' = \{u_1, \dots, u_n\}$ dos bases de $\mathbb{V}$.

La **Matriz de Cambio de Base de $B$ a $B'$**, denotada como $C_{BB'}$, cumple que:

$$[v]_{B'} = C_{BB'} \cdot [v]_B$$

#### Propiedades de $C_{BB'}$:

1.  **Construcción:** Las columnas de $C_{BB'}$ son las coordenadas de los vectores de la base "vieja" $B$ escritos en función de la base "nueva" $B'$.

2.  **Inversibilidad:** Toda matriz de cambio de base es inversible.

3.  **Inversa:** $C_{B'B} = (C_{BB'})^{-1}$.

4.  **Composición:** $C_{BB''} = C_{B'B''} \cdot C_{BB'}$.

> **Nota práctica:** Es fácil construir la matriz de cambio de base $C_{BE}$ (de una base $B$ cualquiera a la canónica $E$), simplemente poniendo los vectores de $B$ como columnas. Para volver (de canónica a $B$), calculamos la inversa: $C_{EB} = (C_{BE})^{-1}$. Si $B$ es ortonogonal, $ C_{EB} = C_{BE}^{-1} = C_{BE}^T$. 

---

### 2.5. Cambio de Base en Transformaciones Lineales

$$[f]_{BE} = [f]_{EE} C_{BE}$$

- Donde $f_{EE}$ está definida por las ecuaciones dadas de $f(X)=(\phi_1(X), ..., \phi_n(X))$
> **Conclusión Clave:** Dos matrices representan la misma transformación lineal en distintas bases si y solo si son **semejantes**. Esto implica que comparten propiedades intrínsecas como el determinante, la traza y los autovalores.

## 3. Proyecciones

Este es el resumen completo de propiedades y teoremas clave sobre proyecciones que necesitas para el examen, completando y formalizando la sección que proporcionaste:

### 3.1. Proyecciones en General
Una matriz $P$ se llama **Proyector** si, al aplicarse dos veces, el resultado es el mismo que aplicarse una sola vez.

* **Idempotencia:** $P$ es una matriz de proyección $\iff P^2 = P$.
* **Autovalores:** Los autovalores de cualquier matriz de proyección $P$ solo pueden ser **0** o **1**.
* **Relación Núcleo/Imagen (Clave):** Para cualquier proyector $P$:
    $$\operatorname{Im}(P) = \operatorname{Nu}(I - P)$$
    $$\operatorname{Im}(I-P) = \operatorname{Nu}(- P) = \operatorname{Nu}(P)$$
    (Los vectores en la imagen son los que quedan fijos: $Px=x$).

#### Armado de proyector en general:
Supongamos tengo un subespacio $M = [v_1 | ... | v_n]$
$$P M = R = [p(v_1) | ... | p(v_n)]$$
$$P = RM^{-1}$$

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
* $Nu(P) \perp Im(P)$

#### 3.2.4. Fórmulas de Construcción

Sea $S = \text{col}(A)$, donde las columnas de $A$ forman una base para $S$.

| Caso | Matriz $A$ (Base de $S$) | Fórmula del Proyector $P_S$ |
| :--- | :--- | :--- |
| **Caso General** | Columnas de $A$ **L.I.** (no necesariamente ortonormales) | $$P_S = RS^{-1}$$ |
| **Caso Especial** | Columnas de $Q$ **ortonormales** ($Q^TQ = I$) | $$P_S = Q Q^T$$ |

**Nota Operativa:** Si $A$ no tiene columnas ortogonales, usualmente se le aplica **Gram-Schmidt** a las columnas de $A$ para obtener la matriz $Q$ de la base ortonormal, y así usar la fórmula simplificada $P_S = QQ^T$.

#### 3.2.5. Proyector a partir de una BON
Si tienes una **Base Ortonormal (B.O.N.)** de $S$ dada por $\{v_1, ..., v_r\}$, el proyector ortogonal sobre $S$ es la suma de proyectores de rango uno:
$$ P_S = \sum_{i=1}^r v_i v_i^T$$
Y su acción sobre un vector $x$ es:
$$ P_S(x) = P_S x = \sum_{i=1}^r (\underbrace{v_i^T x}_{\text{proy. sobre } v_i}) v_i$$

#### 3.2.6. Proyector Complementario
El proyector ortogonal sobre el complemento ortogonal $S^\perp$ se define como:
$$P_{S^\perp} = I - P_S$$

#### 3.2.7. Relación de Subespacios (matriz de proyeccion ortonormal)
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
- $Q^T = Q^{-1}$
- $Q$ unitaria, $||Qx|| = ||x||$

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
---

### 4.3. Condición y propiedades

El **número de condición $(k)$** de una matriz $A$ respecto a una norma vectorial es:

- $ k_*(A) = ||A||_* \cdot ||A^{-1}||_* $

Propiedades importantes del número de condición:
- $k_*(A) \geq \sup_{\substack{H \\ \text{singular}}} \left\{ \frac{||A||_* }{||A - H||_*} \right\}$
- $k_*(A) \leq \inf_{\substack{H \\ \text{singular}}} \left\{ \frac{||A - H||_*}{||A||_*} \right\}$

**Observación**
- $k \in [1, +\inf)$
- Un número de condición grande implica que el sistema es mal condicionado y sensible a errores de redondeo.

#### Error relativo
Se puede acotar el error relativo de computo de $x$ con:
$$\frac{\lVert x - \tilde{x} \rVert}{\lVert x \rVert} \le \kappa(A) \frac{\lVert b - \tilde{b} \rVert}{\lVert b \rVert}$$
---

## 5. Descomposición $A=LU$
$$A = LU$$
- $L$: triangular inferior (con la diagonal principal).
- $U$: triangular superior.

### 5.1 Existencia
Existe LU para matrices cuadradas 
$\iff \text{durante la eliminación gaussiana, todos los pivotes deben ser distintos a cero (sin intercambiar filas)} \\
\iff det(A_k) \neq 0 \text{ para Ak submatriz principal hasta n-1 (todos sus menores principales son distintos a 0)}$

**Unicidad**
La factorización es única si A es invertible.

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
  
**Unicidad** 
La factorización es única cuando $A$ es definida positiva y simétrica.

### 6.2 Algoritmo
1. Calculo $A=LU$
2. Defino D de la siguiente forma:
$$ D_{ij} = \begin{cases}
    L_{ii} & \text{si i=j} \\
    0 & \text{cc.}
\end{cases}
$$
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

**Propiedad**
$Im(A) = Im(Q) \text{ pero la transformación lineal es distinta, es decir, generalmente } Ax \neq Qx$

**Unicidad** 
QR es única si y solo si A tiene columnas LI o tambien $R_{ii} > 0$.

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
- $A \text{ simetrica} \rightarrow det(A) = \prod_i \lambda_i$
- $A \text{ cuadrada } \implies \sum_{i=0}^n \lambda_i = \text{traza(A)}$
- $A \text{ diagonal o triangular } \implies \text{ Los autovalores están en la diagonal}$
- $\lambda \text{ autovalor de } A  \text{ (invertible) } \implies \\
1/\lambda \text{ autovalor de } A^{-1} \text{ y } mg_A(\lambda) = mg_{A^{-1}}(1/\lambda) \text{ (solo porque son el mismo subespacio, la algebraica no se conserva)}$
- $\lambda \text{ autovalor de A } \implies \lambda^k \text{ autovalor de }A^k$
- $e \in \ker(A^T) \iff e \perp \operatorname{col}(A)$
- $D \text{ matriz diagonal } \implies ||D||_2 = max_{\substack{\lambda}} \{ |\lambda| : \lambda \text{ autovalor de D} \}=max{|D_{ii}|}$

### 8.1. Definida positiva
- $A \text{ es definida positiva } \iff \text{Para todo vector no nulo } \mathbf{x} \in \mathbb{R}^n \text{, se cumple que la forma cuadrática } \mathbf{x}^T A \mathbf{x} > 0$.
- $A \text{ es definida positiva } \iff \text{Todos los autovalores de } A \text{ son estrictamente positivos } (\lambda_i > 0)$.
- $A \text{ es definida positiva } \iff A=LU \text{ con } U_{ii} > 0 \text{ para } \forall i \in \{1 .. n\}$
- $A \text{ es definida positiva } \iff \text{Todos los menores principales de } A \text{ (determinantes de submatrices } k \times k \text{ superiores izquierdas) son estrictamente positivos } (\det(A_k) > 0 \text{ para } k=1, \dots, n)$.

### 8.3 Semidefinida positiva
- $A \text{ es semidefinida positiva } \iff \lambda_i \geq 0 \text{ (autovalores de A)}$

### 8.3. Diagonalizable
- $ A\in K^{N\times N} \text { es diagonalizable } \iff \text{ Existen N autovectores LI (vectores columna de A)}$
- $ A\in K^{N\times N} \text { es diagonalizable } \iff \text{ Para todo } \lambda_i \text{ autovalor de A, vale que } mg_A(\lambda_i) = ma_A(\lambda_i)$
- $ A\in K^{N\times N} \text { es diagonalizable } \iff A = PDP^{-1}, D \text{ diagonal }$
- $ A\in K^{N\times N} \text { es diagonalizable } \iff A^m = PD^mP^{-1}$
- $ A \text{ es diagonalizable si es semejante a una matriz diagonal}$
- $ A \text{ con autovalores distintos} \implies \text{ A diagonalizable}$

### 8.4. Multiplicidades
- $ mg_a(\lambda) = dim(E_\lambda)$
- $ ma_a(\lambda) = $ multiplicidad de lambda en el polinomio característico (cuántas veces aparece como autovalor).
- $1 \leq mg_a(\lambda) \leq ma_a(\lambda)$
- $A \in \mathbb{K}^{N \times N} \implies \sum_{i} ma(\lambda_i) = N$

### 8.5. Matrices semejantes
- $A, B \in K^{N\times N} \text{ son semejantes } \iff \exists c \in K^{N \times N} \text{ tal que } A = CBC^{-1}$ 
- $ \text{A y B son semejantes si } A= CBC^{-1}$
