---
title: "Resumen Prácticas 1 y 2 — ALC"
author: "Álgebra Lineal Computacional — UBA"
date: "2C 2025"
geometry: margin=2cm
fontsize: 11pt
header-includes:
  - \usepackage{amsmath}
  - \usepackage{amssymb}
  - \usepackage{booktabs}
  - \usepackage{xcolor}
  - \usepackage{tcolorbox}
  - \usepackage{mdframed}
  - \definecolor{azul}{RGB}{30,90,170}
  - \definecolor{gris}{RGB}{240,240,240}
  - \usepackage{titlesec}
  - \titleformat{\section}{\large\bfseries\color{azul}}{}{0em}{}[\titlerule]
  - \titleformat{\subsection}{\normalsize\bfseries}{}{0em}{}
---

# I. Subespacios Vectoriales

**Definición.** $U \subseteq V$ es subespacio $\Longleftrightarrow$ se cumplen las 3 condiciones:

1. $\mathbf{0} \in U$
2. $\mathbf{u}, \mathbf{v} \in U \Rightarrow \mathbf{u} + \mathbf{v} \in U$
3. $\mathbf{u} \in U,\ \alpha \in \mathbb{K} \Rightarrow \alpha\mathbf{u} \in U$

**Ejemplo típico:** soluciones de $A\mathbf{x} = \mathbf{0}$ siempre forman un subespacio.  
**Contraejemplo típico:** soluciones de $A\mathbf{x} = \mathbf{b}$ con $\mathbf{b} \neq \mathbf{0}$ NO son subespacio (el $\mathbf{0}$ no está).

**Sistema de generadores:**
$$\langle \mathbf{v}_1, \ldots, \mathbf{v}_m \rangle = \{ a_1\mathbf{v}_1 + \cdots + a_m\mathbf{v}_m : a_i \in \mathbb{K} \}$$

## Suma e intersección (Prop. 1.7)

| Dado | $S + T$ | $S \cap T$ |
|:---|:---|:---|
| $S$, $T$ por generadores | unir todos los generadores | calcular ecuaciones de cada uno y unirlas |
| $S$, $T$ por ecuaciones | calcular generadores de cada uno y unirlos | unir las ecuaciones |

**Fórmula de dimensión:**
$$\dim(U + V) = \dim(U) + \dim(V) - \dim(U \cap V)$$

**Suma directa:** $S \oplus T$ cuando $S \cap T = \{\mathbf{0}\}$, equivale a $\dim(S+T) = \dim(S) + \dim(T)$.

**Verificar $S \subseteq T$:** si $S$ está dado por generadores y $T$ por ecuaciones, verificar que cada generador de $S$ satisface las ecuaciones de $T$.

---

# II. Base e Independencia Lineal

**Independencia lineal:** $\{\mathbf{v}_1, \ldots, \mathbf{v}_m\}$ es L.I. $\Longleftrightarrow$ la única solución de $\sum a_i\mathbf{v}_i = \mathbf{0}$ es $a_1 = \cdots = a_m = 0$.

**Práctica:** armar la matriz con los vectores como **filas** y triangular. L.I. $\Longleftrightarrow$ todas las filas quedan no nulas.

**Base:** conjunto L.I. que genera el espacio.

**Prop. 1.1:** si $\{\mathbf{v}_1,\ldots,\mathbf{v}_m\}$ es L.I., cada vector del espacio generado se escribe de forma *única* como C.L.

**Prop. 1.3:** todas las bases de $V$ tienen la misma cantidad de elementos $\to$ **dimensión** de $V$.

**Para $n$ vectores en $V$ con $\dim(V) = n$:**
$$\text{L.I.} \Longleftrightarrow \text{base} \Longleftrightarrow \text{generan } V$$

- Menos de $n$ vectores: no pueden generar $V$
- Más de $n$ vectores: no pueden ser L.I.

## Algoritmos

**Base de subespacio dado por generadores:**

1. Poner los vectores como filas de $A$
2. Triangular (eliminación gaussiana)
3. Las filas no nulas son la base

**Base de subespacio dado por ecuaciones $A\mathbf{x} = \mathbf{0}$:**

1. Triangular $A$
2. Identificar variables libres (una por cada columna sin pivote)
3. Despejar variables dependientes en función de las libres
4. Un vector generador por cada variable libre

**Extensión de base (Prop. 1.6):** dada base $\{\mathbf{v}_1,\ldots,\mathbf{v}_s\}$ de $S \subseteq V$ con $\dim(V) = n$, existen $\mathbf{w}_1,\ldots,\mathbf{w}_{n-s}$ tales que $\{\mathbf{v}_1,\ldots,\mathbf{v}_s, \mathbf{w}_1,\ldots,\mathbf{w}_{n-s}\}$ es base de $V$.

## Rango

$$\text{rango}(A) = \dim(\text{im}(A)) = \text{máx número de filas (= columnas) L.I.}$$

$$\text{rango}(A) = \text{rango}(A^T) = \text{rango}(A^T A) = \text{rango}(A A^T)$$

$$\text{rango}(AB) \leq \min\{\text{rango}(A),\ \text{rango}(B)\}$$

Rango máximo $\Longleftrightarrow$ $A$ inversible $\Longleftrightarrow$ $A\mathbf{x} = \mathbf{b}$ tiene solución única para todo $\mathbf{b}$.

## Cambio de base

Sea $\mathcal{B} = \{\mathbf{b}_1,\ldots,\mathbf{b}_n\}$ base de $V$:

$$C_{\mathcal{B}\mathcal{E}} = \begin{pmatrix} | & & | \\ \mathbf{b}_1 & \cdots & \mathbf{b}_n \\ | & & | \end{pmatrix} \quad \text{(vectores de } \mathcal{B} \text{ como columnas, en coord. canónica)}$$

$$\mathbf{v}_\mathcal{E} = C_{\mathcal{B}\mathcal{E}} \cdot \mathbf{v}_\mathcal{B} \qquad \text{(base } \mathcal{B} \to \text{canónica)}$$

$$\mathbf{v}_\mathcal{B} = C_{\mathcal{E}\mathcal{B}} \cdot \mathbf{v}_\mathcal{E} = (C_{\mathcal{B}\mathcal{E}})^{-1} \cdot \mathbf{v}_\mathcal{E} \qquad \text{(canónica } \to \text{ base } \mathcal{B})$$

$$C(\mathcal{B},\mathcal{B}') = (C_{\mathcal{B}\mathcal{E}})^{-1} \cdot C_{\mathcal{B}'\mathcal{E}}$$

---

# III. Transformaciones Lineales

**Definición.** $f : V \to W$ es T.L. $\Longleftrightarrow$

1. $f(\mathbf{v} + \mathbf{v}') = f(\mathbf{v}) + f(\mathbf{v}')$
2. $f(a \cdot \mathbf{v}) = a \cdot f(\mathbf{v})$

**Consecuencias:** $f(\mathbf{0}) = \mathbf{0}$ siempre; $f(\alpha\mathbf{u} + \beta\mathbf{v}) = \alpha f(\mathbf{u}) + \beta f(\mathbf{v})$.

**Truco clave:** $f$ queda completamente determinada por sus valores en una base.  
Si $\mathcal{B} = \{\mathbf{v}_1,\ldots,\mathbf{v}_n\}$ es base de $V$ y $\mathbf{v} = a_1\mathbf{v}_1 + \cdots + a_n\mathbf{v}_n$, entonces:
$$f(\mathbf{v}) = a_1 f(\mathbf{v}_1) + \cdots + a_n f(\mathbf{v}_n)$$

**Representación matricial:** toda $f : \mathbb{K}^n \to \mathbb{K}^m$ tiene una única matriz $A \in \mathbb{K}^{m \times n}$ con $f(\mathbf{v}) = A\mathbf{v}$.  
Las **columnas** de $A$ son las imágenes de los vectores canónicos: columna $j = f(\mathbf{e}_j)$.

## Imagen y núcleo

$$\text{im}(A) = \{\mathbf{b} \in \mathbb{K}^m : \mathbf{b} = A\mathbf{x} \text{ para algún } \mathbf{x}\} \quad \to \text{ generado por las COLUMNAS de } A$$

$$\ker(A) = \{\mathbf{x} \in \mathbb{K}^n : A\mathbf{x} = \mathbf{0}\} \quad \to \text{ soluciones del sistema homogéneo}$$

**Teorema de la dimensión (Teorema 1.1):**
$$n = \dim(\ker(A)) + \dim(\text{im}(A))$$

## Clasificación

| Tipo | Definición | Condición equivalente |
|:---|:---|:---|
| **Monomorfismo** (inyectivo) | $\ker(f) = \{\mathbf{0}\}$ | $\dim(\ker) = 0$ |
| **Epimorfismo** (sobreyectivo) | $\text{im}(f) = W$ | $\dim(\text{im}) = \dim(W)$ |
| **Isomorfismo** | mono + epi | $\det(A) \neq 0$, $A$ inversible |

Si $f$ es isomorfismo: $f^{-1}$ tiene matriz $A^{-1}$.  
**Para hallar $A^{-1}$:** extender $(A \mid I)$ por Gauss hasta obtener $(I \mid A^{-1})$.

---

# IV. Normas y Número de Condición

## Normas vectoriales

$$\|\mathbf{x}\|_1 = |x_1| + \cdots + |x_n|$$
$$\|\mathbf{x}\|_2 = \sqrt{x_1^2 + \cdots + x_n^2} \qquad \text{(euclídea)}$$
$$\|\mathbf{x}\|_\infty = \max\{|x_1|, \ldots, |x_n|\}$$

**Axiomas de norma:**

1. $\|a\mathbf{v}\| = |a|\cdot\|\mathbf{v}\|$
2. $\|\mathbf{v}\| = 0 \Rightarrow \mathbf{v} = \mathbf{0}$
3. $\|\mathbf{u} + \mathbf{v}\| \leq \|\mathbf{u}\| + \|\mathbf{v}\|$ (desigualdad triangular)

**Desigualdad de Cauchy-Schwarz:**
$$\left|\sum \bar{u}_i v_i\right| \leq \|\mathbf{u}\|_2 \cdot \|\mathbf{v}\|_2$$

**Equivalencia de normas (Prop. 3.2):** en dim finita todas las normas son equivalentes. Para $\mathbf{x} \in \mathbb{K}^n$:
$$\|\mathbf{x}\|_\infty \leq \|\mathbf{x}\|_2 \leq \sqrt{n}\,\|\mathbf{x}\|_\infty$$
$$\tfrac{1}{\sqrt{n}}\|\mathbf{x}\|_1 \leq \|\mathbf{x}\|_2 \leq \|\mathbf{x}\|_1$$
$$\|\mathbf{x}\|_\infty \leq \|\mathbf{x}\|_1 \leq n\,\|\mathbf{x}\|_\infty$$

## Normas matriciales

**Norma subordinada (inducida):**
$$\|A\| = \max_{\mathbf{x} \neq \mathbf{0}} \frac{\|A\mathbf{x}\|}{\|\mathbf{x}\|} = \max_{\|\mathbf{x}\|=1} \|A\mathbf{x}\|$$

**Fórmulas cerradas:**
$$\|A\|_\infty = \max_{1 \leq i \leq n} \sum_{j=1}^n |a_{ij}| \qquad \text{(máxima suma de FILA en valor absoluto)}$$
$$\|A\|_1 = \max_{1 \leq j \leq n} \sum_{i=1}^n |a_{ij}| \qquad \text{(máxima suma de COLUMNA en valor absoluto)}$$

**Propiedades de normas subordinadas:**
$$\|A\mathbf{x}\| \leq \|A\|\cdot\|\mathbf{x}\| \qquad \|AB\| \leq \|A\|\cdot\|B\| \qquad \|A^k\| \leq \|A\|^k \qquad \|I\| = 1$$

**Equivalencia entre normas matriciales para $A \in \mathbb{R}^{n \times n}$:**
$$\tfrac{1}{\sqrt{n}}\|A\|_\infty \leq \|A\|_2 \leq \sqrt{n}\,\|A\|_\infty$$
$$\tfrac{1}{\sqrt{n}}\|A\|_1 \leq \|A\|_2 \leq \sqrt{n}\,\|A\|_1$$

**Norma de Frobenius** (no subordinada):
$$\|A\|_F = \sqrt{\sum_{i,j} |a_{ij}|^2} = \sqrt{\text{tr}(A^T A)} \qquad \|AB\|_F \leq \|A\|_F\cdot\|B\|_F$$

**Matrices ortogonales/unitarias** ($Q^T Q = I$):
$$\|Q\mathbf{v}\|_2 = \|\mathbf{v}\|_2 \qquad \|Q\|_2 = 1 \qquad \|QA\|_2 = \|A\|_2 \qquad \|QA\|_F = \|A\|_F$$

## Número de condición

**Definición (Def. 3.11):**
$$\kappa(A) = \|A\| \cdot \|A^{-1}\|$$

Si $A$ no es inversible: $\kappa(A) = +\infty$.

**Propiedades:**
$$\kappa(I) = 1 \qquad 1 \leq \kappa(A) \qquad \kappa(\alpha A) = \kappa(A) \quad (\alpha \neq 0)$$
$$\kappa(A) = \kappa(A^{-1}) \qquad \kappa(Q) = 1 \text{ para } Q \text{ ortogonal/unitaria}$$

**Fórmula central de error relativo:** dado $A\mathbf{x} = \mathbf{b}$ con error $\Delta\mathbf{b}$ en el dato:
$$\frac{1}{\kappa(A)} \cdot \frac{\|\Delta\mathbf{b}\|}{\|\mathbf{b}\|} \;\leq\; \frac{\|\Delta\mathbf{x}\|}{\|\mathbf{x}\|} \;\leq\; \kappa(A) \cdot \frac{\|\Delta\mathbf{b}\|}{\|\mathbf{b}\|}$$

Con errores tanto en $A$ como en $\mathbf{b}$:
$$\frac{\|\Delta\mathbf{x}\|}{\|\mathbf{x}\|} \;\lesssim\; \kappa(A) \left(\frac{\|\Delta\mathbf{b}\|}{\|\mathbf{b}\|} + \frac{\|\Delta A\|}{\|A\|}\right)$$

**Regla de oro:** si $\kappa(A) \approx 10^k$, podemos perder hasta $k$ dígitos significativos al resolver $A\mathbf{x} = \mathbf{b}$.

**Interpretación geométrica (Teorema 3.1):**
$$\frac{1}{\kappa(A)} = \inf_{B \text{ singular}} \frac{\|A - B\|}{\|A\|}$$

$\kappa(A)$ grande $\Longleftrightarrow$ $A$ está cerca (relativamente) de ser singular.  
**Ojo:** el determinante NO es buen indicador ($\kappa(\alpha I) = 1$ pero $\det(\alpha I) = \alpha^n \to 0$).

## Cómo calcular $\kappa_\infty(A)$

1. $\|A\|_\infty$ = máxima suma de fila de $A$ (en valor absoluto)
2. Calcular $A^{-1}$ vía Gauss sobre $[A \mid I]$
3. $\|A^{-1}\|_\infty$ = máxima suma de fila de $A^{-1}$ (en valor absoluto)
4. $\kappa_\infty(A) = \|A\|_\infty \cdot \|A^{-1}\|_\infty$
