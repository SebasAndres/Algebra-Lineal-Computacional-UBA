## I. Descomposiciones y Propiedades de Matrices

### Matrices Especiales
| Tipo de Matriz | Propiedad Definitoria | Propiedades Clave |
| :--- | :--- | :--- |
| **Hermitiana** ($\mathbb{C}^{n \times n}$) | $\boldsymbol{A}^* = \boldsymbol{A}$ (Simétrica en $\mathbb{R}$) | Autovalores **reales**. Diagonalizable unitariamente. |
| **Unitaria/Ortogonal** ($\mathbb{C}^{n \times n}$/$\mathbb{R}^{n \times n}$) | $\boldsymbol{U}^{-1} = \boldsymbol{U}^*$ (o $\boldsymbol{Q}^{-1} = \boldsymbol{Q}^t$) | Preserva la norma 2: $||\boldsymbol{U}\boldsymbol{x}||_2 = ||\boldsymbol{x}||_2$. |
| **Idempotente** ($\boldsymbol{A}^2 = \boldsymbol{A}$) | $\boldsymbol{A}^2 = \boldsymbol{A}$ | Autovalores son exclusivamente $\boldsymbol{\lambda \in \{0, 1\}}$. |
| **Nilpotente** | $\boldsymbol{B}^k = \boldsymbol{0}$ para algún $k$. | Todos los autovalores son cero. El método iterativo $\boldsymbol{x}^{(k+1)} = \boldsymbol{B}\boldsymbol{x}^{(k)}+\boldsymbol{c}$ converge exactamente en $k$ pasos. |

### Descomposición QR
* **Factorización:** $\boldsymbol{A} = \boldsymbol{Q}\boldsymbol{R}$, donde $\boldsymbol{Q}$ es unitaria/ortogonal y $\boldsymbol{R}$ es triangular superior.
* **Propiedad Clave:** La matriz $\boldsymbol{A}$ y la matriz **$\boldsymbol{R}$ tienen los mismos valores singulares**. Esto se debe a que $||\boldsymbol{A}\boldsymbol{x}||_2 = ||\boldsymbol{Q}\boldsymbol{R}\boldsymbol{x}||_2 = ||\boldsymbol{R}\boldsymbol{x}||_2$.

---

## II. Valores Singulares (SVD) y Normas

### SVD (Descomposición en Valores Singulares)
* **Fórmula:** $\boldsymbol{A} = \boldsymbol{U} \boldsymbol{\Sigma} \boldsymbol{V}^*$
    * $\boldsymbol{\Sigma}$: Matriz diagonal que contiene los **valores singulares** $\boldsymbol{\sigma_i \ge 0}$ ordenados de forma decreciente ($\sigma_1 \ge \sigma_2 \ge \ldots$).
    * $\boldsymbol{U}, \boldsymbol{V}$: Matrices unitarias.
* **Cálculo:** Los valores singulares $\sigma_i$ son las raíces cuadradas de los autovalores de $\boldsymbol{A}^* \boldsymbol{A}$ (o $\boldsymbol{A}\boldsymbol{A}^*$).
* **Rango:** El **rango** de $\boldsymbol{A}$ es igual al número de **valores singulares no nulos** ($\sigma_i \ne 0$).

### Normas Matriciales (Norma 2)
* **Norma 2 (Espectral):** La norma 2 de una matriz es igual a su mayor valor singular.
    $$\boldsymbol{||\boldsymbol{A}||_2 = \sigma_{\max} = \sigma_1}$$
* **Número de Condición (Norma 2):** Mide la sensibilidad de la solución de un sistema $\boldsymbol{A}\boldsymbol{x} = \boldsymbol{b}$ a perturbaciones.
    $$\text{cond}_2(\boldsymbol{A}) = \frac{\sigma_{\max}}{\sigma_{\min}} = \frac{\sigma_1}{\sigma_n}$$
* **Acotación de Autovalores:** El radio espectral de $\boldsymbol{A}$ está acotado por cualquier norma de $\boldsymbol{A}$ (incluyendo la norma 2).
    $$\boldsymbol{\rho(\boldsymbol{A}) \le ||\boldsymbol{A}||}$$

---

## III. Métodos Iterativos (Jacobi y Gauss-Seidel)

### Condición de Convergencia
Ambos métodos, $\boldsymbol{x}^{(k+1)} = \boldsymbol{B} \boldsymbol{x}^{(k)} + \boldsymbol{c}$, convergen si y solo si el **radio espectral** de la matriz de iteración $\boldsymbol{B}$ es menor que 1:
$$\boldsymbol{\rho(\boldsymbol{B}) = \max \{|\lambda|: \lambda \text{ es autovalor de } \boldsymbol{B}\} < 1}$$

### Criterios Suficientes
| Propiedad de $\boldsymbol{A}$ | Convergencia de Jacobi | Convergencia de Gauss-Seidel |
| :--- | :--- | :--- |
| **Estrictamente Diagonal Dominante** | Sí, siempre. | Sí, siempre. |
| **Simétrica y Definida Positiva (SDP)** | Converge (en $2\times 2$) | Sí, siempre. |

### Paralelización
* **Jacobi:** **Es altamente paralelizable** 🚀. La fórmula de iteración utiliza solo los valores de la iteración anterior $\boldsymbol{x}^{(k)}$ para calcular cada componente de $\boldsymbol{x}^{(k+1)}$. Todos los cálculos pueden ejecutarse simultáneamente.
* **Gauss-Seidel:** **No es fácilmente paralelizable**. Utiliza los valores ya actualizados de la iteración actual, lo que impone una dependencia secuencial ($x_i^{(k+1)}$ depende de $x_j^{(k+1)}$ para $j < i$).

---

## IV. Cuadrados Mínimos (Least Squares)

### Problema Fundamental
Encontrar el vector $\hat{\boldsymbol{x}}$ que **minimiza el error** $||\boldsymbol{A}\boldsymbol{x}-\boldsymbol{b}||_2$.

### Solución por Pseudo-Inversa
La solución de cuadrados mínimos de norma 2 mínima se obtiene mediante la **pseudo-inversa** de Moore-Penrose, $\boldsymbol{A}^{\dagger}$:
$$\boldsymbol{\hat{x} = \boldsymbol{A}^{\dagger}\boldsymbol{b}}$$

### Definición de la Pseudo-Inversa (Usando SVD)
Si $\boldsymbol{A} = \hat{\boldsymbol{U}}\hat{\boldsymbol{\Sigma}}\hat{\boldsymbol{V}}^t$ es la SVD reducida (para matrices con rango $r$), la pseudo-inversa es:
$$\boldsymbol{A}^{\dagger} = \hat{\boldsymbol{V}} \hat{\boldsymbol{\Sigma}}^{-1} \hat{\boldsymbol{U}}^t$$

### Compresión de Imágenes (Aplicación de SVD)
* **Aproximación de rango $r$ ($\boldsymbol{A}_r$):** La mejor aproximación a $\boldsymbol{A}$ de rango $r$ (en norma 2) se obtiene truncando la SVD, usando solo los primeros $r$ valores singulares.
* **Error de Aproximación:** El error de esta aproximación está dado por el siguiente valor singular:
    $$\frac{||\boldsymbol{A}-\boldsymbol{A}_r||_2}{||\boldsymbol{A}||_2} = \frac{\sigma_{r+1}}{\sigma_{1}}$$