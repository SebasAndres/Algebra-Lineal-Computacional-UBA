### **I. Propiedades de Matrices y Descomposiciones**

#### **1. Matrices Especiales y Sus Propiedades**
* **Matriz de Markov (Estocástica):**
    * **Definición:** Matriz cuadrada $P$ con entradas no negativas ($P_{ij} \ge 0$) tal que la suma de los elementos de cada columna es 1 ($\sum_i P_{ij} = 1$).
    * **Propiedades Clave:**
        * Siempre tiene al autovalor $\lambda=1$.
        * Si $P$ es diagonalizable y el único autovalor con $|\lambda|=1$ es $\lambda=1$, entonces existe una matriz límite $P^{\infty} = \lim_{k \to \infty} P^k$.
        * El estado límite $v^{(\infty)}$ es el autovector asociado a $\lambda=1$, normalizado tal que la suma de sus componentes sea 1 (vector de probabilidad).
        * La multiplicidad algebraica de $\lambda=1$ suele ser 1 en problemas prácticos irreducibles, implicando que el autoespacio tiene dimensión 1 ($mg_P(1)=1$).
    * **Cálculo de $P^{\infty}$:** Las columnas de $P^{\infty}$ son todas iguales al vector límite $v^{(\infty)}$.

* **Matriz Ortogonal ($Q \in \mathbb{R}^{n \times n}$):**
    * **Definición:** $Q^T = Q^{-1}$, o equivalentemente $Q^T Q = Q Q^T = I$.
    * **Propiedades:**
        * Preserva la norma 2: $||Qx||_2 = ||x||_2$.
        * Las columnas (y filas) forman una base ortonormal de $\mathbb{R}^n$.
        * Si en la SVD de una matriz $A = U \Sigma V^T$, se cumple que $\Sigma = I$, entonces $A$ es ortogonal (producto de ortogonales).

* **Matriz Hermitiana ($A^* = A$) / Simétrica ($A^T = A$):**
    * Todos sus autovalores son reales.
    * Es diagonalizable ortogonalmente (o unitariamente).

#### **2. Relaciones entre Matrices Reales**
* **Autovalores de $A^T A$ y $A A^T$:**
    * Si $A \in \mathbb{R}^{m \times n}$, las matrices $A^T A$ y $A A^T$ tienen los **mismos autovalores no nulos**. Estos autovalores son el cuadrado de los valores singulares de $A$ ($\sigma_i^2$).
* **Norma Matricial Inducida (Norma 2):**
    * $||A||_2 = \sigma_1$ (mayor valor singular).
    * Se cumple que $||Ax||_2 \ge \sigma_n ||x||_2$ (donde $\sigma_n$ es el menor valor singular). Si $||Ax||_2 \ge ||x||_2$ para todo $x$, entonces $\sigma_n \ge 1$.

---

### **II. Descomposiciones Matriciales**

#### **1. Descomposición en Valores Singulares (SVD)**
* **Forma:** $A = U \Sigma V^T$.
* **Rango Completo:** Si $A$ es cuadrada y de rango completo, todos sus valores singulares son no nulos ($\sigma_i > 0$).
* **Caracterización de Matrices Ortogonales:** Una matriz cuadrada real $A$ es ortogonal si y solo si todos sus valores singulares son iguales a 1 ($\sigma_i = 1$).
    * *Demo:* Si $\max_{||x||_2=1} ||AA^T x||_2 = 1 \implies \sigma_1^2=1 \implies \sigma_1=1$. Si además $||Ax||_2 \ge ||x||_2 \implies \sigma_n \ge 1$. Entonces todos los $\sigma_i = 1$.

#### **2. Descomposición de Schur**
* **Forma:** $A = U T U^*$.
* **Utilidad:** Permite calcular potencias de matrices, aunque en la práctica para matrices de Markov pequeñas suele ser más directo calcular autovalores y autovectores o multiplicar matrices.

---

### **III. Cuadrados Mínimos Lineales y No Lineales**

#### **1. Planteo General**
* Dado un conjunto de datos $(x_i, y_i, z_i)$, se busca la función $f(x,y)$ que minimice el error cuadrático medio: $\min \sum (z_i - f(x_i, y_i))^2$.
* El problema se reduce a resolver el sistema sobredeterminado $Ax = z$ en el sentido de cuadrados mínimos, buscando $\hat{x}$ tal que $A^T A \hat{x} = A^T z$ (Ecuaciones Normales).

#### **2. Modelos Lineales**
* Ejemplo: $z = ax + by + c$.
* Se plantea la matriz $A$ con las variables independientes (incluyendo una columna de 1s para el término independiente $c$).
* Se resuelve el sistema normal $A^T A \hat{x} = A^T z$.

#### **3. Modelos No Lineales (Linealización)**
* Ejemplo: $z = y + \frac{3}{ax+b}$.
* **Estrategia:** Transformar la ecuación para que sea lineal en los parámetros desconocidos ($a, b$).
    * Despeje: $ax+b = \frac{3}{z-y}$.
    * Cambio de variable: Definir una nueva variable respuesta $\Phi_i = \frac{3}{z_i - y_i}$.
    * Nuevo sistema lineal: $ax_i + b = \Phi_i$.
* **Comparación de Modelos:** Para decidir qué modelo es mejor, se debe calcular el **residuo original** $\sum (z_i - f(x_i, y_i))^2$ para los parámetros hallados en cada caso. **Cuidado:** No comparar el error del sistema linealizado directamente con el del sistema original, siempre volver a la variable original $z$.

---

### **IV. Métodos Iterativos para Sistemas Lineales**

#### **1. Método de Jacobi y Gauss-Seidel**
* **Jacobi ($B_J$):** $B_J = -D^{-1}(L+U)$.
    * Si $D=I$, entonces $B_J = -(L+U)$938].
* **Convergencia:** Converge $\iff \rho(B) < 1$.

#### **2. Métodos de Relajación y SOR**
* **Matriz de Iteración General ($B_{\omega}$):** Se define como una combinación lineal.
    * Ejemplo (Jacobi Relajado): $B_{\omega} = \omega B_J + (1-\omega)I$.
    * Ejemplo (SOR): $B(\omega) = (D+\omega L)^{-1}((1-\omega)D - \omega U)$
* **Relación con el Sistema Original:**
    * Para probar que el método resuelve $Ax=b$, se debe demostrar que la ecuación de punto fijo $x = B_{\omega}x + \omega b$ (u otro término independiente) es equivalente a $Ax=b$.
* **Análisis de Convergencia en Función de $\omega$:**
    * Calcular el polinomio característico de $B_{\omega}$ en función de $\omega$ y $\lambda$: $\det(B_{\omega} - \lambda I) = 0$.
    * Hallar las raíces $\lambda(\omega)$.
    * Imponer la condición $|\lambda(\omega)| < 1$ para encontrar el rango válido de $\omega$.
* **Velocidad de Convergencia:**
    * Un método converge más rápido que otro si su radio espectral es menor ($\rho(B_1) < \rho(B_2)$).
    * Puede suceder que para ciertos rangos de $\omega$ el método relajado sea mejor, o que nunca lo sea (dependiendo de la matriz $A$).
* **Propiedad del Determinante (SOR):** $\det(B(\omega)) = (1-\omega)^n$. Condición necesaria para convergencia: $\omega \in (0, 2)$.

---

### **V. Consejos Prácticos para la Resolución**
1.  **Polinomios Característicos de Grado 3:** Si aparece un polinomio de grado 3 en una matriz de Markov, recordar que $\lambda=1$ es siempre raíz. Usar Ruffini para factorizarlo.
2.  **Cálculo de Potencias de Matrices:** Para calcular $v^{(k)} = P^k v^{(0)}$, a veces es computacionalmente más barato hacer iteraciones vector-matriz sucesivas ($v^{(i+1)} = P v^{(i)}$) que calcular la potencia de la matriz $P^k$ primero.
3.  **Comparación de Errores:** Al linealizar modelos no lineales, prestar atención a denominadores pequeños que pueden amplificar el error en el modelo original.
4.  **Cotas de Valores Singulares:** Usar la propiedad $||Ax||_2 \ge \sigma_n ||x||_2$ para acotar inferiormente el menor valor singular.