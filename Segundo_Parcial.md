## ¿Cómo diagonalizo una matríz?

##### <u>Procedimiento:</u>
> 1. Determinar si A es diagonalizable (ej. Teo. Espectral).
> 2. Hallar autovalores $\lambda_1, ..., \lambda_n$ de A. 
> - Definir $D = \begin{pmatrix} \lambda_1 & 0 & ... & 0 & 0 \\
0 & \lambda_2 & ... & 0 & 0 \\
... & ... & ... & ... & ... \\
... & ... & ... & ... & \lambda_n \\
\end{pmatrix}$
> 3. Hallar autovectores $v_1, ..., v_n$ de A. 
> - Definir $Q = (v_1| ... | v_n)$.
> 4. Escribir $A = QDQ^*$


## ¿Cómo calculo una descomposicion Shur?

Teo: Toda matriz $A$ es unitariamente semejante a una matríz triangular superior ($\exists U \text{ unitaria y } T \text{ triangular: } A=UTU^*$).

Cualquier matríz cuadrada $A$ puede ser escrita como $A = QUQ^*$, con $Q$ matríz unitaria ($Q* = Q^{-1}$) y $D$ diagonal.

##### <u>Procedimiento:</u>
> 1. Encontrar un autovector
> 2. Completar una b.o.n para $Q_1$ con el autovector
> 3. Calcular $Q_1^* A Q_1$
> 4. Repetir los pasos anteriores para la submatriz del resultado, pisando los valores en las submatrices de $Q$ y $A$.


## ¿Cómo calculo una descomposicion SVD?

La descomposición en valores singulares de una matríz $A \in \mathbb{C} ^{m\times n}$ es un producto de la forma
$$A = U 
\Sigma V^*$$

Con:
- $U \in \mathbb{C}^{m \times m}$: Las columnas $u_1, ..., u_m$ vienen dadas por la relación $Av_j = \sigma_j u_j$ con $j=1...n$.
- $V \in \mathbb{C}^{n \times n}$: Las columnas son los autovectores (de $A^*A$). 
- $\Sigma \in \mathbb{C}^{m \times n}$: Diagonal real y no negativa.

Pueden pasar dos casos: 
- $m>n$: En tal situacion se completan las columnas de $U$ para tener una b.o.n en $\mathbb{C}^m$.
- $m<n$: Hay varios $v_j$ asociados a un autovalor 0. Si $\sigma_j=0$ para algun $j\leq \min(n,m)$ entonces se puede elegir $u_j$ completando la ortonormalidad de las columnas de $U$.

##### <u>Procedimiento:</u>
Sea $A\in\mathbb{C}^{m\times n}$...
> 1. Calcular $A^*A$
> 2. Calcular los autovalores y autovectores de $A^*A$
> 3. Formar las matrices $U, \Sigma, V^T$:
> - $V$ son los autovectores normalizados de $A^*A$.
> - $\Sigma$ es la matriz diagonal $\mathbb{C}^{m\times n}$ con los $\sigma_i = \sqrt{\lambda_i}$.
> - $U$ la calculamos con la relación $Av_i = \sigma_i u_i$ para cada $v_i$ $(i = 1...n)$ y luego completando una base ortonormal de para $u_n ... u_m$. 

Calcular SVD para $A^*$ sale de...
$$A^* = (U \Sigma V^*)^* = V \Sigma^t U^*$$

#### Propiedades para $A=U\Sigma V^*$
- $Im(A) = < u_1, ..., u_r >$
- $Im(A^t) = < v_1, ..., v_r >$
- $Nu(A) = < v_{r+1}, ..., v_{n} >$
- $Nu(A^t) = < u_{r+1}, ..., u_{m} > = Nu(A^tA)$

## Procesos de Markov.

#### Proceso de Markov
$v^{(k+1)} = Av^{(k)}$

#### Matríz de transición
Una matriz de transición $A$ cumple las siguientes propiedades:
- $A_{ij} \geq 0$
- $\sum_i A_ij = 1$, $\forall j$
- 1 es autovalor de $A$.
- Si $\lambda$ es autovalor entonces $|\lambda| < 1$.
- Si $\lambda \neq 1$ es autovalor y $x$ es autovector asociado entonces $\sum_i x_i = 0$
- $mg(\lambda) = ma(\lambda)$

Obs: No toda matríz de Markov es diagonalizable, pero si lo es $A = PDP^{-1} \rightarrow v^{(k+1)} = A^k v^{(0)} = PD^{k}P^{-1}v^{(0)}$

#### Estado de equilibrio
Un vector $v$ se dice estado de equilibrio si $Av=v$.
Es un autovector asociado para $\lambda=1$.

#### Estado límite
Dado $v^{(0)}$, es el estado $v^* = \lim_{k \rightarrow \inf} A^k v^0$ 
$\exists A^{\inf} \iff \text{no hay autovalores de modulo 1 que no sean el 1}$ 

Si $P$ tiene otros autovalores de módulo 1 además de $\lambda=1$, entonces no converge una única distribución estacionaria global, sino que puede haber oscilaciones o ciclos en la cadena de Markov.

> Si v^(0) es ortogonal a los autovectores asociados a los autovalores complejos de módulo 1, entonces v^(k) sí puede converger. 
> Si v^(0) tiene componente en dirección de esos autovectores, entonces v^(k) no converge (por ejemplo, oscila).

Conclusión clave
La existencia del estado límite v^(∞) puede depender de v^(0) cuando la matriz P tiene autovalores de módulo 1 distintos de λ=1. En cambio, si P^∞ existe (es decir, solo λ=1 está en la circunferencia unitaria), entonces el estado límite existe para cualquier v^(0) y es igual al vector estacionario π. 

## Cuadrados mínimos
Queremos aproximar una solución de $Ax=b$.

$$||Ax - b|| \rightarrow 0$$

Para resolverlo usamos ecuaciones normales
$$A^tAx= A^tb$$
$$x=(A^tA)^{-1}A^tb$$

Props:
- Si existen soluciones para $Ax=b \rightarrow z=A^tb$ es solucion. Si hay infinitas soluciones, $z=A^tb$ es la solucion de norma 2 mínima.
- $Ax=b$ tiene solucion $\iff AA^tb = b$
- Las ecuaciones normales no están bien condicionadas por eso las soluciones usan QR o SVD ($\text{cond}(A^tA) = \text{cond}_2(A)²$ y $\text{cond}_2(A) = \frac{\sigma_1}{\sigma_n}$).
- $A \text{ tiene columnas LI} \iff \text{ Cuadrados mínimos tiene única solución}$.

## Métodos Iterativos. Convergencia.

$$A = D + L + U$$

Definimos un método iterativo de la siguiente forma...
$$X = TX + c$$

Teo: $x^{n} \text{ converge } \iff \rho(T) < 1$

Con: $\rho(T) = \max{ |\lambda| : \lambda \text{ autovalor de T} }$

A veces se escribe $B=T$, como la matríz de iteración.

#### Método Jacobi
$x^{(n+1)} = -D^{-1}(L+U)x^{(n)} + D^{-1}b$

#### Método Gauss-Seidel
$x^{(n+1)} = -(D+L)^{-1}Ux^{(n)} + (D+L)^{-1}b$

#### Propiedades
- $T = -M^{⁻1}N \land \lambda \text{ autovalor de T} \iff det(\lambda M + N) = 0$ 
- El método converge en $n$ pasos si $T^n = 0$.
- $A$ matriz cuadrada y tridiagonal ($|a_{ij}=0| si |j-i|>1$) con $a_{ii}\neq 0$ para todo $i=1...n$. Entonces $\rho(B_{GS}) = \rho(B_j)^2$
- Para toda norma subordinada $||.||$ vale que:
    $$\rho(B) = lim_{n \rightarrow \inf} || B^n ||^{1/n}$$
- $x^* = Bx^* + c$
- $err_k = x_k - x^* = B_{gs} * err_{k-1} = B_{gs} * (x_{k-1} - x^*)$

#### Radio espectral
El radio espectral $\rho(T)$ de la matriz de iteración $T$ determina si un método iterativo converge y qué tan rápido lo hace. 
- Si $\rho(B)<1$: el método converge .
- Si $\rho(B)=1$: puede no converger o converger muy lentamente.
- Si $\rho(B)>1$ : el método diverge .
Cuanto menor sea $\rho(B)$, más rápida será la convergencia .

### Propiedades

- $A$ es hermitana $\iff$ $A=A^*$. 

- Las matrices hermitianas tienen autovalores reales y son diagonalizables ortogonalmente con $A=UDU*$

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

- Pseudoinversa: $$A^{+} = (A^tA)^{-1}A^t$$


