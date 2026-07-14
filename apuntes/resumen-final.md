# Resumen Final — Álgebra Lineal Computacional

Unifica `resumen-1er-parcial.md` y `resumen-2do-parcial.md`, reorganizado según los 5 bloques temáticos del final (28/07). Pensado para repasar rápido cuando te perdés en un concepto o algoritmo — no para leer de punta a punta.

⚠️ = temas que marcaste como más difíciles (Householder, Schur). Tienen la sección más desarrollada a propósito.

## Índice
0. [Tabla rápida: existencia y costo de descomposiciones](#0-tabla-rápida-existencia-y-costo-de-descomposiciones)
1. [Subespacios, Bases, TL, Normas, Número de Condición](#1-subespacios-bases-tl-normas-número-de-condición)
2. [LU, Cholesky, Ortogonalidad, Proyectores, QR (⚠️ Householder)](#2-lu-cholesky-ortogonalidad-proyectores-qr)
3. [Autovalores, Diagonalización, Markov](#3-autovalores-diagonalización-markov)
4. [⚠️ Schur, SVD, Pseudoinversa](#4-schur-svd-pseudoinversa)
5. [Mínimos Cuadrados, Métodos Iterativos](#5-mínimos-cuadrados-métodos-iterativos)
6. [Apéndice: errores comunes / V-F rápidas](#6-apéndice-errores-comunes--v-f-rápidas)

---

## 0. Tabla rápida: existencia y costo de descomposiciones

| Descomposición | Forma | Requisito | ¿Existe siempre? | Costo | Notas |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **LU** | $A=LU$ | Cuadrada, menores principales $\det(A_k)\neq0$ para $k=1..n-1$ | NO ❌ (sí con permutación $PA=LU$) | $O(n^3/3)$ | Único si $A$ inversible |
| **Cholesky** | $A=\hat L\hat L^T$ | Cuadrada, simétrica y **definida positiva** | NO ❌ | $O(n^3/3)$ | La más restrictiva; única si existe |
| **QR** | $A=QR$ | $m\geq n$, cualquier $A$ | SÍ ✅ | Householder: $O(n^3)$ real (naive armando $H_k$: $O(n^4)$) · GS: $O(n^3)$ | Única si columnas LI o $R_{ii}>0$ |
| **Schur** | $A=UTU^*$ | Cuadrada | SÍ ✅ | — | $T$ triangular superior, $U$ unitaria |
| **SVD** | $A=U\Sigma V^*$ | Ninguno — **cualquier matriz** $m\times n$ | SÍ ✅ | — | La más general de todas |

**Regla de costo general:** multiplicar $(p\times q)\cdot(q\times r)$ cuesta $O(p\cdot q\cdot r)$.

---

## 1. Subespacios, Bases, TL, Normas, Número de Condición

### 1.1 Dimensión y subespacios
- **Producto matriz-vector = combinación lineal de columnas:** si $B=[b_1|b_2|\dots|b_n]$ (columnas $b_i$), entonces
  $$B\begin{pmatrix}a_1\\a_2\\\vdots\\a_n\end{pmatrix} = a_1b_1+a_2b_2+\dots+a_nb_n$$
  Demostración: $(a_1,\dots,a_n)^T=\sum a_ie_i$, y $Be_i=b_i$ (columna $i$), entonces por linealidad $B\sum a_ie_i=\sum a_i(Be_i)=\sum a_ib_i$. Es **la** forma correcta de pensar $Bx$ (en vez de "fila por columna"). Consecuencia directa: $Bx=0$ tiene solución no trivial $\iff$ las columnas de $B$ son LD (ver 1.2).
- **Notación:** $f:\mathbb V\to\mathbb W$ es una TL entre $\mathbb K$-espacios vectoriales; $\mathbb V$ = **dominio** (donde vive el vector de entrada), $\mathbb W$ = **codominio** (espacio de llegada declarado). $\operatorname{Im}(f)\subseteq\mathbb W$ es lo que $f$ efectivamente alcanza (ver "Ojo" más abajo).
- **Teorema de la dimensión (rango-nulidad):** $\dim(\operatorname{Im}(f)) + \dim(\ker(f)) = \dim(\mathbb V)$.
  - Vale **siempre**, para toda TL (solo requiere $\dim(\mathbb V)$ finita) — no presupone mono/epi/iso. Al revés: mono/epi/iso son casos particulares que se derivan de él. Ej: $f:\mathbb R^3\to\mathbb R^3$, $f(x,y,z)=(x,y,0)$ — ni mono ni epi, pero $\dim(\ker f)=1$, $\dim(\operatorname{Im} f)=2$, y $1+2=3$ igual se cumple.
- **Mono/epi/isomorfismo:**
  - $f$ inyectiva (**monomorfismo**) $\iff \ker(f)=\{0\}$.
  - $f$ sobreyectiva (**epimorfismo**) $\iff \operatorname{Im}(f)=\mathbb W$.
  - $f$ biyectiva (**isomorfismo**) $\iff$ mono y epi a la vez $\Rightarrow \dim(\mathbb V)=\dim(\mathbb W)$ (espacios isomorfos).
  - Si $\dim(\mathbb V)=\dim(\mathbb W)$ (dim. finita): alcanza con probar una sola (mono $\iff$ epi $\iff$ iso).
  - **Ojo (codominio ≠ imagen):** el codominio $\mathbb W$ es el espacio de llegada **declarado** (parte de la definición de $f:\mathbb V\to\mathbb W$); la imagen es el espacio de llegada **real** (lo que $f$ efectivamente alcanza). Siempre $\operatorname{Im}(f)\subseteq\mathbb W$; epimorfismo es cuando esa inclusión es **igualdad**, no una propiedad trivial de $\operatorname{Im}(f)$ en sí misma (todo elemento de $\operatorname{Im}(f)$ tiene preimagen por definición, eso no dice nada). Ej: $f:\mathbb R^2\to\mathbb R^3$, $f(x,y)=(x,y,0)$ — codominio $\mathbb R^3$, imagen el plano $\{(x,y,0)\}\subsetneq\mathbb R^3$: no es epimorfismo.
- $\dim(A+B) = \dim(A)+\dim(B)-\dim(A\cap B)$.
- $rk(A) = \dim(\operatorname{Im}(A))$, $\ker(A) = \{x : Ax=0\}$.
- **Suma:** $S+T=\langle$ generadores de $S$ y $T$ juntos $\rangle$.
- **Intersección:** resolver $S\alpha - T\beta = 0$ (sistema homogéneo con $[S\,|\!-T]$); los vectores $S\alpha$ que resultan forman una base de $S\cap T$.
- **Suma directa:** $L\oplus M \iff L\cap M=\{0\} \land L+M=\mathbb K^n$.
- $S^\perp = \{v : \langle v,s\rangle=0\ \forall s\in S\}$.
- Para extraer una base de generadores: ponerlos como filas, triangular, descartar filas que dan 0 o CL de otras.

### 1.2 Inversibilidad
$A$ **singular** $\iff \lambda=0$ autovalor $\iff \det(A)=0 \iff \dim(\ker A)>0 \iff$ columnas/filas LD.
$A$ **inversible** $\iff \det(A)\neq0 \iff$ columnas/filas LI $\iff$ columnas forman base de $\mathbb K^n$.
- Diagonal dominante estricta $\implies$ inversible.

### 1.3 Transformaciones lineales y cambio de base
- Extraer una base dado un conjunto de generadores. **Opcion 1:** Poner los vectores como filas y triangular la matriz, luego la base viene del sistema resultante. **Opcion 2:** Poner los vectores como columnas y triangular la matríz, luego las columnas que quedan con pivotes (diagonal o valores != 0) en las filas corresponden a columnas cuyos vectores originales forman una base.
- Extender una base para completar el subespacio. Se puede probar agregando los vectores como fila y al final un canónico (ir probando), triangular el sistema y ver si la fila se hace cero, si no se hace cero entonces los generadores anteriores + el e_i forman una base.
  - **Ojo:** no alcanza con "mirar" la fila del $e_i$ agregado — hay que completar la triangulación (reducirla contra **todos** los pivotes ya existentes que caigan en su misma columna). Si $e_i$ es redundante, esa reducción la va a llevar a fila nula sin importar que esté al final; si hace falta, puede requerir reordenar filas (pivoteo) para que el pivote quede visible. No es que el último $e_i$ "sobreviva" siempre por estar en el fondo del sistema.
- $C_{BE} @ [(x,y,z)]_B = [(x,y,z)]_E \iff [(x,y,z)]_B = C_{EB} @ [(x,y,z)]_E $
- **Requisitos para que $f$ sea TL** (equivalentes, cualquiera de las dos formas alcanza):
  - Forma separada: $f(u+v)=f(u)+f(v)$ (aditividad) **y** $f(\alpha u)=\alpha f(u)$ (homogeneidad), $\forall u,v\in\mathbb V, \forall\alpha\in\mathbb K$.
  - Forma combinada (la más rápida de chequear en un ejercicio): $f(\alpha u+v)=\alpha f(u)+f(v)$.
  - Consecuencias (no son lo que hay que probar, salen gratis si $f$ es TL): $f(0)=0$, $f(-v)=-f(v)$.
- $[f]_{BB'}$: columnas = coordenadas de $f(v_i)$ en base $B'$. $[f(v)]_{B'} = [f]_{BB'}[v]_B$.
- $C_{BB'}$: columnas = vectores de $B$ escritos en $B'$. $C_{B'B}=(C_{BB'})^{-1}$. $C_{BB''}=C_{B'B''}C_{BB'}$.
- Práctico: $C_{BE}$ = poner vectores de $B$ como columnas. $C_{EB}=(C_{BE})^{-1}$; si $B$ ortonormal, $C_{EB}=C_{BE}^T$.
- $[f]_{BE} = [f]_{EE}\,C_{BE}$.
- **Clave:** dos matrices representan la misma TL en distintas bases $\iff$ son **semejantes** $\iff \exists C: A=CBC^{-1}$ $\implies$ comparten det, traza, autovalores.

### 1.4 Normas vectoriales y matriciales
- Norma: $\|av\|=|a|\|v\|$; $\|v\|=0\Rightarrow v=0$; desigualdad triangular.
- $\|v\|_1=\sum|v_i|$, $\|v\|_\infty=\max|v_i|$, $\|v\|_p=(\sum|v_i|^p)^{1/p}$.
- Cauchy-Schwarz: $|x^*y|\leq\|x\|\|y\|$.
- Norma matricial inducida: $\|A\| = \max_{\|x\|=1}\|Ax\|$.
- $\|A\|_\infty$ = máx suma de **filas** (en módulo). $\|A\|_1$ = máx suma de **columnas**.
- $\|A\|_2 = \sigma_1(A)$ (mayor valor singular) — ver Tema 4.
- $\|Ax\|\leq\|A\|\|x\|$, $\|AB\|\leq\|A\|\|B\|$.
- Relaciones: $\tfrac{1}{\sqrt n}\|A\|_\infty<\|A\|_2<\sqrt n\|A\|_\infty$; $\tfrac1{\sqrt n}\|A\|_2<\|A\|_1<\sqrt n\|A\|_2$; $\|A\|_\infty<\|A\|_1<n\|A\|_\infty$.

### 1.5 Número de condición
- $\kappa_*(A) = \|A\|_*\|A^{-1}\|_*$, rango $[1,\infty)$.
- $\kappa$ grande $\implies$ mal condicionado, sensible a errores de redondeo.
- $\dfrac1{\kappa(A)} \leq \inf_{H \text{ singular}}\dfrac{\|A-H\|}{\|A\|}$ — el condicionamiento acota la distancia relativa a una matriz singular.
- Error relativo: $\dfrac{\|x-\tilde x\|}{\|x\|} \leq \kappa(A)\dfrac{\|b-\tilde b\|}{\|b\|}$.
- **Ojo (Práctica 2, Ej 24):** $\det(A)\to0$ **no implica** mal condicionamiento (ej: matriz escalar $D_n=\varepsilon I$, $\det\to0$ pero $\kappa(D_n)=1$).
- **cond mal condicionado ≠ convergencia lenta de iterativos** — son propiedades independientes (ver Tema 5).

---

## 2. LU, Cholesky, Ortogonalidad, Proyectores, QR

### 2.1 LU
$A=LU$, $L$ triangular inferior (diagonal 1), $U$ triangular superior.
- **Existencia (sin pivoteo):** todos los pivotes $\neq0$ en Gauss $\iff \det(A_k)\neq0$ para $k=1,\dots,n-1$ (menores principales top-left). No hace falta para $k=n$.
- Un cero en $a_{11}$ ya mata la existencia. Un cero en otra posición de la diagonal no implica que algún menor sea 0.
- **Única** si $A$ inversible.
- **Algoritmo:** triangular con operaciones elementales → $U$; en $L_{ij}$ guardar el inverso (con signo) del multiplicador usado para anular la celda $ij$.
- **$PA=LU$**: existe siempre para cuadrada, aplicando permutación de filas primero.

### 2.2 Cholesky
$A=\hat L\hat L^T$. Existe $\iff A$ **simétrica y definida positiva**. Única si existe.
**Algoritmo:**
1. $A=LU$.
2. $D_{ii}=U_{ii}$ (diagonal), resto 0.
3. $D_1$: raíz de cada entrada de $D$.
4. $\hat L = L D_1$.
5. $A=\hat L\hat L^T$.

### 2.3 Definida positiva — condiciones equivalentes
$A$ DP $\iff$
- $x^TAx>0\ \forall x\neq0$
- todos los $\lambda_i>0$
- $A=LU$ con $U_{ii}>0\ \forall i$
- todos los menores principales $\det(A_k)>0$

**Ojo:** entradas positivas ≠ definida positiva. Contraejemplo: $\begin{bmatrix}1&2\\2&1\end{bmatrix}$, $\det=-3<0$.
Cadena: $\exists\det(A_k)\leq0 \iff \exists\lambda_i\leq0 \iff$ no DP $\iff$ no tiene Cholesky.
Semidefinida positiva (SDP): $\lambda_i\geq0$.

### 2.4 Proyectores (generales)
- $P$ proyector $\iff P^2=P$. Autovalores solo pueden ser $0$ o $1$.
- $\operatorname{Im}(P)=\operatorname{Nu}(I-P)$, $\operatorname{Im}(I-P)=\operatorname{Nu}(P)$ (en la imagen, $Px=x$).
- $V = \operatorname{Im}(P) \oplus \operatorname{Nu}(P)$ siempre.
- **Armado:** si $M=[v_1|\dots|v_n]$ y quiero $P(v_i)$ dado (columnas de $R$): $PM=R \Rightarrow P=RM^{-1}$.

### 2.5 Proyección ortogonal
- $P_S$ es ortogonal $\iff P_S^2=P_S \land P_S^T=P_S$ (simetría extra). $\operatorname{Nu}(P)\perp\operatorname{Im}(P)$.
- $x = P_S(x) + P_{S^\perp}(x)$, descomposición única; $x-P_S(x)\in S^\perp$ (vector error).
- $A$ con columnas LI (no ortonormales): $P_S = A(A^TA)^{-1}A^T$.
- $Q$ con columnas ortonormales ($Q^TQ=I$): $P_S = QQ^T$.
- Desde BON $\{v_1,\dots,v_r\}$: $P_S=\sum v_iv_i^T$; $P_S(x)=\sum(v_i^Tx)v_i$.
- Complementario: $P_{S^\perp}=I-P_S$. $\operatorname{Im}(P_S)=S$, $\operatorname{Nu}(P_S)=S^\perp$.
- Matrices $u_iu_i^T$ siempre tienen rango 1. $A=uv^T \Rightarrow \operatorname{Im}(A)=\langle u\rangle$, $\operatorname{Nu}(A)=\langle v\rangle^\perp$.

### 2.6 Gram-Schmidt
$a_1,\dots,a_n$ columnas de $A$. $v_1=a_1$; para $i\geq1$: $v_{i+1}=a_{i+1}-\sum_{j\leq i} p_{v_j}(a_{i+1})$, con $p_a(b)=\frac{a^Tb}{a^Ta}a$. Normalizar: $v_i'=v_i/\|v_i\|$.
- $Q$ ortonormal $\Rightarrow Q^T=Q^{-1}$, $\|Qx\|=\|x\|$ (preserva norma).

### 2.7 QR
Existe para toda $A\in\mathbb R^{m\times n}$, $m\geq n$ (con rango completo de columnas; si no, algunas filas de $R$ son nulas). Única si columnas LI o $R_{ii}>0$.
- $\operatorname{Im}(A)=\operatorname{Im}(Q)$ pero **$Ax\neq Qx$ en general** (son TL distintas).
- **Vía Gram-Schmidt:** $Q=[q_1|\dots|q_n]$ con $q_i=v_i/\|v_i\|$. $R=Q^TA$, triangular superior con $r_{kk}=\|v_k\|$.

### 2.8 ⚠️ QR vía Householder

**Idea:** en cada paso, reflejar la columna actual sobre el eje $e_1$ del subespacio restante, usando una matriz de reflexión (Householder), hasta triangular $A$ por completo.

**Algoritmo:**
1. $v_k$ = columna $k$ de la submatriz restante $A[k:n,\,k:n]$. $w_k = (\|v_k\|,0,\dots,0)^T$.
2. $u_k = \dfrac{v_k-w_k}{\|v_k-w_k\|}$.
3. $H_k = I - 2u_ku_k^T$ (reflector, actúa solo sobre el bloque restante).
4. $A \leftarrow H_kA$ (anula todo debajo de la diagonal en la columna $k$).
5. Repetir hasta $A$ triangular superior $= R$.
6. $Q = H_1^T H_2^T \cdots H_{n-1}^T$ (o equivalentemente $Q=H_{n-1}\cdots H_1$ según convención, tal que $A=QR$).

**Propiedades del reflector $H_u=I-2uu^T$** (clásicas de examen):
- $H_u$ es simétrica y ortogonal: $H_u^T=H_u=H_u^{-1}$.
- $u$ es autovector de $H_u$ con autovalor $-1$ ($H_uu = u-2u(u^Tu)=u-2u=-u$).
- Cualquier $v\perp u$ es autovector con autovalor $+1$ ($H_uv=v$).
- $H_u$ es diagonalizable (es simétrica → diagonalizable ortogonalmente; autovalores $\{-1,1,\dots,1\}$).
- Interpretación geométrica: si $x=w+v$ con $w\parallel u$ y $v\perp u$, entonces $H_ux = -w+v$ — refleja la componente en la dirección de $u$ y deja fija la componente ortogonal.
- $\|Ax\|_2=\|QRx\|_2=\|Rx\|_2$ porque $Q$ ortogonal preserva norma-2 ⇒ $\|A\|_2=\|R\|_2$.
- Caso especial $u=e_i$: $H_{e_i}$ es la identidad con un $-1$ en la posición $i$ (refleja solo esa coordenada).

**⚠️ Por qué Householder es eficiente (la pregunta de "ganancia computacional" que suele aparecer):**

Si armás $H_k$ **explícitamente** como matriz $n\times n$ y hacés productos matriz-matriz genéricos en cada paso:
```
Por cada uno de los n pasos:
  H_k = I - 2 u_k u_kᵗ         → O(n²) construir
  A = H_k · A                  → O(n³) producto genérico
  Q = Q · H_kᵗ                 → O(n³) producto genérico
Total: n × O(n³) = O(n⁴)
```
La versión real **nunca arma $H_k$**: usa la identidad de rango 1 $A - 2u(u^TA)$ y $Q - 2(Qu)u^T$, que cuesta $O((n-k)^2)$ por paso en vez de $O(n^3)$:
```
Por cada uno de los n pasos:
  A[k:,k:] -= 2 u_k (u_kᵗ A[k:,k:])   → O((n-k)²)
  Q -= 2 (Q u_k) u_kᵗ                 → O(n²)
Total: O(n³/3) triangular + O(n³) acumular Q = O(n³)
```
**Conclusión:** la ganancia real es de un factor $n$: $O(n^3)$ en vez de $O(n^4)$ — mismo orden que LU. Si te preguntan "V/F: no hay ganancia computacional al usar Householder en vez de multiplicar matrices completas" → **Falso**, la ganancia viene de no armar $H_k$ explícitamente y usar productos de rango 1.

---

## 3. Autovalores, Diagonalización, Markov

### 3.1 Propiedades básicas de autovalores
- $A$ cuadrada: $\det(A)=\prod\lambda_i$, $\operatorname{tr}(A)=\sum\lambda_i$.
- Triangular/diagonal: autovalores = diagonal.
- $A$ inversible, $\lambda$ autovalor $\Rightarrow 1/\lambda$ autovalor de $A^{-1}$ (mismo autoespacio, mg se conserva, ma no necesariamente).
- $\lambda$ autovalor de $A$ $\Rightarrow \lambda^k$ autovalor de $A^k$.
- $A$ y $A^T$ tienen los mismos autovalores (no necesariamente los mismos autovectores).
- Proyector: únicos autovalores posibles son $0,1$.
- $e\in\ker(A^T) \iff e\perp\operatorname{col}(A)$.

### 3.2 Diagonalización — condiciones equivalentes
$A\in\mathbb K^{n\times n}$ diagonalizable $\iff$
- columnas de $A$ (autovectores) forman base $\iff$ existen $n$ autovectores LI
- $mg_A(\lambda_i)=ma_A(\lambda_i)\ \forall\lambda_i$
- $A=PDP^{-1}$ ($D$ diagonal con autovalores) $\iff A^m=PD^mP^{-1}$
- $A$ semejante a una diagonal

**Procedimiento:** hallar $\lambda_i$ → armar $D=\operatorname{diag}(\lambda_i)$ → hallar autovectores $v_i$ → $Q=[v_1|\dots|v_n]$ → $A=QDQ^{-1}$.

**Multiplicidades:**
- $mg(\lambda)=\dim(\ker(A-\lambda I))$ — triangular $A-\lambda I$ y medir el núcleo.
- $ma(\lambda)$ = multiplicidad en el polinomio característico.
- $1\leq mg(\lambda)\leq ma(\lambda)$; $\sum ma(\lambda_i)=n$.

**Diagonalizable sobre $\mathbb R$:** todos los $\lambda_i\in\mathbb R$ **y** $mg=ma$ para todos. Autovalores complejos $\Rightarrow$ no diagonalizable sobre $\mathbb R$.
- Autovalores todos distintos $\Rightarrow$ diagonalizable.
- **Toda simétrica real es diagonalizable ortogonalmente** ($A=QDQ^T$).

### 3.3 Matrices hermitianas / normales — Teorema Espectral
- $A$ hermitiana $\iff A=A^*$. Diagonal real. Autovalores reales. Autovectores de autovalores distintos son ortogonales.
- $A$ normal $\iff A^*A=AA^*$. Hermitianas $\subset$ normales. Simétrica real $\Rightarrow$ hermitiana $\Rightarrow$ normal.
- **Teorema espectral:** $A$ normal $\iff$ existe BON de autovectores $\iff A=QD Q^*$, $D$ diagonal.
- Caso particular vía Schur: si $A$ es normal, la $T$ de la descomposición de Schur ya sale diagonal.

### 3.4 Semejanza
$A,B$ semejantes $\iff \exists C: A=CBC^{-1}$. Comparten $\det$, traza, autovalores (no autovectores).

### 3.5 Markov
- Proceso: $v^{(k+1)}=Av^{(k)}$.
- Matriz de transición: $A_{ij}\geq0$, $\sum_i A_{ij}=1$ (columnas suman 1). $1$ siempre es autovalor. $|\lambda|\leq1$ para todo autovalor.
- Diagonalizable (si lo es): $v^{(k)}=A^kv^{(0)}=PD^kP^{-1}v^{(0)}$.
- **Estado de equilibrio:** $v^*$ con $Av^*=v^*$ (autovector de $\lambda=1$). Siempre existe.
- **Unicidad:** $\iff \dim(E_{\lambda=1})=1$.
- **Convergencia desde cualquier $v^{(0)}$:** $\iff \lambda=1$ es el único autovalor de módulo 1 (todos los demás $|\lambda|<1$) $\iff \exists A^\infty$.
  - Si hay otros $|\lambda|=1$: converge solo si $v^{(0)}$ no tiene componente en esos autoespacios; si tiene, oscila.
  - **Ciclos:** un subciclo determinístico de longitud $k$ genera autovalores $e^{2\pi i j/k}$ (raíces $k$-ésimas de la unidad), todos de módulo 1 — rompen la convergencia general aunque el equilibrio sea único.
- **Irreducible** (grafo fuertemente conexo) $\Rightarrow \dim(E_{\lambda=1})=1$ (equilibrio único), pero hace falta además **aperiodicidad** (sin ciclos) para converger desde cualquier $v^{(0)}$.

---

## 4. ⚠️ Schur, SVD, Pseudoinversa

### 4.1 Descomposición de Schur
**Toda matriz cuadrada** $A$ es unitariamente semejante a una triangular superior: $A=UTU^*$, $U$ unitaria ($U^*=U^{-1}$), $T$ triangular superior. Siempre existe (a diferencia de LU/Cholesky).

**Procedimiento (constructivo):**
1. Encontrar un autovector $v_1$ de $A$ (autovalor $\lambda_1$).
2. Completar una BON $Q_1$ que tenga a $v_1/\|v_1\|$ como primera columna.
3. Calcular $Q_1^*AQ_1$ — queda de la forma $\begin{pmatrix}\lambda_1 & *\\ 0 & A_2\end{pmatrix}$ ($A_2$ de tamaño $n-1$).
4. Repetir el proceso sobre la submatriz $A_2$, "pisando" (embebiendo) el resultado en las posiciones correspondientes de $U$ y $T$.

**Truco típico de examen — usar Schur para calcular $A^k$:**
Si $A=UTU^*$, descomponer $T=D+S$ con $D$ = diagonal de $T$ (los autovalores) y $S$ = parte estrictamente triangular superior de $T$ (nilpotente: $S^n=0$, y en la práctica $S^j=0$ para $j$ chico si $T$ es pequeña). Como $D$ y $S$ **no conmutan en general**, no se puede usar binomio de Newton directo salvo casos donde sí conmutan o donde se explota la nilpotencia de $S$ para truncar una serie. Revisar el ejercicio de práctica resuelto (Práctica 5, Ej 1) para el desarrollo completo caso por caso — es el ejercicio más representativo de esto.

**Usos de Schur que aparecieron en finales:**
- Splitting de un método iterativo a partir de $A=UTU^*$: $M=UD̂U^*$ (parte diagonal), $N=U T̂U^*$ (parte estrictamente triangular), y probar convergencia via $\|M^{-1}N\|_2<1$ — **practicar a fondo** (ver plan, Final 17/12/2025 Ej 3).
- Autovalores de $T$ (los mismos que $A$, en la diagonal) usados indirectamente para deducir autovalores de matrices construidas a partir de $A$ (ej. $A^TA$).

### 4.2 SVD
$A=U\Sigma V^*$, $A\in\mathbb C^{m\times n}$. Existe **para cualquier matriz**, sin restricciones.
- $U\in\mathbb C^{m\times m}$: columnas $u_j$ dadas por $Av_j=\sigma_ju_j$ (completar BON si $m>n$).
- $V\in\mathbb C^{n\times n}$: columnas = autovectores de $A^*A$.
- $\Sigma$: diagonal real $\geq0$.

**Procedimiento:**
1. Calcular $A^*A$ (o $AA^*$, lo que sea más chico).
2. Autovalores/autovectores de $A^*A$ → $V$ (autovectores normalizados), $\sigma_i=\sqrt{\lambda_i}$.
3. $U$: $u_i = Av_i/\sigma_i$ para $\sigma_i\neq0$; completar el resto con BON.

**Identidades:** $A^*A = V\Sigma^*\Sigma V^*$, $AA^* = U\Sigma\Sigma^*U^*$. Si $A$ cuadrada: $\Sigma^*\Sigma=\Sigma^2$.

**Propiedades (muy preguntadas):**
- $\|A\|_2=\sigma_1$ (el mayor). $\max_{\|x\|=1}\|Ax\|_2=\sigma_1$, $\min_{\|x\|=1}\|Ax\|_2=\sigma_n$ (si $A$ cuadrada/inyectiva).
- $A$ inversible (cuadrada, rango completo): $\operatorname{cond}_2(A)=\sigma_1/\sigma_n$.
- $A$ ortogonal $\iff$ todos los $\sigma_i=1$ ($\Sigma=I$).
- Con $r=rk(A)$: $\operatorname{Im}(A)=\langle u_1,\dots,u_r\rangle$, $\operatorname{Im}(A^T)=\langle v_1,\dots,v_r\rangle$, $\ker(A)=\langle v_{r+1},\dots,v_n\rangle$, $\ker(A^T)=\langle u_{r+1},\dots,u_m\rangle=\ker(A^TA)$.
- **Eckart-Young:** la mejor aproximación de rango $k$ en norma 2/Frobenius es $B_k=\sum_{i=1}^k\sigma_iu_iv_i^T$, con error $\|A-B_k\|_2=\sigma_{k+1}$.
- $A^*=(U\Sigma V^*)^* = V\Sigma^TU^*$ (SVD de la traspuesta/adjunta).

### 4.3 Pseudoinversa (Moore-Penrose)
$$A^+ = V\hat\Sigma^{-1}\hat U^* \quad\text{(SVD reducida, solo valores singulares no nulos)}$$
Caso columnas LI (rango completo): $A^+=(A^TA)^{-1}A^T$.
- Resuelve cuadrados mínimos: $\hat x = A^+b$ minimiza $\|Ax-b\|_2$ (ver Tema 5).
- Si $Ax=b$ tiene solución (o infinitas), $A^+b$ da la de **norma mínima**.
- 4 propiedades de Moore-Penrose (unicidad): $AA^+A=A$; $A^+AA^+=A^+$; $(AA^+)^*=AA^+$; $(A^+A)^*=A^+A$.

---

## 5. Mínimos Cuadrados, Métodos Iterativos

### 5.1 Cuadrados mínimos
Buscamos minimizar $\|Ax-b\|_2^2$ cuando $Ax=b$ no tiene solución exacta.
- **Ecuaciones normales:** $A^TAx=A^Tb \Rightarrow x=(A^TA)^{-1}A^Tb$.
- Mal condicionadas: $\operatorname{cond}(A^TA)=\operatorname{cond}_2(A)^2$ — por eso en la práctica se resuelve vía **QR** o **SVD**, no armando $A^TA$ directamente.
- $Ax=b$ tiene solución $\iff AA^+b=b \iff b\in\operatorname{col}(A)$.
- $A$ con columnas LI $\iff$ solución única de cuadrados mínimos.
- Si $Q$ tiene columnas ortonormales (de una QR de $A$): resolver $Rx=Q^Tb$ es más estable que las ecuaciones normales.

### 5.2 Métodos iterativos — marco general
$A = D+L+U$ (diagonal, triangular inferior estricta, triangular superior estricta). Forma general: $x^{(n+1)}=Tx^{(n)}+c$.

**Teorema clave:** el método converge para cualquier $x^{(0)} \iff \rho(T)<1$, con $\rho(T)=\max\{|\lambda| : \lambda \text{ autovalor de } T\}$.
- Splitting $A=M-N$, $T=-M^{-1}N$: $\lambda$ autovalor de $T$ $\iff \det(\lambda M+N)=0$.
- Converge en $n$ pasos exactos $\iff T^{n+1}=0$.
- $\rho(B) = \lim_{n\to\infty}\|B^n\|^{1/n}$ para toda norma subordinada.
- $\rho(T)\leq \|T\|$ para cualquier norma matricial inducida — **por eso puede haber $\|T\|\geq1$ y el método converger igual** (lo que importa es $\rho$, no una norma en particular). Buscar la norma en la que sí dé $<1$ es un ejercicio típico.
- **Mal condicionamiento de $A$ ≠ convergencia lenta del iterativo** — son conceptos distintos, no confundir.

### 5.3 Jacobi
$$x^{(n+1)} = -D^{-1}(L+U)x^{(n)} + D^{-1}b$$
Paralelizable (cada componente depende solo de $x^{(n)}$, no de actualizaciones parciales).
```python
def metodoJacobi(A, b, epochs):
    D, L, U = splitMatrix(A)
    D_inv = alc.inv(D)
    Bj = -D_inv @ (L + U)
    c = D_inv @ b
    x = np.random.random((n, 1))
    for i in range(epochs):
        x = Bj @ x + c
    return x
```

### 5.4 Gauss-Seidel
$$x^{(n+1)} = -(D+L)^{-1}Ux^{(n)} + (D+L)^{-1}b$$
No paralelizable (usa valores ya actualizados dentro de la misma iteración).
- Si $A$ tridiagonal: $\rho(B_{GS})=\rho(B_J)^2$ (Gauss-Seidel converge más rápido en ese caso).

### 5.5 SOR (Successive Over-Relaxation)
Generaliza Gauss-Seidel con parámetro $\omega$: matriz de iteración $B(\omega)$.
- $\det(B(\omega)) = (1-\omega)^n$ ⇒ condición **necesaria** de convergencia: $\omega\in(0,2)$.
- $\omega=1$ recupera Gauss-Seidel.

### 5.6 Método del gradiente
Para $A$ simétrica definida positiva, resolver $Ax=b$ equivale a minimizar $f(x)=\tfrac12x^TAx - b^Tx$ (forma cuadrática con mínimo único si $A$ DP).
- Dirección de máximo descenso: $-\nabla f(x) = b-Ax = $ residuo $r^{(k)}$.
- Paso óptimo en esa dirección: $t^* = \dfrac{(r^{(k)})^Tr^{(k)}}{(r^{(k)})^TAr^{(k)}}$.
- Iteración: $x^{(k+1)} = x^{(k)} + t^*r^{(k)}$.
- Converge pero con trayectoria en "zigzag" sobre las curvas de nivel (a diferencia del gradiente conjugado, que corrige esto — no cubierto con ejercicios en la guía práctica, repasar del apunte si hace falta).

### 5.7 Radio espectral — resumen
- $\rho(B)<1$: converge. $\rho(B)=1$: puede no converger o converger muy lento. $\rho(B)>1$: diverge.
- Cuanto menor $\rho(B)$, más rápida la convergencia.

---

## 6. Apéndice: errores comunes / V-F rápidas

- **Entradas positivas ≠ definida positiva.** Chequear autovalores o menores principales, no el signo de las entradas.
- **$\det(A)\to0$ no implica mal condicionamiento** (ej. matriz escalar $\varepsilon I$).
- **Mal condicionamiento ≠ convergencia lenta de métodos iterativos** — son propiedades distintas de $A$.
- **$\|T\|\geq1$ no implica que el método iterativo diverja** — lo que decide es $\rho(T)<1$, no una norma matricial particular.
- **$Ax\neq Qx$ en general** aunque $\operatorname{Im}(A)=\operatorname{Im}(Q)$ en la QR — son TL distintas.
- **Un ciclo en una cadena de Markov genera autovalores de módulo 1 ≠ 1** — el equilibrio puede ser único y aun así no converger desde todo $v^{(0)}$.
- **No toda matriz de Markov es diagonalizable** — no asumir diagonalización sin chequear $mg=ma$.
- **Householder sin la optimización de rango 1 es $O(n^4)$, no $O(n^3)$** — la ganancia viene de no armar $H_k$ explícitamente.
- **LU no existe para toda matriz** (Schur y SVD sí); un cero en la diagonal durante la eliminación no siempre mata la factorización, pero un cero en $a_{11}$ sí.
- **Cholesky y LU son las únicas descomposiciones "no universales"** de la tabla del punto 0 — todo lo demás (QR, Schur, SVD) existe siempre bajo condiciones mínimas de forma.
