# Resumen Final — Álgebra Lineal Computacional

Unifica `resumen-1er-parcial.md` y `resumen-2do-parcial.md`, reorganizado según los 5 bloques temáticos del final (28/07). Formato pensado para **repaso rápido**: cada bullet empieza con la propiedad/fórmula en negrita (lo que hay que recordar); las líneas `·` debajo son la justificación o un ejemplo — se pueden saltear si ya te acordás del resultado.

## Índice
0. [Tabla rápida: existencia y costo de descomposiciones](#0-tabla-rápida-existencia-y-costo-de-descomposiciones)
1. [Subespacios, Bases, TL, Normas, Número de Condición](#1-subespacios-bases-tl-normas-número-de-condición)
2. [LU, Cholesky, Ortogonalidad, Proyectores, QR ( Householder)](#2-lu-cholesky-ortogonalidad-proyectores-qr)
3. [Autovalores, Diagonalización, Markov](#3-autovalores-diagonalización-markov)
4. [ Schur, SVD, Pseudoinversa](#4-schur-svd-pseudoinversa)
5. [Mínimos Cuadrados, Métodos Iterativos](#5-mínimos-cuadrados-métodos-iterativos)
6. [Apéndice: errores comunes / V-F rápidas](#6-apéndice-errores-comunes--v-f-rápidas)
   - 6.1 [Numpy — funciones útiles](#61-numpy--funciones-útiles)

---

## 0. Tabla rápida: existencia y costo de descomposiciones

| Descomposición | Forma | Requisito | ¿Existe siempre? | Costo | Notas |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **LU** | $A=LU$ | Cuadrada, $\det(A_k)\neq0$ para $k=1..n-1$ | NO ❌ (sí con $PA=LU$) | $O(n^3/3)$ | Única si $A$ inversible |
| **Cholesky** | $A=\hat L\hat L^T$ | Simétrica **definida positiva** | NO ❌ | $O(n^3/3)$, mitad de flops que LU | La más restrictiva; única si existe |
| **QR** | $A=QR$ | $m\geq n$, cualquier $A$ | SÍ ✅ | Householder: $O(n^3)$ (naive $O(n^4)$) · GS: $O(n^3)$, menos estable | Única si columnas LI o $R_{ii}>0$ |
| **Schur** | $A=UTU^*$ | Cuadrada | SÍ ✅ | $O(n^3)$ (QR iterativo, ${\sim}10n^3$ en la práctica) | $T$ triangular superior, $U$ unitaria |
| **SVD** | $A=U\Sigma V^*$ | Ninguno — **cualquier matriz** | SÍ ✅ | $O(mn\cdot\min(m,n))$; $O(n^3)$ si cuadrada | La más general y la más cara |
| **Inversa** ($A^{-1}$) | — | Cuadrada, $\det\neq0$ | NO ❌ | $O(n^3)$ (mismo orden que LU, más constante) | Casi nunca hace falta calcularla |

**Costo de multiplicar** $(p\times q)\cdot(q\times r)$: $O(pqr)$.

**Costos básicos:** matriz-vector $O(n^2)$ · matriz-matriz $O(n^3)$ · sistema **triangular** (sustitución) $O(n^2)$ · $Ax=b$ genérico vía LU/QR: $O(n^3)$ factorizar + $O(n^2)$ por cada $b$ adicional reusando la factorización.

**"Resolver, no invertir"**
- $A^{-1}$ explícita cuesta $O(n^3)$ (mismo orden que resolver, constante mayor — equivale a resolver $n$ sistemas $Ax_i=e_i$).
- **Nunca** $x=A^{-1}b$ para resolver $Ax=b$: más lento y menos estable que LU/QR + sustitución.
- Muchos $b$ con la misma $A$: factorizar una vez ($O(n^3)$) + resolver cada uno ($O(n^2)$ c/u) sigue siendo más estable que invertir, incluso cuando el costo total se empareja ($\geq n$ lados derechos).
- $(A^TA)^{-1}$ (ecuaciones normales) es doblemente caro: $O(mn^2)$ armar $A^TA$ + $O(n^3)$ invertir, y $\operatorname{cond}(A^TA)=\operatorname{cond}_2(A)^2$ → preferir QR/SVD (ver 5.1).

---

## 1. Subespacios, Bases, TL, Normas, Número de Condición

### 1.1 Dimensión y subespacios
- **$Bx$ = combinación lineal de columnas:** $B=[b_1|\dots|b_n]\Rightarrow Bx=\sum x_ib_i$.
  · De $x=\sum x_ie_i$, $Be_i=b_i$, linealidad. Es *la* forma correcta de pensar $Bx$ (no "fila por columna"). Consecuencia: $Bx=0$ no trivial $\iff$ columnas LD.
- **Dominio/codominio:** $f:\mathbb V\to\mathbb W$, $\mathbb V$=dominio, $\mathbb W$=codominio declarado; $\operatorname{Im}(f)\subseteq\mathbb W$ es lo que $f$ realmente alcanza.
- **Teorema rango-nulidad:** $\dim(\operatorname{Im} f)+\dim(\ker f)=\dim(\mathbb V)$ — vale **siempre**, para toda TL (mono/epi/iso se derivan de acá, no al revés).
  · Ej: $f(x,y,z)=(x,y,0)$ en $\mathbb R^3$: ni mono ni epi, $\dim\ker=1$, $\dim\operatorname{Im}=2$, $1+2=3$ igual.
  ·  Reparte la dimensión del **dominio**, no del codominio: $A\in\mathbb R^{m\times n}\Rightarrow$ TL $\mathbb R^n\to\mathbb R^m$, la suma da $n$ (columnas), no $m$.
  · Consecuencia (cuadrados mínimos): $A\in\mathbb R^{m\times n}$, $m\geq n$, $rg(A)=n\Rightarrow\ker A=\{0\}\Rightarrow$ columnas LI $\Rightarrow A^TA$ inversible.
- **Mono / epi / iso:**

  | | condición |
  | :--- | :--- |
  | Mono (inyectiva) | $\ker(f)=\{0\}$ |
  | Epi (sobreyectiva) | $\operatorname{Im}(f)=\mathbb W$ |
  | Iso (biyectiva) | mono + epi $\Rightarrow\dim\mathbb V=\dim\mathbb W$ |

  · Si $\dim\mathbb V=\dim\mathbb W$ (dim. finita): alcanza probar una sola (mono $\iff$ epi $\iff$ iso).
  ·  Codominio ≠ imagen: todo elemento de $\operatorname{Im}(f)$ tiene preimagen por definición (trivial); epi es cuando $\operatorname{Im}(f)=\mathbb W$ **completo**. Ej: $f:\mathbb R^2\to\mathbb R^3$, $(x,y)\mapsto(x,y,0)$ — imagen es un plano $\subsetneq\mathbb R^3$, no epi.
- $\dim(A+B)=\dim(A)+\dim(B)-\dim(A\cap B)$.
- $rk(A)=\dim(\operatorname{Im}(A))$, $\ker(A)=\{x:Ax=0\}$.
- **Suma:** $S+T=\langle$ generadores de $S$ y $T$ juntos $\rangle$.
- **Intersección:** resolver $S\alpha-T\beta=0$ (sistema homogéneo, $[S\,|\!-T]$); los $S\alpha$ resultantes son base de $S\cap T$.
- **Suma directa:** $L\oplus M\iff L\cap M=\{0\}\ \land\ L+M=\mathbb K^n$.
- $S^\perp=\{v:\langle v,s\rangle=0\ \forall s\in S\}$.
- **Base desde generadores:** ponerlos como filas, triangular, descartar filas que dan 0 o son CL de otras.

### 1.2 Inversibilidad
$A$ **singular** $\iff\lambda=0$ autovalor $\iff\det(A)=0\iff\dim(\ker A)>0\iff$ columnas/filas LD.
$A$ **inversible** $\iff\det(A)\neq0\iff$ columnas/filas LI $\iff$ columnas base de $\mathbb K^n$.
- Diagonal dominante estricta $\implies$ inversible.

### 1.3 Transformaciones lineales y cambio de base
- **Extraer base de generadores:** (1) filas + triangular, base sale del sistema resultante; o (2) columnas + triangular, las columnas con pivote (≠0 en la diagonal) son las que forman base entre los vectores originales.
- **Extender una base:** agregar candidatos (ej. canónicos $e_i$) como fila al final, triangular contra los pivotes existentes; si la fila se anula, es redundante.
  ·  No alcanza con "mirar" la fila agregada sin terminar de reducirla contra todos los pivotes previos en su columna — puede requerir pivoteo (reordenar filas) para revelarlo.
- $C_{BE}\,[(x,y,z)]_B=[(x,y,z)]_E\iff[(x,y,z)]_B=C_{EB}\,[(x,y,z)]_E$.
- **$f$ es TL** (dos formas equivalentes, cualquiera alcanza):
  - $f(u+v)=f(u)+f(v)$ **y** $f(\alpha u)=\alpha f(u)$, ó
  - $f(\alpha u+v)=\alpha f(u)+f(v)$ (más rápida de chequear).
  · Consecuencias gratis (no hace falta probarlas): $f(0)=0$, $f(-v)=-f(v)$.
- $[f]_{BB'}$: columnas = coordenadas de $f(v_i)$ en base $B'$. $[f(v)]_{B'}=[f]_{BB'}[v]_B$.
- $C_{BB'}$: columnas = vectores de $B$ en $B'$. $C_{B'B}=(C_{BB'})^{-1}$. $C_{BB''}=C_{B'B''}C_{BB'}$.
- Práctico: $C_{BE}$ = vectores de $B$ como columnas; $C_{EB}=(C_{BE})^{-1}$ (si $B$ ortonormal, $=C_{BE}^T$).
- $[f]_{BE}=[f]_{EE}\,C_{BE}$.
- **Clave:** dos matrices representan la misma TL en bases distintas $\iff$ **semejantes** $\iff\exists C:A=CBC^{-1}$ $\implies$ comparten det, traza, autovalores.

### 1.4 Normas vectoriales y matriciales
- **Axiomas norma vectorial:** (1) $\|v\|\geq0$, $=0\iff v=0$; (2) $\|\alpha v\|=|\alpha|\|v\|$; (3) $\|u+v\|\leq\|u\|+\|v\|$.
- **Axiomas norma matricial:** los mismos 3 (viendo $A$ como vector) **+ submultiplicatividad** $\|AB\|\leq\|A\|\|B\|$.
  · No sale gratis de los primeros 3 — ej. norma del máximo de entradas cumple los 3 pero no submult., no es norma matricial. Toda norma **inducida** ($\|A\|=\max_{\|x\|=1}\|Ax\|$) sí es automáticamente submultiplicativa.
- $\|v\|_1=\sum|v_i|$, $\|v\|_\infty=\max|v_i|$, $\|v\|_p=(\sum|v_i|^p)^{1/p}$.
- **Cauchy-Schwarz:** $|x^*y|\leq\|x\|\|y\|$. Sin valor absoluto ($\|\cdot\|_2$): $\langle u,v\rangle\leq\|u\|_2\|v\|_2$.
- $\|A\|_\infty$ = máx suma de **filas** (mód.). $\|A\|_1$ = máx suma de **columnas**.
  · Prueba: $|(Ax)_i|=|\sum_ja_{ij}x_j|\leq\sum_j|a_{ij}||x_j|$ da la cota superior; la inferior sale eligiendo $x_j=\operatorname{sign}(a_{pj})$.
- $\|A\|_2=\sigma_1(A)$ (mayor valor singular, ver Tema 4).
- $\|Ax\|\leq\|A\|\|x\|$, $\|AB\|\leq\|A\|\|B\|$.
- Relaciones: $\tfrac1{\sqrt n}\|A\|_\infty<\|A\|_2<\sqrt n\|A\|_\infty$ · $\tfrac1{\sqrt n}\|A\|_2<\|A\|_1<\sqrt n\|A\|_2$ · $\|A\|_\infty<\|A\|_1<n\|A\|_\infty$.

### 1.5 Número de condición
- $\kappa_*(A)=\|A\|_*\|A^{-1}\|_*$, rango $[1,\infty)$. Grande $\implies$ mal condicionado, sensible a errores de redondeo.
- $\dfrac1{\kappa(A)}\leq\inf_{H\text{ sing.}}\dfrac{\|A-H\|}{\|A\|}$ — acota la distancia relativa a una matriz singular.
- Error relativo: $\dfrac{\|x-\tilde x\|}{\|x\|}\leq\kappa(A)\dfrac{\|b-\tilde b\|}{\|b\|}$.
-  $\det(A)\to0$ **no implica** mal condicionamiento (ej: $D_n=\varepsilon I$, $\det\to0$ pero $\kappa=1$).
-  Mal condicionamiento ≠ convergencia lenta de iterativos — independientes (ver Tema 5).

---

## 2. LU, Cholesky, Ortogonalidad, Proyectores, QR

### 2.1 LU
$A=LU$, $L$ triangular inferior (diag 1), $U$ triangular superior.
- **Existencia (sin pivoteo):** todos los pivotes $\neq0$ $\iff\det(A_k)\neq0$ para $k=1,\dots,n-1$ (menores principales top-left, no hace falta para $k=n$).
- Cero en $a_{11}$ mata la existencia; cero en otra posición de la diagonal no implica que algún menor sea 0.
- Única si $A$ inversible.
- **Algoritmo:** triangular con operaciones elementales → $U$; en $L_{ij}$ guardar el inverso (con signo) del multiplicador que anuló la celda $ij$.
- **$PA=LU$:** existe siempre para cuadrada, permutando filas primero.

### 2.2 Cholesky
$A=\hat L\hat L^T$. Existe $\iff A$ **simétrica y definida positiva**. Única si existe.
**Algoritmo:** (1) $A=LU$ · (2) $D_{ii}=U_{ii}$, resto 0 · (3) $D_1=\sqrt D$ (entrada a entrada) · (4) $\hat L=LD_1$ · (5) $A=\hat L\hat L^T$.

### 2.3 Definida positiva — condiciones equivalentes
$A$ DP $\iff$ $x^TAx>0\ \forall x\neq0$ $\iff$ todos $\lambda_i>0$ $\iff$ $A=LU$ con $U_{ii}>0\ \forall i$ $\iff$ todos los menores principales $\det(A_k)>0$.

 Entradas positivas ≠ definida positiva. Contraejemplo: $\begin{bmatrix}1&2\\2&1\end{bmatrix}$, $\det=-3<0$.
Cadena: $\exists\det(A_k)\leq0\iff\exists\lambda_i\leq0\iff$ no DP $\iff$ no tiene Cholesky. Semidefinida positiva (SDP): $\lambda_i\geq0$.

### 2.4 Proyectores (generales)
- $P$ proyector $\iff P^2=P$. Autovalores solo $0$ o $1$.
- $\operatorname{Im}(P)=\operatorname{Nu}(I-P)$, $\operatorname{Im}(I-P)=\operatorname{Nu}(P)$ (en la imagen, $Px=x$).
- $V=\operatorname{Im}(P)\oplus\operatorname{Nu}(P)$ siempre.
- **Armado:** $M=[v_1|\dots|v_n]$, quiero $P(v_i)$=columnas de $R$: $PM=R\Rightarrow P=RM^{-1}$.

### 2.5 Proyección ortogonal
- $P_S$ ortogonal $\iff P_S^2=P_S\land P_S^T=P_S$. $\operatorname{Nu}(P)\perp\operatorname{Im}(P)$.
- $x=P_S(x)+P_{S^\perp}(x)$, única; $x-P_S(x)\in S^\perp$ (error).
- $A$ con columnas LI (no ortonorm.): $P_S=A(A^TA)^{-1}A^T$.
- $Q$ con columnas ortonormales ($Q^TQ=I$): $P_S=QQ^T$.
- Desde BON $\{v_1,\dots,v_r\}$: $P_S=\sum v_iv_i^T$; $P_S(x)=\sum(v_i^Tx)v_i$.
- Complementario: $P_{S^\perp}=I-P_S$. $\operatorname{Im}(P_S)=S$, $\operatorname{Nu}(P_S)=S^\perp$.
- $u_iu_i^T$ siempre rango 1. $A=uv^T\Rightarrow\operatorname{Im}(A)=\langle u\rangle$, $\operatorname{Nu}(A)=\langle v\rangle^\perp$.

### 2.6 Gram-Schmidt
$a_1,\dots,a_n$ columnas de $A$. $v_1=a_1$; $v_{i+1}=a_{i+1}-\sum_{j\leq i}p_{v_j}(a_{i+1})$, con $p_a(b)=\frac{a^Tb}{a^Ta}a$. Normalizar: $v_i'=v_i/\|v_i\|$.
- $Q$ ortonormal $\Rightarrow Q^T=Q^{-1}$, $\|Qx\|=\|x\|$.

### 2.7 QR
Existe para toda $A\in\mathbb R^{m\times n}$, $m\geq n$ (rango completo; si no, algunas filas de $R$ nulas). Única si columnas LI o $R_{ii}>0$.
- $\operatorname{Im}(A)=\operatorname{Im}(Q)$ pero **$Ax\neq Qx$ en general** (TL distintas).
- **Vía Gram-Schmidt:** $Q=[q_1|\dots|q_n]$, $q_i=v_i/\|v_i\|$. $R=Q^TA$, triangular superior, $r_{kk}=\|v_k\|$.

### 2.8  QR vía Householder

**Idea:** en cada paso, reflejar la columna actual sobre $e_1$ del subespacio restante con una matriz de reflexión, hasta triangular $A$ por completo.

**Algoritmo:**
1. $v_k$ = columna $k$ de $A[k{:}n,\,k{:}n]$. $w_k=(\|v_k\|,0,\dots,0)^T$.
2. $u_k=\dfrac{v_k-w_k}{\|v_k-w_k\|}$.
3. $H_k=I-2u_ku_k^T$ (reflector, actúa solo sobre el bloque restante).
4. $A\leftarrow H_kA$ (anula debajo de la diagonal en la columna $k$).
5. Repetir hasta $A$ triangular superior $=R$.
6. $Q=H_1^T H_2^T\cdots H_{n-1}^T$ (equiv. $Q=H_{n-1}\cdots H_1$ según convención, tal que $A=QR$).

**Propiedades del reflector $H_u=I-2uu^T$** (clásicas de examen):
- Simétrica y ortogonal: $H_u^T=H_u=H_u^{-1}$.
- $u$ es autovector con autovalor $-1$ ($H_uu=u-2u=-u$); cualquier $v\perp u$ es autovector con autovalor $+1$.
- Diagonalizable (simétrica → diagonalizable ortogonalmente; autovalores $\{-1,1,\dots,1\}$).
- Geométricamente: $x=w+v$ ($w\parallel u$, $v\perp u$) $\Rightarrow H_ux=-w+v$ — refleja la componente en $u$, deja fija la ortogonal.
- $\|Ax\|_2=\|QRx\|_2=\|Rx\|_2$ ($Q$ preserva norma-2) $\Rightarrow\|A\|_2=\|R\|_2$.
- Caso $u=e_i$: $H_{e_i}=I$ con un $-1$ en la posición $i$.

** Por qué Householder es eficiente** (pregunta típica de "ganancia computacional"):

| Enfoque | Por paso | Total ($n$ pasos) |
| :--- | :--- | :--- |
| Naive: armar $H_k$ explícita ($n\times n$) y hacer $A=H_kA$, $Q=QH_k^T$ como productos genéricos | $O(n^2)$ construir $H_k$ + $O(n^3)$ por cada producto genérico | $O(n^4)$ |
| Real: identidad de rango 1, $A-2u(u^TA)$ y $Q-2(Qu)u^T$, sin armar $H_k$ | $O((n-k)^2)$ triangular + $O(n^2)$ actualizar $Q$ | $O(n^3)$ |

**Conclusión:** la ganancia es un factor $n$: $O(n^3)$ en vez de $O(n^4)$ — mismo orden que LU. V/F "no hay ganancia al usar Householder en vez de multiplicar matrices completas" → **Falso**: la ganancia viene de no armar $H_k$ y usar productos de rango 1.

### 2.9 Resolver $Ax=b$ con cada descomposición

Una vez factorizada $A$, resolver es **sustitución triangular** ($O(n^2)$, ver 0). La receta cambia según la descomposición:

**LU** ($PA=LU$): factorizar → $Ly=Pb$ (progresiva) → $Ux=y$ (regresiva).
```python
P, L, U = scipy.linalg.lu(A)
y = solve_triangular(L, P @ b, lower=True)
x = solve_triangular(U, y, lower=False)
```

**Cholesky** ($A=\hat L\hat L^T$, requiere sim. DP): factorizar → $\hat Ly=b$ (progresiva) → $\hat L^Tx=y$ (regresiva). Mismo esquema que LU, mitad de operaciones (ver 2.2).
```python
Lhat = np.linalg.cholesky(A)
y = solve_triangular(Lhat, b, lower=True)
x = solve_triangular(Lhat.T, y, lower=False)
```

**QR** ($A=QR$, cuadrada rango completo): factorizar → $y=Q^Tb$ ($Q^{-1}=Q^T$, nunca invertir) → $Rx=y$ (regresiva). Más caro que LU pero más estable — preferir si $A$ mal condicionada (ver 5.1).
```python
Q, R = np.linalg.qr(A)
y = Q.T @ b
x = solve_triangular(R, y, lower=False)
```

**SVD** ($A=U\Sigma V^T$, sirve con $A$ singular o rectangular): si $A$ inversible, $x=V\Sigma^{-1}U^Tb$; en general $x=A^+b=\hat V\hat\Sigma^{-1}\hat U^Tb$ (pseudoinversa, ver 4.3) — exacta si existe, o de norma mínima entre las de cuadrados mínimos si no.
```python
U, s, Vt = np.linalg.svd(A)
x = Vt.T @ ((U.T @ b) / s)   # válido si ningún s_i es 0
```

**Cuál usar:**

| Situación | Descomposición |
| :--- | :--- |
| $A$ cuadrada, bien condicionada, sin estructura | LU |
| $A$ simétrica definida positiva | Cholesky (mitad de costo que LU) |
| $A$ mal condicionada / importa estabilidad | QR |
| $A$ singular, rectangular, o rango incompleto | SVD (pseudoinversa) |

**Nunca:** $A^{-1}$ explícita + $x=A^{-1}b$ (ver 0, "resolver, no invertir").

---

## 3. Autovalores, Diagonalización, Markov

### 3.1 Propiedades básicas de autovalores
- $\det(A)=\prod\lambda_i$, $\operatorname{tr}(A)=\sum\lambda_i$.
- Triangular/diagonal: autovalores = diagonal.
- $A$ inversible, $\lambda$ autovalor $\Rightarrow1/\lambda$ autovalor de $A^{-1}$ (mismo autoespacio; mg se conserva, ma no necesariamente).
- $\lambda$ autovalor de $A\Rightarrow\lambda^k$ autovalor de $A^k$.
- $A$ y $A^T$: mismos autovalores (no necesariamente los mismos autovectores).
- Proyector: únicos autovalores posibles $0,1$.
- $e\in\ker(A^T)\iff e\perp\operatorname{col}(A)$.

### 3.2 Diagonalización — condiciones equivalentes
$A\in\mathbb K^{n\times n}$ diagonalizable $\iff$
- columnas (autovectores) forman base $\iff$ existen $n$ autovectores LI
- $mg_A(\lambda_i)=ma_A(\lambda_i)\ \forall\lambda_i$
- $A=PDP^{-1}$ ($D$ diagonal de autovalores) $\iff A^m=PD^mP^{-1}$
- $A$ semejante a una diagonal

**Procedimiento:** $\lambda_i\to D=\operatorname{diag}(\lambda_i)\to$ autovectores $v_i\to Q=[v_1|\dots|v_n]\to A=QDQ^{-1}$.

**Multiplicidades:** $mg(\lambda)=\dim(\ker(A-\lambda I))$ (triangular $A-\lambda I$, medir núcleo); $ma(\lambda)$=mult. en el pol. característico; $1\leq mg\leq ma$; $\sum ma(\lambda_i)=n$.

**Diagonalizable sobre $\mathbb R$:** todos $\lambda_i\in\mathbb R$ **y** $mg=ma$ para todos (complejos $\Rightarrow$ no diagonalizable sobre $\mathbb R$).
- Autovalores todos distintos $\Rightarrow$ diagonalizable.
- **Toda simétrica real es diagonalizable ortogonalmente** ($A=QDQ^T$).

### 3.3 Matrices hermitianas / normales — Teorema Espectral
- Hermitiana ($A=A^*$): diagonal real, autovalores reales, autovectores de autovalores distintos son ortogonales.
- Normal ($A^*A=AA^*$): hermitianas $\subset$ normales; simétrica real $\Rightarrow$ hermitiana $\Rightarrow$ normal.
- **Teorema espectral:** $A$ normal $\iff\exists$ BON de autovectores $\iff A=QDQ^*$, $D$ diagonal.
- Vía Schur: si $A$ normal, la $T$ de Schur ya sale diagonal.

### 3.4 Semejanza
$A,B$ semejantes $\iff\exists C:A=CBC^{-1}$. Comparten $\det$, traza, autovalores (no autovectores).

### 3.5 Markov
- Proceso: $v^{(k+1)}=Av^{(k)}$.
- Matriz de transición: $A_{ij}\geq0$, $\sum_iA_{ij}=1$ (columnas suman 1) $\Rightarrow\lambda=1$ siempre es autovalor, $|\lambda|\leq1\ \forall\lambda$.
- Diagonalizable: $v^{(k)}=A^kv^{(0)}=PD^kP^{-1}v^{(0)}$.
- **Estado de equilibrio** $v^*$: $Av^*=v^*$ (autovector de $\lambda=1$). Siempre existe.
- **Unicidad** $\iff\dim(E_{\lambda=1})=1$.
- **Convergencia desde cualquier $v^{(0)}$** $\iff\lambda=1$ es el único autovalor de módulo 1 (resto $|\lambda|<1$) $\iff\exists A^\infty$.
  · Si hay otros $|\lambda|=1$: converge solo si $v^{(0)}$ no tiene componente en esos autoespacios; si tiene, oscila.
  · **Ciclos:** un subciclo determinístico de longitud $k$ genera autovalores $e^{2\pi ij/k}$ (raíces $k$-ésimas de 1), todos de módulo 1 — rompen la convergencia general aunque el equilibrio sea único.
- **Irreducible** (grafo fuertemente conexo) $\Rightarrow\dim(E_{\lambda=1})=1$ (equilibrio único); hace falta además **aperiodicidad** para converger desde cualquier $v^{(0)}$.

---

## 4.  Schur, SVD, Pseudoinversa

### 4.1 Descomposición de Schur
**Toda matriz cuadrada** es unitariamente semejante a triangular superior: $A=UTU^*$, $U$ unitaria, $T$ triangular superior. Siempre existe (a diferencia de LU/Cholesky).

**Procedimiento (constructivo):**
1. Autovector $v_1$ de $A$ (autovalor $\lambda_1$).
2. Completar BON $Q_1$ con $v_1/\|v_1\|$ como primera columna.
3. $Q_1^*AQ_1=\begin{pmatrix}\lambda_1&*\\0&A_2\end{pmatrix}$ ($A_2$ de tamaño $n-1$).
4. Repetir sobre $A_2$, embebiendo el resultado en $U$ y $T$.

**Truco de examen — $A^k$ vía Schur:** $A=UTU^*$, $T=D+S$ ($D$=diagonal de autovalores, $S$=parte estrictamente triangular superior, nilpotente: $S^n=0$). $D$ y $S$ **no conmutan en general** → no vale binomio de Newton directo salvo que conmuten o se explote la nilpotencia para truncar. (Ver Práctica 5, Ej 1 para el desarrollo caso por caso.)

**Usos que aparecieron en finales:**
- Splitting de un iterativo desde $A=UTU^*$: $M=UD̂U^*$ (diagonal), $N=UT̂U^*$ (estrictamente triangular), probar convergencia vía $\|M^{-1}N\|_2<1$ — practicar (Final 17/12/2025 Ej 3).
- Autovalores de $T$ (= los de $A$, en la diagonal) usados para deducir autovalores de matrices derivadas (ej. $A^TA$).

### 4.2 SVD
$A=U\Sigma V^*$, existe **para cualquier matriz** $\mathbb C^{m\times n}$, sin restricciones.
- $U\in\mathbb C^{m\times m}$: $u_j$ dados por $Av_j=\sigma_ju_j$ (completar BON si $m>n$).
- $V\in\mathbb C^{n\times n}$: autovectores de $A^*A$.
- $\Sigma$: diagonal real $\geq0$.

**Procedimiento:** (1) $A^*A$ (o $AA^*$, el más chico) → (2) autovalores/autovectores → $V$, $\sigma_i=\sqrt{\lambda_i}$ → (3) $u_i=Av_i/\sigma_i$ ($\sigma_i\neq0$), completar resto con BON.

**Identidades:** $A^*A=V\Sigma^*\Sigma V^*$, $AA^*=U\Sigma\Sigma^*U^*$. Si $A$ cuadrada: $\Sigma^*\Sigma=\Sigma^2$.

**Propiedades (muy preguntadas):**
- $\|A\|_2=\sigma_1$; $\max_{\|x\|=1}\|Ax\|_2=\sigma_1$, $\min_{\|x\|=1}\|Ax\|_2=\sigma_n$ (si cuadrada/inyectiva).
- $A$ inversible cuadrada: $\operatorname{cond}_2(A)=\sigma_1/\sigma_n$.
- $A$ ortogonal $\iff$ todos $\sigma_i=1$ ($\Sigma=I$).
- Con $r=rk(A)$: $\operatorname{Im}(A)=\langle u_1,\dots,u_r\rangle$, $\operatorname{Im}(A^T)=\langle v_1,\dots,v_r\rangle$, $\ker(A)=\langle v_{r+1},\dots,v_n\rangle$, $\ker(A^T)=\langle u_{r+1},\dots,u_m\rangle=\ker(A^TA)$.
- **Eckart-Young:** mejor aproximación de rango $k$ (norma 2/Frobenius): $B_k=\sum_{i=1}^k\sigma_iu_iv_i^T$, error $\|A-B_k\|_2=\sigma_{k+1}$.
- $A^*=(U\Sigma V^*)^*=V\Sigma^TU^*$.

### 4.3 Pseudoinversa (Moore-Penrose)
$$A^+=V\hat\Sigma^{-1}\hat U^*\quad\text{(SVD reducida, solo }\sigma_i\neq0\text{)}$$
Rango completo de columnas: $A^+=(A^TA)^{-1}A^T$.
- Resuelve cuadrados mínimos: $\hat x=A^+b$ minimiza $\|Ax-b\|_2$ (ver Tema 5).
- Si $Ax=b$ tiene solución (o infinitas), $A^+b$ da la de **norma mínima**.
- 4 propiedades (unicidad): $AA^+A=A$; $A^+AA^+=A^+$; $(AA^+)^*=AA^+$; $(A^+A)^*=A^+A$.

---

## 5. Mínimos Cuadrados, Métodos Iterativos

### 5.1 Cuadrados mínimos

**Planteo:** $A\in\mathbb R^{m\times n}$, $b\in\mathbb R^m$; cuando $Ax=b$ no tiene solución exacta ($b\notin\operatorname{Col}(A)$), buscar $\hat x$ que minimice $\|Ax-b\|_2^2$ (MSE).

**① Existencia — siempre**, sin condición sobre $A$ o $b$: $A\hat x$ debe ser la **proyección ortogonal** de $b$ sobre $\operatorname{Col}(A)$, que siempre existe y es única *como vector* para cualquier subespacio. Lo que puede no ser único es el $\hat x$ que la produce (③).

**② Ecuaciones normales** ($\hat x$ minimizador $\iff A^TA\hat x=A^Tb$):
1. $A\hat x$=proyección de $b$ $\iff b-A\hat x\perp\operatorname{Col}(A)$ (propiedad que define proyección ortogonal).
2. $\operatorname{Col}(A)=\operatorname{gen}\{a_1,\dots,a_n\}$; ⊥ a todo el subespacio $\iff$ ⊥ a cada generador: $a_i^T(b-A\hat x)=0\ \forall i$.
3. Apilar esas $n$ ecuaciones = multiplicar por $A^T$: $A^T(b-A\hat x)=0$.
4. $\Rightarrow A^TA\hat x=A^Tb$.

 $\operatorname{Col}(A)=\{Ax:x\}$ (columnas) ≠ $\{A^Ty:y\}$ (espacio fila; coinciden solo si $A$ simétrica). No es "$b-Ax\perp A^T$" (una matriz no es subespacio): es $b-Ax\perp\operatorname{Col}(A)$, y $A^T(b-Ax)=0$ es su forma matricial.

**③ Unicidad — depende del rango:** todo $\hat x$ que cumple las ecs. normales da el mismo $A\hat x$, pero $\hat x$ es único $\iff$ columnas de $A$ LI $\iff rg(A)=n$ $\iff A^TA$ inversible. Ahí: $\hat x=(A^TA)^{-1}A^Tb$.
- Si $rg(A)<n$: infinitos $\hat x$ minimizan (misma $A\hat x$). El de **norma mínima** es $\hat x=A^+b$ (pseudoinversa, 4.3) — por eso $A^+$ generaliza $(A^TA)^{-1}A^T$ también sin rango completo.

**En la práctica:**
- Mal condicionadas: $\operatorname{cond}(A^TA)=\operatorname{cond}_2(A)^2$ → resolver vía **QR** o **SVD**, no armando $A^TA$.
- $Ax=b$ tiene solución exacta $\iff AA^+b=b\iff b\in\operatorname{Col}(A)$ (error mínimo = 0).
- $Q$ ortonormal (de QR de $A$): resolver $Rx=Q^Tb$ es más estable que las ecuaciones normales.

### 5.2 Interpolación polinomial — matriz de Vandermonde

Dados $n+1$ puntos $(x_0,y_0),\dots,(x_n,y_n)$ con $x_i$ **distintos**: buscar $p$ de grado $\leq n$ tal que $p(x_i)=y_i\ \forall i$ (interpola **exactamente**).

$p(x)=a_0+a_1x+\dots+a_nx^n$ + condiciones $p(x_i)=y_i$ → sistema $Va=y$ con
$$V=\begin{pmatrix}1&x_0&x_0^2&\cdots&x_0^n\\1&x_1&x_1^2&\cdots&x_1^n\\\vdots&&&&\vdots\\1&x_n&x_n^2&\cdots&x_n^n\end{pmatrix}\in\mathbb R^{(n+1)\times(n+1)}$$

- **Cuadrada:** $n+1$ ecuaciones para $n+1$ incógnitas (coeficientes).
- $\det(V)=\prod_{i<j}(x_j-x_i)\neq0\iff x_i$ distintos $\Rightarrow V$ inversible (bajo la hipótesis del problema) $\Rightarrow$ polinomio interpolador existe y es **único**.
  · Prueba (inducción en $n$): fijando $x_0,\dots,x_{n-1}$, $\det(V_n)$ es polinomio en $x_n$ de grado $n$ (coef. principal, columna $x_n^n$, única vez). Si $x_n=x_i$ ($i<n$) hay dos filas iguales $\Rightarrow\det=0$ → cada $(x_n-x_i)$ es raíz → $\det(V_n)=c\prod_{i<n}(x_n-x_i)$, y por cofactores $c=\det(V_{n-1})$. Repitiendo (inducción) da el producto telescópico sobre todos los pares $i<j$. Caso base $n=1$: $\det\begin{pmatrix}1&x_0\\1&x_1\end{pmatrix}=x_1-x_0$. ✓
  ·  No es hipótesis extra: la única hipótesis es "$x_i$ distintos" (dato del enunciado); la invertibilidad de $V$ **se deduce**, no se asume aparte (coherente con 1.2: filas LD $\iff\det=0$).
- **Es la misma matriz que en cuadrados mínimos con polinomios** (Ej 3) — lo que distingue "ajuste" de "interpolación" es la **forma** de $V$:
  - Cuadrados mínimos: $m$ puntos, grado $d$ con $d+1<m$ → $V$ rectangular ($m\times(d+1)$, sobredeterminado) → sin solución exacta en general, resolver por ecs. normales.
  - Interpolación: grado = puntos$-1$ → $V$ cuadrada → si inversible, solución única con **residuo cero**.
  - Interpolar = caso particular/degenerado de cuadrados mínimos con $V$ cuadrada e inversible (Práctica 6, Ej 6) — ahí minimizar y resolver exactamente coinciden (mínimo = 0).
  - Ejemplo (Práctica 6, Ej 3): 4 puntos, grado 3 → 4 coef. para 4 puntos → $V$ cuadrada → interpola exacto (a diferencia de grado 1 o 2, donde $V$ es $4\times2$/$4\times3$ rectangular y hay error).
- **Resultado final:** resolver $Va=y$ (Gauss/LU, nunca $V^{-1}$ explícita) da $a=(a_0,\dots,a_n)$ → reemplazar en $p(x)=a_0+a_1x+\dots+a_nx^n$.
  · Ejemplo (Práctica 6, Ej 3, grado 3): $x=(-1,0,2,3)$, $y=(-1,3,11,27)$ → $a=(3,2,-1,1)$ → $p(x)=3+2x-x^2+x^3$. Verificación: $p(-1)=-1$, $p(0)=3$, $p(2)=11$, $p(3)=27$ ✓ (pasa exacto por los 4 puntos, $V$ cuadrada).
- **Alternativa — Lagrange** (sin resolver sistema):
  $$p(x)=\sum_{i=0}^n y_i\,L_i(x), \qquad L_i(x)=\prod_{j\neq i}\frac{x-x_j}{x_i-x_j}$$
  · Por qué funciona: $L_i(x_k)=0$ si $k\neq i$ (numerador tiene factor $(x_k-x_k)$), $L_i(x_i)=1$ → $p(x_k)=y_k$.
  · Mismos denominadores $x_i-x_j$ que en $\det(V)$ — prueba constructiva alternativa de existencia/unicidad.
  · Trade-off: evita resolver el sistema, pero no da los coeficientes $a_i$ en base canónica — para eso, Vandermonde/Gauss.
- **Error de interpolación:** si $f\in C^{n+1}$, existe $\xi$ tal que
  $$f(x)-p_n(x)=\frac{f^{(n+1)}(\xi)}{(n+1)!}\prod_{i=0}^n(x-x_i)$$
  (análogo al resto de Taylor: error controlado por la derivada $(n+1)$-ésima y la distancia de $x$ a los nodos.)
- ** Fenómeno de Runge** (Práctica 6, Ej 7): interpolar $f(x)=\frac1{1+25x^2}$ en $n+1$ nodos **equiespaciados** en $[-1,1]$ — aumentar $n$ **empeora** el error cerca de los bordes.
  · Por qué: en $\prod_i(x-x_i)$, cerca del centro los factores se compensan (signos mezclados, magnitud chica); cerca de un borde casi todos tienen el mismo signo y magnitud grande (no hay nodos "más allá") → el producto crece fuerte en los extremos, y empeora con $n$ (además $f^{(n+1)}$ de Runge también crece).
  · Consecuencia: para $n$ alto, oscila violentamente cerca de $\pm1$ aunque ajuste bien en el centro — es **intrínseco** al esquema (nodos equiespaciados + grado alto), no un problema de cómo se resuelve el sistema.
  · Salida: interpolar **a trozos** con polinomios de grado bajo (motiva trapecios compuesto, 5.12), o usar nodos de Chebyshev (concentrados en los bordes).
  · $V$ suele estar muy mal condicionada para $n$ grande — síntoma relacionado (distinto) del mismo problema de fondo.

### 5.3 Evaluación de polinomios — Horner (Práctica 6, Ej 1)

$p(x)=a_nx^n+\dots+a_1x+a_0$ — el costo de evaluar depende de cómo se escriba (relevante: el interpolador de 5.2 se evalúa muchas veces).

| Forma | Costo |
| :--- | :--- |
| Directa (cada $x_0^i$ recalculada) | $O(n^2)$ productos |
| Horner: $p(x)=a_0+x(a_1+x(a_2+\dots+x(a_{n-1}+xa_n)\dots))$ | $O(n)$: 1 producto + 1 suma por nivel |

Por eso NumPy (`np.polyval`, `np.poly1d`) evalúa así en vez de término a término.

### 5.4 Linealización de modelos no lineales (Práctica 6, Ej 8, 9, 10)

Cuando el modelo es **no lineal en los parámetros**, un cambio de variable a veces lo convierte en cuadrados mínimos lineal estándar (5.1).

- **Ej 10 — ley de potencia $Y=aX^b$:** $\ln Y=\ln a+b\ln X$ — lineal en $(\ln a,b)$. $Y'=\ln Y$, $X'=\ln X$, ajustar $Y'=c+bX'$ (recta), recuperar $a=e^c$.
- **Ej 8 — $g(x)=e^{p(x)}$, $p$ polinomio:** $\ln g(x)=p(x)$ — lineal en los coeficientes; ajustar/interpolar $p$ sobre $\ln(y_i)$ con Vandermonde (5.2). Ej: 6 puntos (años 1950–2000), $p$ grado 5 → $V$ **cuadrada** → interpolación exacta de $\ln(y_i)$, no ajuste aproximado.
- **Ej 9 — $f(t)=at^2+b$ (sin término lineal):** columnas de diseño $[t^2,1]$, no potencias consecutivas → `np.polyfit(t, altura, 2)` no sirve (arma Vandermonde completa con $t^1$ de más). Truco: $u=t^2$, ajustar $f(u)=au+b$ como recta (`np.polyfit(u, altura, 1)`).
-  Linealizar cambia **en qué escala** se minimiza el error: ajustar $\ln Y$ vs $\ln X$ minimiza $\sum(\ln y_i-(\ln a+b\ln x_i))^2$ (error relativo/multiplicativo), **no** $\sum(y_i-ax_i^b)^2$ (error absoluto) — son problemas de optimización distintos, aproximación práctica estándar pero no equivalente.

### 5.5 Cuadrados mínimos con funciones base generales (Práctica 6, Ej 11, 12, 13)

Generalización de 5.1/5.2: columnas de $A$ = evaluaciones de **cualquier familia** $f_1,\dots,f_m$ en $x_1,\dots,x_n$:
$$A=\begin{pmatrix}f_1(x_1)&\cdots&f_m(x_1)\\\vdots&&\vdots\\f_1(x_n)&\cdots&f_m(x_n)\end{pmatrix}\in\mathbb R^{n\times m}$$
Minimizar $\|A\alpha-y\|_2$ con las mismas ecuaciones normales $A^TA\alpha=A^Ty$. $f(x)=\sum_i\alpha_if_i(x)$ — el polinomio es el caso $f_i(x)=x^{i-1}$.

- **Ej 11:** implementación genérica dado $\{f_1,\dots,f_m\}$.
- **Ej 12 (humedad, periodicidad anual):** $f_1(t)=1$, $f_2(t)=\sin(t\tfrac{2\pi}{366})$, $f_3(t)=\sin(t\tfrac{4\pi}{366})$.
- **Ej 13 (regresión lineal múltiple):** $Y=\beta_0+\beta_1X_1$ o $+\beta_2X_2$, $A=[\mathbf1|X_1|X_2]$.
  · $\beta_i$ = cambio promedio en $Y$ por unidad de $X_i$, **manteniendo constantes las demás variables** ($\beta_1$ múltiple ≠ $\beta_1$ simple en general). Ítem (d): efecto de 10g extra ($\beta_2$=coef. de peso) es $10\beta_2$ cm, a edad gestacional constante.

### 5.6 Métodos iterativos — marco general
$A=D+L+U$ (diagonal, triangular inf./sup. estrictas). Forma general: $x^{(n+1)}=Tx^{(n)}+c$.

**Teorema clave:** converge para cualquier $x^{(0)}\iff\rho(T)<1$, $\rho(T)=\max\{|\lambda|:\lambda$ autovalor de $T\}$.
- Splitting $A=M-N$, $T=-M^{-1}N$: $\lambda$ autovalor de $T\iff\det(\lambda M+N)=0$.
- Converge en $n$ pasos exactos $\iff T^{n+1}=0$.
- ** $\lim_{n\to\infty}A^n=0\iff\rho(A)<1$** (para **cualquier** $A$ cuadrada, no hace falta que sea diagonalizable — sale de la forma de Jordan: cada bloque $J_\lambda^n\to0\iff|\lambda|<1$, porque $J_\lambda^n$ tiene entradas $\binom nk\lambda^{n-k}$ que van a 0 si $|\lambda|<1$ pese al crecimiento polinomial de $\binom nk$). Es la base del "Teorema clave" de arriba, y la herramienta directa para ejercicios tipo "probar que $\lim B^n=0$" (Práctica 7, Ej 8) — **no confundir con nilpotencia** ($B^k=0$ para algún $k$ finito): son cosas distintas, $\rho(B)<1$ no implica que $B$ sea nilpotente.
- $\rho(B)=\lim_{n\to\infty}\|B^n\|^{1/n}$ para toda norma subordinada — **atención**: esto es un límite (Gelfand), no vale $\rho(B)=\|B^n\|^{1/n}$ para un $n$ finito particular.
- $\rho(T)\leq\|T\|$ para cualquier norma inducida —  por eso puede haber $\|T\|\geq1$ y converger igual (importa $\rho$, no una norma particular). Buscar la norma con $<1$ es ejercicio típico.
-  Mal condicionamiento de $A$ ≠ convergencia lenta — no confundir.

### 5.7 Jacobi
$$x^{(n+1)}=-D^{-1}(L+U)x^{(n)}+D^{-1}b$$
Paralelizable (cada componente depende solo de $x^{(n)}$).
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

### 5.8 Gauss-Seidel
$$x^{(n+1)}=-(D+L)^{-1}Ux^{(n)}+(D+L)^{-1}b$$
No paralelizable (usa valores ya actualizados en la misma iteración).
- Si $A$ tridiagonal: $\rho(B_{GS})=\rho(B_J)^2$ (GS converge más rápido).

### 5.9 SOR (Successive Over-Relaxation)
Generaliza Gauss-Seidel con parámetro $\omega$: matriz de iteración $B(\omega)$.
- $\det(B(\omega))=(1-\omega)^n$ ⇒ condición **necesaria**: $\omega\in(0,2)$.
- $\omega=1$ recupera Gauss-Seidel.

### 5.10 Método del gradiente
$A$ simétrica DP: resolver $Ax=b$ ⟺ minimizar $f(x)=\tfrac12x^TAx-b^Tx$ (mínimo único si $A$ DP).
- Dirección de máximo descenso: $-\nabla f(x)=b-Ax=r^{(k)}$ (residuo).
- Paso óptimo: $t^*=\dfrac{(r^{(k)})^Tr^{(k)}}{(r^{(k)})^TAr^{(k)}}$.
- Iteración: $x^{(k+1)}=x^{(k)}+t^*r^{(k)}$.
- Converge en zigzag sobre curvas de nivel (a diferencia del gradiente conjugado, que lo corrige — no cubierto en la práctica).

### 5.11 Radio espectral — resumen
$\rho(B)<1$: converge · $\rho(B)=1$: puede no converger o muy lento · $\rho(B)>1$: diverge. Cuanto menor $\rho(B)$, más rápida la convergencia.
- **$\lim_{n\to\infty}B^n=0\iff\rho(B)<1$** — la razón de fondo de todo lo anterior (ver 5.6). Vale para cualquier $B$, diagonalizable o no.

### 5.12 Integración numérica — trapecios compuesta (Práctica 6, Ej 14)

**Motivación:** interpolar con grado alto puede ser malo (Runge, 5.2) — alternativa: interpolar **a trozos** con polinomios de bajo grado (lineales), y usar esos tramos para aproximar la integral.

$f:[a,b]\to\mathbb R$, partición uniforme $x_i=a+ih$, $h=\frac{b-a}n$.

1. Recta $p_i$ que interpola $f$ en $x_i,x_{i+1}$: $p_i(x)=f(x_i)+\frac{f(x_{i+1})-f(x_i)}h(x-x_i)$.
2. Integral exacta (trapecio de bases $f(x_i)$, $f(x_{i+1})$, altura $h$): $\int_{x_i}^{x_{i+1}}p_i\,dx=h\cdot\frac{f(x_i)+f(x_{i+1})}2$.
3. Sumar: cada nodo interior aparece en 2 trapecios (peso $h$ total), extremos $a,b$ una vez (peso $h/2$):
   $$I\sim\frac{b-a}n\left[\frac{f(a)+f(b)}2+\sum_{i=1}^{n-1}f(x_i)\right]$$

**Implementación (vectorizada, sin `for`):** evaluar $f$ en `x=np.linspace(a,b,n+1)`, sumar con pesos `[0.5,1,...,1,0.5]*h`, o `h*(f(x).sum() - (f(a)+f(b))/2)`.

-  Esta es la versión *compuesta* — no confundir con trapecio *simple* (un solo trapecio en $[a,b]$, mucho menos preciso).
- **Error:** $O(h^2)$ global (cada trapecio $O(h^3)$, hay $O(1/h)$ trapecios) → duplicar $n$ reduce el error a un cuarto.

---

## 6. Apéndice: errores comunes / V-F rápidas

- **Entradas positivas ≠ definida positiva.** Chequear autovalores o menores principales, no el signo de las entradas.
- **$\det(A)\to0$ no implica mal condicionamiento** (ej. $\varepsilon I$).
- **Mal condicionamiento ≠ convergencia lenta de iterativos** — propiedades distintas de $A$.
- **$\rho(B)<1$ (i.e. $\lim B^n=0$) ≠ $B$ nilpotente** ($B^k=0$ para algún $k$ finito) — son propiedades distintas; $\rho(B)<1$ solo dice que la sucesión de potencias tiende a 0, no que se anule de golpe.
- **$\|T\|\geq1$ no implica que el iterativo diverja** — decide $\rho(T)<1$, no una norma particular.
- **$Ax\neq Qx$ en general** aunque $\operatorname{Im}(A)=\operatorname{Im}(Q)$ en QR — TL distintas.
- **Un ciclo en Markov genera autovalores de módulo 1 ≠ 1** — equilibrio único no implica convergencia desde todo $v^{(0)}$.
- **No toda matriz de Markov es diagonalizable** — chequear $mg=ma$.
- **Householder sin la optimización de rango 1 es $O(n^4)$, no $O(n^3)$** — la ganancia viene de no armar $H_k$ explícitamente.
- **LU no existe para toda matriz** (Schur y SVD sí); cero en la diagonal durante la eliminación no siempre mata la factorización, pero cero en $a_{11}$ sí.
- **Cholesky y LU son las únicas descomposiciones "no universales"** de la tabla del punto 0 — QR, Schur, SVD existen siempre bajo condiciones mínimas de forma.

### 6.1 Numpy — funciones útiles

**`np.tril(A, k=0)` / `np.triu(A, k=0)`:** parte triangular inferior/superior (en o por debajo/encima de la diagonal $k$), resto en cero. `k=0`=diagonal principal, `k>0` sube, `k<0` baja.
```python
A = np.array([[1, 2, 3],
              [4, 5, 6],
              [7, 8, 9]])
np.tril(A)   # [[1,0,0],[4,5,0],[7,8,9]]
np.triu(A)   # [[1,2,3],[0,5,6],[0,0,9]]
```

**`np.diag`:** doble uso según dimensión del input.
- Input 2D (matriz) → **extrae** la diagonal: `np.diag(A)` → `[1,5,9]`.
- Input 1D (vector) → **construye** matriz diagonal: `np.diag([1,5,9])`.

**Conjugar:**
- Elemento a elemento: `np.conj(A)` (o `A.conj()`/`A.conjugate()`).
- Conjugada transpuesta / adjunta hermítica ($A^H=(\bar A)^T$, ver 3.3, 4.1-4.2): `A.conj().T`.
- Con matrices reales, `np.conj(A)` no cambia nada — solo importa con `dtype=complex`.
