

Ej 7/8/25

1.a. definir proyector p / im(f) = im(p) ; nu(p) = nu(f)

Visualmente como f(x)=Ax, se ve que las primeras dos columnas son LD. 

Entonces:
- c4 = c2 + c3 
- c2 = - c1

Luego :

im(f) = <(-1,1,1,1), (0,-1,0,1)>

Para el nucleo se puede usar:

c2 + c1 = 0 --> A(1,1,0,0) = 0 --> (1,1,0,0) in Nu(f)

c4 - c2 - c3 = 0 --> A(0,-1,-1,1) = 0 --> (0,-1,-1,1) in Nu(f)

Nu(f) = <(1,1,0,0), (0,-1,-1,1)>


:: definir un proyector es definir una TL que cumple p^2 = p <-> p(v) = v si v \in Im(p).
Es un si solo si, entonces si p(v)=v , v\in Im(p) --> p^2 = p

Para definir una TL tengo que definirla sobre la base.

Si quiero un proyector ortogonal, ortogonalizo con GS sobre la base para la cual defino el proyector


Si quiero armar un proyector ortogonal sobre <u1,u2> --> (u_1 | u_2) ( u_1^T
								       u_2^T )

Proyectores ortogonales:
P<u1>(x) = (u1^t x/u1^t u1) u1 = u1^x u1 = u1 u1^t x 

Las matrices u_i u_i^t siempre tienen rango 1


A = uv^t tiene imagen <u> y el nucleo son los ortogonales a v (nu(A) = <v>^\perp, im(a)=<u>).


El tip en la matríz diagonal es que f(x) = lambda * x en proyecciones es con lambda \in {0, 1}

--> si tomo B como Im(p) ++ Nu(p) (para este caso en particular) 


Si B = {v1, v2, v3, v4} --> [v1]_B = e_1 

[f]_{BB}(v_B) = (f(v))_B

---


Recordar $Ax = x1 col_1(A) + ... + x_n col_n(A)$ para $x=(x_1,...x_n)$

im([0,1],[0,0]) = < (0,0), (1,0 )> = <(1,0)>

nu(B) = (1,0) 

El nucleo y la imagen pueden ser el mismo

Be_1 = col_1(B)

