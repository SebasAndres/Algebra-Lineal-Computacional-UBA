import numpy as np
import numpy.typing as npt


def cumpleRequisitosMetodoPotencia(
    eigvals: npt.NDArray[np.complexfloating],
    eigvecs: npt.NDArray[np.complexfloating],
) -> bool:
    """
    Tengo que validar que los autovectores forman una base (A diagonalizable)
    y que hay un unico autovalor de modulo maximo: |l1| > |l2| >= ... >= |lN|.
    """
    if np.linalg.matrix_rank(eigvecs) < eigvecs.shape[1]:
        return False

    moduli = np.sort(np.abs(eigvals))[::-1]
    return moduli[0] > moduli[1]


def estimarAutovaloryAutovectorDominante(A: npt.NDArray[np.floating], epochs: int = 10_000) -> float:
    """ Método de la potencia usando el cociente de Reileigh"""
    N, _ = A.shape
    x_k = np.random.random(N)
    x_k /= np.linalg.norm(x_k)

    rk = 0.0
    for _ in range(epochs):
        x_k = A @ x_k
        x_k /= np.linalg.norm(x_k)
    
    rk = x_k.T @ A @ x_k / (x_k.T @ x_k)
    return rk, x_k



def parteA(nombre: str, M: npt.NDArray[np.floating]) -> tuple[npt.NDArray, npt.NDArray]:
    eigvals, eigvecs = np.linalg.eig(M)
    cumple = cumpleRequisitosMetodoPotencia(eigvals, eigvecs)
    print(f"{nombre}: autovalores = {eigvals}")
    print(f"{nombre}: cumple hipotesis del metodo de la potencia = {cumple}")
    return eigvals, eigvecs


def parteB(eigvals: npt.NDArray, eigvecs: npt.NDArray) -> None:
    """ S = span de los autovectores de A que NO son el dominante """
    dominante = int(np.argmax(np.abs(eigvals)))
    print("S = span{")
    for i in range(len(eigvals)):
        if i != dominante:
            print(f"  v (lambda={eigvals[i]:.3f}) = {eigvecs[:, i]}")
    print("}")


def parteC(eigvals: npt.NDArray, eigvecs: npt.NDArray) -> float:
    """
    Busca alpha tal que v0=(-1,alpha,-2) tenga coeficiente 0 en la
    direccion del autovector dominante, para que el metodo converja
    al segundo autovalor de mayor modulo.
    """
    dominante = int(np.argmax(np.abs(eigvals)))

    base = np.array([-1.0, 0.0, -2.0])
    direccion = np.array([0.0, 1.0, 0.0])

    coefBase = np.linalg.solve(eigvecs, base)
    coefDireccion = np.linalg.solve(eigvecs, direccion)

    return -coefBase[dominante] / coefDireccion[dominante]


def main():
    A = np.array([
        [-6.0, 9.0, 3.0],
        [0.0, 8.0, -2.0],
        [0.0, -1.0, 7.0],
    ])
    B = np.array([
        [5.0, 9.0, 6.0],
        [-3.0, -7.0, -6.0],
        [0.0, 0.0, -1.0],
    ])

    eigvalsA, eigvecsA = parteA("A", A)
    parteA("B", B)

    parteB(eigvalsA, eigvecsA)

    eigvalsB, eigvecsB = np.linalg.eig(B)
    alpha = parteC(eigvalsB, eigvecsB)
    print(f"alpha = {alpha}")

    l1, v1 = estimarAutovaloryAutovectorDominante(A)
    print(f"lambda1 estimado (metodo potencia) = {l1}")
    print(f"v1 estimado (metodo potencia) = {v1}")


if __name__ == '__main__':
    main()
