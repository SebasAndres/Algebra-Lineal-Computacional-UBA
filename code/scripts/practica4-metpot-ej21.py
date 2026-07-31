import numpy as np


def metPot(A: np.ndarray, tol: float = 1e-10, max_iter: int = 10000) -> tuple[float, np.ndarray]:
    n = A.shape[0]
    v = np.random.rand(n)
    v = v / np.linalg.norm(v)
    r_prev = 0.0
    for _ in range(max_iter):
        w = A @ v
        v = w / np.linalg.norm(w)
        r = v @ A @ v
        if abs(r - r_prev) < tol:
            break
        r_prev = r
    return r, v


def metPotInversa(A: np.ndarray, tol: float = 1e-10, max_iter: int = 10000) -> tuple[float, np.ndarray]:
    n = A.shape[0]
    v = np.random.rand(n)
    v = v / np.linalg.norm(v)
    r_prev = 0.0
    for _ in range(max_iter):
        w = np.linalg.solve(A, v)
        v = w / np.linalg.norm(w)
        r = v @ A @ v
        if abs(r - r_prev) < tol:
            break
        r_prev = r
    return r, v
