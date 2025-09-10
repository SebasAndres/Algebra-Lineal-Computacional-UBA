import pytest 
import numpy as np

from funciones import *

SAMPLE_VECTOR = [1,2,3]
SAMPLE_MATRIZ_2_2 = np.array([[1, 2], [3, 4]])

def test_norma_2():
    assert norma(
        x=SAMPLE_VECTOR,
        p=2
    ) == np.linalg.norm(SAMPLE_VECTOR, 2)


def test_norma_5():
    assert norma(
        x=SAMPLE_VECTOR,
        p=5
    ) == np.linalg.norm(SAMPLE_VECTOR, 5)


def test_normalizado():
    """Valida si la funcion normaliza cumple su trabajo"""
    assert np.allclose(normaliza(SAMPLE_VECTOR, 2), SAMPLE_VECTOR / np.linalg.norm(SAMPLE_VECTOR, 2))


def test_norma_inf():
    """Valida si la funcion norma toma la norma infinito"""
    assert np.linalg.norm(SAMPLE_VECTOR, np.inf) == \
        norma(SAMPLE_VECTOR, 'inf')


def test_norma_uno():
    """
    Valida que funcione la norma 1
    """
    assert norma1Exacta(
        SAMPLE_MATRIZ_2_2
    ) == 6

def test_norma_inf():
    """
    Valida que funcione la norma 1
    """
    assert normaInfExacta(
        SAMPLE_MATRIZ_2_2
    ) == 7

# def test_norma_mat_mc():
#     """
#     Valida si la funcion normaMatMc estima la norma matricial con
#     met. montecarlo
#     """
#    matriz = np.array([[1, 2], [3, 4]])
#    assert np.allclose(
#        normaMatMC(matriz, p=2, q=2, Np=1000), 
#        np.linalg.norm(matriz, ord=2)
#    )

# def test_nro_cond():
#     """
#     Valida que el nro de condicion se estime bien
#     """
#     condicionEstimada = condMc(SAMPLE_MATRIZ_2_2, p=2)
#     condicionNumpy  = np.linalg.cond(SAMPLE_MATRIZ_2_2, p=2)
#     assert np.allclose(
#         condicionEstimada,
#         condicionNumpy
#     )