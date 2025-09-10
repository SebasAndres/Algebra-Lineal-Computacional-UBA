import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def proyectarPts(T, wz):
    assert(T.shape == (3,3)) # chequeo de matriz 3x3
    assert(T.shape[1] == wz.shape[0]+1) # multiplicacion matricial valida   
    wz = np.vstack((wz, np.ones((1, wz.shape[1]))))
    xy = T @ wz     
    return xy[:-1, :]

def pointsGrid(corners):
    # crear 10 lineas horizontales
    [w1, z1] = np.meshgrid(np.linspace(corners[0,0], corners[1,0], 46),
                        np.linspace(corners[0,1], corners[1,1], 10))

    [w2, z2] = np.meshgrid(np.linspace(corners[0,0], corners[1,0], 10),
                        np.linspace(corners[0,1], corners[1,1], 46))

    w = np.concatenate((w1.reshape(1,-1),w2.reshape(1,-1)),1)
    z = np.concatenate((z1.reshape(1,-1),z2.reshape(1,-1)),1)
    wz = np.concatenate((w,z))
                         
    return wz
          
def vistform(T, wz, titulo=''):
    # transformar los puntos de entrada usando T
    xy = proyectarPts(T, wz)
    if xy is None:
        print('No fue implementada correctamente la proyeccion de coordenadas')
        return

    # calcular los limites para ambos plots
    minlim = np.min(np.concatenate((wz, xy), 1), axis=1)
    maxlim = np.max(np.concatenate((wz, xy), 1), axis=1)

    bump = [np.max(((maxlim[0] - minlim[0]) * 0.05, 0.1)),
            np.max(((maxlim[1] - minlim[1]) * 0.05, 0.1))]
    limits = [[minlim[0]-bump[0], maxlim[0]+bump[0]],
               [minlim[1]-bump[1], maxlim[1]+bump[1]]]             

    fig, (ax1, ax2) = plt.subplots(1, 2)         
    fig.suptitle(titulo)
    grid_plot(ax1, wz, limits, 'w', 'z')    
    grid_plot(ax2, xy, limits, 'x', 'y')    
    
def grid_plot(ax, ab, limits, a_label, b_label):
    ax.plot(ab[0,:], ab[1,:], '.')
    ax.set(aspect='equal',
           xlim=limits[0], ylim=limits[1],
           xlabel=a_label, ylabel=b_label)
    plt.show()


def main():
    print('Ejecutar el programa')
    # generar el tipo de transformacion dando valores a la matriz T
    T = np.array([
        [1, 0.4, 0],
        [0, 1, 0],
        [0, 0, 1]
    ])
    corners = np.array([[0,0],[100,100]])
    # corners = np.array([[-100,-100],[100,100]]) # array con valores positivos y negativos
    wz = pointsGrid(corners)
    vistform(T, wz, 'Deformar coordenadas')
    
    
if __name__ == "__main__":
    main()
