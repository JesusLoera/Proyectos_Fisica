import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from Methods import FFT

"""
El atractor de Lorentz
"""


def Lorentz_ED(xyz, t, sigma, rho, beta):
    x, y, z = xyz
    return [sigma * (y - x), x * (rho - z) - y, x * y - beta * z]


"""
Discretizamos la variable dependiente en N=2**m terminos,
el cual es un requerimento de la Transformada Rapida de Fourier
"""
m = 9
N = 2 ** m


"""
Condiciones iniciales y parametros del Atractor de Lorentz
"""

# Asignamos valores a los parámetros
sigma, Rho, beta = 10, np.linspace(24, 25, 100), 8.0 / 3.0

# Condición inicial y valores de t sobre los que calcular
xyz0 = [1.0, 1.0, 1.0]
xi = 0.0
xf = 25.0
t = np.linspace(xi, xf, N)
iter = 1

for rho in Rho:

    # Resolvemos las ecuaciones
    xyz = odeint(Lorentz_ED, xyz0, t, args=(sigma, rho, beta))

    # Graficamos las soluciones
    from mpl_toolkits.mplot3d.axes3d import Axes3D

    fig, (ax1) = plt.subplots(1, 1, subplot_kw={"projection": "3d"})
    ax1.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], "r", alpha=0.5)
    ax1.set_xlabel("$x$", fontsize=16)
    ax1.set_ylabel("$y$", fontsize=16)
    ax1.set_zlabel("$z$", fontsize=16)
    plt.title("Diagrama de fase")
    plt.savefig("images_numpy_24_25/fase/fase" + str(iter) + ".png")
    plt.close()

    """
    Calculamos la FFT de X(t) de las ecns de Lorentz
    """
    G_Lorentz_X = np.fft.fft(xyz[:, 0])

    """
    Creamos un arreglo con los espectros de potencia de la FFT
    """
    S_fft = []
    for i in range(N):
        s = (1.0 / N) * abs(G_Lorentz_X)
        S_fft = np.append(S_fft, s)

    """
    Graficamos la función Fx_fft
    """
    plt.plot(t, xyz[:, 0])
    plt.title("X(t)")
    plt.savefig("images_numpy_24_25/eje_x/eje_x" + str(iter) + ".png")
    plt.close()

    """"
    Ahora Graficamos S vs frecuencia (Espectro)
    """
    Freq = np.arange(0, N // 2)
    plt.plot(Freq, S_fft[: N // 2])
    plt.title("Espectro de Frecuencias")
    plt.savefig("images_numpy_24_25/espectro/espectro" + str(iter) + ".png")
    plt.close()

    iter = iter + 1