import numpy as np
import matplotlib.pyplot as plt

"""
DFT, FFT, IDFT, IFFT
"""


def F(X):
    Re_F = []
    Im_F = []
    for x in X:
        re_f = np.sin(50.0 * 2.0 * np.pi * x) + 0.5 * np.sin(80.0 * 2.0 * np.pi * x)
        Re_F = np.append(Re_F, re_f)
        im_f = 0
        Im_F = np.append(Im_F, im_f)
    return [Re_F, Im_F]


def F_fft(X):
    Re_F = []
    for x in X:
        re_f = np.sin(50.0 * 2.0 * np.pi * x) + 0.5 * np.sin(80.0 * 2.0 * np.pi * x)
        Re_F = np.append(Re_F, re_f)
    return Re_F


def DFT(F, N):
    Im_F = F[0]
    Re_F = F[1]
    Im_G = []
    Re_G = []
    G = []
    Freq = []

    for k in range(N):
        im_g = 0
        re_g = 0
        for n in range(N):
            re_g = (
                re_g
                + np.cos((2 * np.pi * k * n) / (N)) * Re_F[n]
                + np.sin((2 * np.pi * k * n) / (N)) * Im_F[n]
            )
            im_g = (
                im_g
                + np.cos((2 * np.pi * k * n) / (N)) * Im_F[n]
                - np.sin((2 * np.pi * k * n) / (N)) * Re_F[n]
            )
        Im_G = np.append(Im_G, im_g)
        Re_G = np.append(Re_G, re_g)
        G = np.append(G, [re_g, im_g])
        Freq = np.append(Freq, k)

    return G, Re_G, Im_G, Freq


"""
def I_DFT(F, N):
    Im_F = F[0]
    Re_F = F[1]
    Im_G = []
    Re_G = []
    G = []

    for k in range(N):
        im_g = 0
        re_g = 0
        for n in range(N):
            re_g = (
                re_g
                + np.cos((2 * np.pi * k * n) / (N)) * Re_F[n]
                + np.sin((2 * np.pi * k * n) / (N)) * Im_F[n]
            )
            im_g = (
                im_g
                + np.cos((2 * np.pi * k * n) / (N)) * Im_F[n]
                + np.sin((2 * np.pi * k * n) / (N)) * Re_F[n]
            )
        Im_G = np.append(Im_G, im_g)
        Re_G = np.append(Re_G, re_g)
        G = np.append(G, [re_g, im_g])

    Im_G = (1.0 / N) * Im_G
    Re_G = (1.0 / N) * Re_G
    G = (1.0 / N) * G

    return G, Re_G, Im_G
"""


def FFT(f):
    N = len(f)
    if N <= 1:
        return f

    # division
    even = FFT(f[0::2])
    odd = FFT(f[1::2])

    # store combination of results
    G = np.zeros(N).astype(np.complex64)

    # only required to compute for half the frequencies
    # since u+N/2 can be obtained from the symmetry property
    for u in range(N // 2):
        G[u] = even[u] + np.exp(-2j * np.pi * u / N) * odd[u]  # conquer
        G[u + N // 2] = even[u] - np.exp(-2j * np.pi * u / N) * odd[u]  # conquer

    return G


"""
Discretizamos la variable dependiente en N=2**m terminos,
el cual es un requerimento de la Transformada Rapida de Fourier
"""
m = 10
N = 2 ** m

"""
Ahora creamos un arreglo con los 2**m valores discretos en [xi,xf]
"""
xi = 0.0
xf = 1.0
X = np.linspace(xi, xf, N)

Fx = F(X)
Fx_fft = F_fft(X)
G, Re_G, Im_G, Freq = DFT(Fx, N)

"""
Creamos un arreglo con los espectros de potencia de la DFT
"""
S = []
for i in range(N):
    s = (1.0 / N) * ((Re_G[i]) ** 2 + (Im_G[i]) ** 2) ** (0.5)
    S = np.append(S, s)

"""
Graficamos la función Re_f(x)
"""
plt.plot(X, Fx[0])
plt.title("Re_Fx")
plt.show()

"""
Graficamos la función Im_f(x)
"""
plt.plot(X, Fx[1])
plt.title("Im_Fx")
plt.show()

""""
Ahora Graficamos S vs frecuencia (Espectro)
"""
plt.plot(Freq[: N // 2], S[: N // 2])
plt.title("Espectro de Frecuencias")
plt.show()

"""
Ahora probamos FFT
"""
G_fft = FFT(Fx_fft)

"""
Creamos un arreglo con los espectros de potencia de la FFT
"""
S_fft = []
for i in range(N):
    s = (1.0 / N) * abs(G_fft)
    S_fft = np.append(S_fft, s)

"""
Graficamos la función Fx_fft
"""
plt.plot(X, Fx_fft)
plt.title("Fx_fft")
plt.show()

""""
Ahora Graficamos S vs frecuencia (Espectro)
"""
plt.plot(Freq[: N // 2], S_fft[: N // 2])
plt.title("Espectro de Frecuencias")
plt.show()

"""
Ahora probamos la transformada Inversa de Fourier

Gk = [Re_G, Im_G]
F_idft, Re_F_idft, Im_F_idft = I_DFT(Fx, N)

Ahora Graficamos F vs F

plt.plot(X, Re_F_idft)
plt.title("Fx")
plt.show()

"""
