import numpy as np
import matplotlib.pyplot as plt


def F(X):
    Re_F = []
    Im_F = []
    for x in X:
        re_f = np.sin(50.0 * 2.0 * np.pi * x) + 0.5 * np.sin(80.0 * 2.0 * np.pi * x)
        Re_F = np.append(Re_F, re_f)
        im_f = 0
        Im_F = np.append(Im_F, im_f)
    return [Re_F, Im_F]


def DFT(F, N):
    Im_F = F[0]
    Re_F = F[1]
    Im_G = []
    Re_G = []
    G = []

    for j in range(len(Re_F)):
        im_g = 0
        re_g = 0
        for i in range(len(Re_F)):
            re_g = (
                re_g
                + np.cos((2 * np.pi * j * i) / (N)) * Re_F[i]
                + np.sin((2 * np.pi * j * i) / (N)) * Im_F[i]
            )
            im_g = (
                im_g
                + np.cos((2 * np.pi * j * i) / (N)) * Im_F[i]
                - np.sin((2 * np.pi * j * i) / (N)) * Re_F[i]
            )
        # re_g = re_g / (N ** (0.5))
        # im_g = im_g / (N ** (0.5))
        Im_G = np.append(Im_G, im_g)
        Re_G = np.append(Re_G, re_g)
        G = np.append(G, [re_g, im_g])

    return G, Re_G, Im_G


"""
Discretizamos la variable dependiente en N=2**m terminos,
el cual es unre querimento de la Transformada Rapida de Fourier
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
G, Re_G, Im_G = DFT(Fx, N)

"""
Creamos un arreglo para la frecuencia
"""
Freq = np.linspace(xi, 1.0 / (2 * xf), N)


"""
Creamos un arreglo con los espectros de potencia de la DFT
"""
S = []
for i in range(N):
    s = (2.0 / N) * ((Re_G[i]) ** 2 + (Im_G[i]) ** 2) ** (0.5)
    S = np.append(S, s)

print(len(S))
print(len(Freq))

""""
Ahora Graficamos S vs f
"""
plt.plot(Freq, S)
plt.title("Espectro de Frecuencias")
plt.show()

"""
Ahora comparamos con la fft de numpy
"""
Sp = np.fft.fft(Fx[0])
Freq_py = np.fft.fftfreq(X.shape[-1])

P = []
for i in range(len(Sp)):
    p = (2.0 / N) * abs(Sp[i])
    P = np.append(P, p)

plt.plot(Freq[: (N / 2)], P[: (N / 2)])
plt.title("Espectro de Frecuencias")
plt.show()
