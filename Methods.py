import numpy as np
import matplotlib.pyplot as plt

"""
DFT, FFT, IDFT, IFFT
"""


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