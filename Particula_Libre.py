import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from math import e, pi 

"Graficamos la funcion phi(k) para a = 1/(2*pi)"

k = np.arange(-10, 10, 0.1)
phik = e**((-pi*k**2)/2)

plt.plot(k, phik)
plt.title("Φ(k)")
plt.show()

"Graficamos la funcion psi(x, t) para a = 1/(2*pi) y t=0 y con m = masa del electron"
"psixt = (2*a/pi)**(1/4)*(1/(1+2j*h*a*t)/m)**(1/2)*e**((-a*x**2)/(1+(2j*h*a*t)/m))"

a = 1/(2*pi)
t = 0
m = 9.11*10**(-31)
h = 1.054571817*10**(-34) 
x = np.arange(-10, 10, 0.1)

psixt = (2*a/pi)**(1/4)*e**(-a*x**2)

plt.plot(x, psixt)
plt.title("Ψ(x,t = 0)")
plt.show()

