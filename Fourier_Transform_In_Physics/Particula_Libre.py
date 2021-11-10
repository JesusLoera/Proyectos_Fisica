import numpy as np
import matplotlib.pyplot as plt
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

"Graficamos la funcion rho(x, t) para a = 1/(2*pi) y t=0"
"cuando t -> 0, θ -> 0"

w = a**(1/2)
x = np.arange(-10, 10, 0.1)

rhoxt = ((2/pi)**(1/2))*w*e**(-2*w**2*x**2)

plt.plot(x, rhoxt)
plt.title("ρ(x,t = 0)")
plt.show()


"Graficamos la funcion rho(x, t) para a = 1/(2*pi) y t=10^5 s y con m = masa del electron"

t = 10**5
x = np.arange(-10, 10, 0.001)

rhoxt = ((2/pi*a)**(1/2))*(m/(h*t))*e**(-(x**2*m**2/(2*h**2*a*t**2)))


plt.plot(x, rhoxt)
plt.title("ρ(x,t = 10^5)")
plt.show()


"Graficamos la funcion rho(x, t) t=10^7 s"

t = 10**7
x = np.arange(-10, 10, 0.001)

plt.plot(x, rhoxt)
plt.title("ρ(x,t = 10^7)")
plt.show()


"Graficamos la funcion rho(x, t) t=10^9 s"

t = 10**9
x = np.arange(-10, 10, 0.001)

plt.plot(x, rhoxt)
plt.title("ρ(x,t = 10^9 s)")
plt.show()