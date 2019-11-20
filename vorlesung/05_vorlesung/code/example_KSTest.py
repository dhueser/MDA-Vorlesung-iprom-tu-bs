import numpy as np
import matplotlib.pyplot as plt
import scipy.special
def cumulatives(x):
    xbar = np.mean(x)
    s = np.std(x, ddof=1)
    J = len(x)
    xsort = np.sort(x)
    H = np.arange(1,J+1)/J;
    z = (xsort - xbar)/s
    # The cumulative of the unit normal distribution is given 
    # by Phi(z) = 1/2[1 + erf(z/sqrt(2))]
    P = 0.5 * (1 + scipy.special.erf(z/np.sqrt(2)))
    return P, H, xsort, xbar, s
#
# Hauptprogramm
#
# generiere synthetische Stichprobe
stichprobenumfang = 51
mu_1 = 42
sigma_1 = 5
stichprobe = np.random.normal(mu_1, sigma_1, stichprobenumfang)
# 
# Aufgabenstellung Hypothesentest
alpha = 0.05
print('pruefe auf einem Signifikanzniveau von alpha = ', alpha*100, '% die Hypothese')
print('H0: Die Stichprobe ist normalverteilt')
#
# Durchfuehrung des Kolmogorow-Smirnow-Tests:
#
# Berechnung der kumulativen Normalverteilung P und der empirischen
# kumulativen Verteilung der Stichprobe H
P, H, xsort, xbar, s = cumulatives(stichprobe)
#
# Berechnung der Pruef- oder Testgroesse
Testgroesse = np.max(np.abs(H - P))
#
# Berechnung des Schwellwertes 
K_alphaJ = np.sqrt(-np.log(alpha/2)/(2*stichprobenumfang))
#
# Antwortsatz
if (Testgroesse > K_alphaJ):
    print("{:1.2f}".format(Testgroesse), ' > ', "{:1.2f}".format(K_alphaJ), ' ==> Die Nullhypothese H0 wird verworfen')
else:
    print("{:1.2f}".format(Testgroesse), ' < oder = ', "{:1.2f}".format(K_alphaJ), ' ==> Die Nullhypothese H0 wird akzeptiert')
#
# graphische Darstellung
plt.figure(1)
plt.plot(xsort, P, 'r-')
plt.plot(xsort, H, 'kd')
plt.grid()
plt.xlabel('Zufallsgroesse')
plt.ylabel('Wahrscheinlichkeit')
plt.show()