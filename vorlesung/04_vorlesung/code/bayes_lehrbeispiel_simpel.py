import numpy as np
import matplotlib.pyplot as plt


x = np.array([32., 38., 41., 45., 48., 53., 58.])
J = len(x)
xbar = float(x.mean())
s = float(x.std(ddof=1))
sse = float(np.sum((x-xbar)**2))
y0, u0 = 44., 5.
sigma_min, sigma_max = 2., 20.

mu_prime = np.linspace(25,65,15000)
sigma_prime = np.linspace(sigma_min,sigma_max,15000)

d_mu = mu_prime[1] - mu_prime[0]
d_sigma = sigma_prime[1] - sigma_prime[0]

mu_grid, sigma_grid = np.meshgrid(mu_prime, sigma_prime)

pi_0 = np.exp(-0.5 * ((mu_grid - y0) / u0)**2) / (np.sqrt(2*np.pi) * u0)

J = len(x)
S_beob = 0
for j in range(0, J):
    S_beob += (x[j] - mu_grid)**2

likeli = np.exp(-0.5 * S_beob / sigma_grid**2) / (np.sqrt(2*np.pi) * sigma_grid)**J 

pi_1 = likeli * pi_0 

Crezi = np.sum(pi_1) * d_mu * d_sigma

print(Crezi)

pi_1 /= Crezi

pi1_margi = np.sum(pi_1, axis=0) * d_sigma

P_1 = np.cumsum(pi1_margi) * d_mu


y_1 = np.sum(mu_prime*pi1_margi) * d_mu
print('y_1', y_1, 'np.sum(pi1_margi) * d_mu', np.sum(pi1_margi) * d_mu)
var_1 = np.sum(((mu_prime - y_1)**2) * pi1_margi) * d_mu
u_1 = np.sqrt(var_1)
print('u_1', u_1)
g_1 = np.exp(-0.5 * (mu_prime - y_1)**2 / var_1) / (np.sqrt(2*np.pi) * u_1)

i_links = np.argmin(np.abs(P_1-0.025))
ger_li = np.polyfit(P_1[i_links-3:i_links+4], mu_prime[i_links-3:i_links+4], 2)
y_min = (ger_li[0] * 0.025 + ger_li[1]) * 0.025 + ger_li[2]
i_rechts = np.argmin(np.abs(P_1-0.975))
ger_re = np.polyfit(P_1[i_rechts-3:i_rechts+4], mu_prime[i_rechts-3:i_rechts+4], 2)
y_max = (ger_re[0] * 0.975 + ger_re[1]) * 0.975 + ger_re[2]

print(y_min, y_max)

plt.figure()
plt.plot(P_1[i_links-3:i_links+4], mu_prime[i_links-3:i_links+4], 'o')
ppp = np.linspace(P_1[i_links-3],P_1[i_links+3],20)
plt.plot(ppp, (ger_li[0] * ppp + ger_li[1]) * ppp + ger_li[2], 'r-')
plt.plot(0.025, y_min, 'kd', mfc='None')

if 0:
    plt.figure(figsize=(7.2,4.8),layout='constrained')
    plt.imshow(
        pi_1,
        extent=[mu_prime[0],mu_prime[-1], sigma_min, sigma_max],
        origin='lower',
        aspect=2
    )
    cb=plt.colorbar()
    cb.set_label(r'$\pi_1(\mu,\sigma\mid\mathbf{x})\;/\;\Omega^{-2}$')
    plt.xlabel(r'$\mu\;/\;\Omega$')
    plt.ylabel(r'$\sigma\;/\;\Omega$')

plt.figure()
plt.plot(mu_prime, pi1_margi, 'kx', ms=4)
plt.plot(mu_prime, g_1, 'r+', ms=4)
plt.grid()
plt.xlabel(r'$\mu\;/\;\Omega$')
plt.ylabel(r'Posteriordichte $/\;\Omega^{-1}$')


with open('P1chatGPT.npy','rb') as fp:
    mu_cgpt = np.load(fp)
    P1_cgpt = np.load(fp) 


plt.figure()
plt.plot(mu_prime, P_1, 'kx', ms=4)
plt.plot(mu_cgpt, P1_cgpt, 'b+', ms=4)
plt.plot([mu_prime[0], mu_prime[-1]], [0.025, 0.025], 'r--')
plt.plot([mu_prime[0], mu_prime[-1]], [0.975, 0.975], 'r--')
plt.grid()
plt.xlabel(r'$\mu\;/\;\Omega$')
plt.ylabel(r'$P_1 \; /\;\Omega^{-1}$')

plt.show()



