import numpy as np
import matplotlib.pyplot as plt

flagpltconsumer = 0
savethefig = 1

gaussian = lambda x, mu, sigm : \
  np.exp(-0.5*((x - mu) / sigm)**2) / (sigm * np.sqrt(2*np.pi))

step = 0.001
# ----
y0 = 1500
s0 = 0.12
T_L = 1499.80
T_U = 1500.20
A_L = 1499.82
A_U = 1500.18
s1 = (T_U - T_L) / (2 * 1.96)
mininf = T_L - 3*s0
maxinf = T_U + 3*s0

#  s0 = s1;
y = np.arange(T_L, T_U+step, step)
prior = gaussian(y, y0, s0)

# inspection device - instrument uncertainty
s_insp = 0.04


yall = np.arange(mininf,maxinf+step,step)
resistor_distri = gaussian(yall, y0, s0)
p_T = gaussian(np.array([T_L, T_U]), y0, s0)
p_A = gaussian(np.array([A_L, A_U]), y0, s0)

if (flagpltconsumer == 1):
    x1 = T_U + 0.03
    x = np.arange(x1-5*s_insp,step+x1+5*s_insp,step)
    xconsumerrisk = np.arange(x1-5*s_insp,step+A_U, step)
    p_consumerrisk = gaussian(xconsumerrisk, x1, s_insp) * gaussian(x1, y0, s0)
else:
    x1 = A_U - 0.03
    x = np.arange(x1-5*s_insp,step+x1+5*s_insp,step)
    xproducerrisk = np.arange(A_U, step+x1+5*s_insp,step)
    p_producerrisk = gaussian(xproducerrisk, x1, s_insp) * gaussian(x1, y0, s0)

print(x1)

p_insp = gaussian(x, x1, s_insp) * gaussian(x1, y0, s0)
p_x1 = gaussian(x1, x1, s_insp) * gaussian(x1, y0, s0)
p_AU = gaussian(A_U, x1, s_insp) * gaussian(x1, y0, s0)

plt.figure(301);
plt.fill_between(y, np.zeros(len(y)),\
 prior, color='yellow', alpha=0.4)
plt.plot(yall, resistor_distri, 'k-', label='prior '+r'$p(y)$')
plt.plot([T_L, T_L], [0, p_T[0]], 'r-')
plt.plot([T_U, T_U], [0, p_T[1]], 'r-')
plt.plot([A_L, A_L], [0, p_A[0]], 'k--')
plt.plot([A_U, A_U], [0, p_AU], 'k--')
plt.plot(x, p_insp, 'b-', label='posterior '+r'$p(y) \cdot p(x|y=x_1)$')

if (flagpltconsumer == 1):
    plt.fill_between(xconsumerrisk, np.zeros(len(xconsumerrisk)),\
      p_consumerrisk, color='red', alpha=0.3)
else:
    plt.fill_between(xproducerrisk, np.zeros(len(xproducerrisk)),\
      p_producerrisk, color='red', alpha=0.3)
plt.plot([x1, x1], [0, p_x1], 'b-')
plt.grid()
plt.xlim(-0.5+1.5e3, 0.5+1.5e3)

if (flagpltconsumer == 1):
    plt.ylim(0, 6)
else:
    plt.ylim(0, 16)

plt.xlabel(r'$R \; / \; \Omega$')
plt.ylabel(r'$p \; / \; \frac{1}{\Omega}$')
plt.legend()

if (savethefig == 1):
    if (flagpltconsumer == 1):
        plt.savefig('Konsumentenrisiko_x0p03.svg')
    else:
        plt.savefig('Produzentenrisiko_x0p03.svg')


plt.show()
