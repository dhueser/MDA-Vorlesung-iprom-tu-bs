import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import stats
Jz = 1000
nu = 45
R = np.random.rand(1, Jz)
Z = scipy.stats.t.ppf(R[0], nu)
plt.figure(1)
plt.plot(R[0], Z,'.')
plt.show()