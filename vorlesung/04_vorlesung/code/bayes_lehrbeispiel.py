"""Reproduzierbare Abbildungen für Vorlesung 4, Teil 2.

Lehrdaten und Priorannahmen sind eigens für dieses Beispiel festgelegt.
Alle Zahlenwerte werden in Ohm eingesetzt. Die Dichte von mu hat daher
Einheit 1/Ohm, die gemeinsame Dichte von (mu,sigma) Einheit 1/Ohm**2.
Benötigt: Python, NumPy, SciPy, Matplotlib.
Aufruf: python bayes_lehrbeispiel.py

Die Integration über mu erfolgt analytisch bedingt auf sigma; über sigma
wird eine Gauss-Legendre-Quadratur verwendet. Damit werden die
Verteilungen über ihren gesamten Definitionsbereich normiert, nicht
anhand einer endlichen Bildfläche oder einer Summe von Pixelwerten.
"""
from pathlib import Path
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ndtr
from scipy.stats import norm, t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
MEDIA = ROOT / '04_vorlesung' / 'media'
MEDIA.mkdir(parents=True, exist_ok=True)
x = np.array([32., 38., 41., 45., 48., 53., 58.])
J = len(x)
xbar = float(x.mean())
s = float(x.std(ddof=1))
sse = float(np.sum((x-xbar)**2))
y0, u0 = 44., 5.
sigma_min, sigma_max = 2., 20.


def conditional(sigma):
    """Posterior-Normalverteilung von mu bei festem sigma."""
    v = 1. / (1./u0**2 + J/sigma**2)
    mean = v*(y0/u0**2 + J*xbar/sigma**2)
    return mean, np.sqrt(v)


def log_weight(sigma):
    """Logarithmische marginale Posteriordichte von sigma bis auf Konstante.

    Integral der Normal-Likelihood mal Normalprior von mu. Die konstante
    Uniformprior von sigma kürzt sich bei der Normierung heraus.
    """
    variance = u0**2 + sigma**2/J
    return (-(J-1)*np.log(sigma)-sse/(2*sigma**2)
            - .5*np.log(variance)-.5*(xbar-y0)**2/variance)


shift = float(np.max(log_weight(np.linspace(sigma_min,sigma_max,1000))))
normalizer, quad_error = quad(lambda z: np.exp(log_weight(z)-shift),
                            sigma_min,sigma_max,epsabs=1e-12,epsrel=1e-12)


def mixture(n):
    nodes, weights = leggauss(n)
    sigma = sigma_min + (nodes+1)*(sigma_max-sigma_min)/2
    weights = weights*(sigma_max-sigma_min)/2
    probabilities = weights*np.exp(log_weight(sigma)-shift)/normalizer
    means, widths = conditional(sigma)
    return sigma, probabilities, means, widths


sigma_nodes, weights, means, widths = mixture(320)
assert abs(weights.sum()-1) < 1e-10


def pdf(mu):
    mu = np.asarray(mu)
    return np.sum(weights*norm.pdf(mu[...,None],means,widths),axis=-1)


def cdf(mu):
    mu = np.asarray(mu)
    return np.sum(weights*ndtr((mu[...,None]-means)/widths),axis=-1)


mean = float(np.sum(weights*means))
variance = float(np.sum(weights*(widths**2+(means-mean)**2)))
std = float(np.sqrt(variance))
low,high = [brentq(lambda z:float(cdf(z))-alpha,-100,200,xtol=1e-12)
            for alpha in [.025,.975]]
scale = s/np.sqrt(J)
ci = xbar + t.ppf([.025,.975],J-1)*scale
# Independent adaptive integration checks the complete marginal density.
area,error = quad(lambda z:float(pdf(z)),-np.inf,np.inf,epsabs=1e-10)
assert abs(area-1)<1e-9
assert abs(float(cdf(high)-cdf(low))-.95)<1e-10
assert abs(t.cdf((ci[1]-xbar)/scale,J-1)
           -t.cdf((ci[0]-xbar)/scale,J-1)-.95)<1e-12
_,w160,m160,s160 = mixture(160)
assert abs(np.sum(w160*m160)-mean)<1e-10
assert abs(np.sum(w160*(s160**2+(m160-mean)**2))-variance)<1e-10
results = {'data_ohm':x.tolist(),'J':J,'sample_mean_ohm':xbar,
           'sample_sd_ohm':s,'prior_mu_mean_ohm':y0,'prior_mu_sd_ohm':u0,
           'prior_sigma_uniform_ohm':[sigma_min,sigma_max],
           'posterior_mean_ohm':mean,'posterior_sd_ohm':std,
           'credible_interval_95_ohm':[low,high],
           'confidence_interval_95_ohm':ci.tolist(),
           'posterior_density_integral':area,
           'credible_interval_probability':float(cdf(high)-cdf(low)),
           'quadrature_nodes':320}
(ROOT/'bayes_lehrbeispiel_ergebnisse.json').write_text(
    json.dumps(results,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
print(json.dumps(results,indent=2,ensure_ascii=False))

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'pdf.fonttype':42,'savefig.bbox':'tight'})
mu_grid=np.linspace(25,65,900)
sigma_grid=np.linspace(sigma_min,sigma_max,450)
mu_cond,sd_cond=conditional(sigma_grid)
sigma_density=np.exp(log_weight(sigma_grid)-shift)/normalizer
joint=sigma_density[:,None]*norm.pdf(mu_grid[None,:],mu_cond[:,None],sd_cond[:,None])
blue='#1464a0'; red='#ba4536'

fig,ax=plt.subplots(figsize=(7.2,4.8),layout='constrained')
img=ax.pcolormesh(mu_grid,sigma_grid,joint,shading='auto',cmap='viridis',rasterized=True)
ax.set(xlabel=r'$\mu\;/\;\Omega$',ylabel=r'$\sigma\;/\;\Omega$')
cb=fig.colorbar(img,ax=ax)
cb.set_label(r'$p_1(\mu,\sigma\mid\mathbf{x})\;/\;\Omega^{-2}$')
fig.savefig(MEDIA/'understand_bayes_mean_posteriormatrix_ueberarbeitet.pdf')
plt.close(fig)

fig,ax=plt.subplots(figsize=(5.5,4.3),layout='constrained')
ax.plot(mu_grid,pdf(mu_grid),color=blue,label='Posterior-Dichte',lw=2)
ax.plot(mu_grid,t.pdf((mu_grid-xbar)/scale,J-1)/scale,color=red,
        ls='--',lw=1.8,label=r'$t$-Referenzkurve $g_t$')
fill=np.linspace(low,high,400)
ax.fill_between(fill,pdf(fill),color=blue,alpha=.13)
for b in [low,high]: ax.axvline(b,color=blue,lw=1.2)
for b in ci: ax.axvline(b,color=red,ls='--',lw=1.2)
ax.set(xlim=(25,65),ylim=(0,None),xlabel=r'$\mu\;/\;\Omega$',
       ylabel=r'Dichte $/\;\Omega^{-1}$')
ax.legend(fontsize=9,loc='upper right',frameon=False)
fig.savefig(MEDIA/'understand_bayes_mean_posteriormarginal_ueberarbeitet.pdf')
plt.close(fig)

fig,ax=plt.subplots(figsize=(5.5,4.3),layout='constrained')
ax.plot(mu_grid,cdf(mu_grid),color=blue,lw=2)
for value,alpha in [(low,.025),(high,.975)]:
    ax.hlines(alpha,25,value,color=red,linestyle='--',lw=1.2)
    ax.vlines(value,0,alpha,color=red,linestyle='--',lw=1.2)
    ax.plot(value,alpha,'o',color=red,ms=4)
ax.set(xlim=(25,65),ylim=(0,1),xlabel=r'$m\;/\;\Omega$',ylabel=r'$P_1(m)$')
ax.set_yticks([.025,.25,.5,.75,.975])
ax.set_yticklabels(['0,025','0,25','0,5','0,75','0,975'])
ax.tick_params(axis='y',labelsize=9)
fig.savefig(MEDIA/'understand_bayes_mean_cumposterior_ueberarbeitet.pdf')
plt.close(fig)
