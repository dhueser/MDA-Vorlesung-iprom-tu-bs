import numpy as np 
from scipy.special import erfinv
from scipy.optimize import leastsq

funcinvP = lambda p, x : np.sqrt(2.0) * p[1] * erfinv(2.0 * x - 1.0) + p[0]
funcGauss = lambda p, x : np.exp(-0.5 * ((x - p[0]) / p[1])**2) \
                          / (np.sqrt(2.0*np.pi) * p[1])

def gen_sample(param_in, prob_in):
    n_gr = len(prob_in)
    measurand_gr = funcinvP(param_in, prob_in)
    measurand_mid = 0.5 * (measurand_gr[0:n_gr-1] + measurand_gr[1:n_gr])
    d_measurand = measurand_gr[1:n_gr] - measurand_gr[0:n_gr-1]
    return measurand_gr, measurand_mid, d_measurand

if __name__ == "__main__":
    n_sample = 87
    n_MC = 31
    n_spalt = 250
    L = 4000 # micron
    a = 20 # micron
    lambda_light = 0.55 # micron
    num_det = 64
    intens_0 = 240
    d_slit = 5.75 # micron
    d_nom = 5.5
    pix_noise = 0.1
    x_min = L * np.tan(np.arcsin(lambda_light/ d_slit))
    x_det = np.linspace(-2.5*x_min, 2.5*x_min, num_det) 
    sigma_dx = 0.1 * (x_det[1] - x_det[0])

    d_p = 1 / (n_sample+1)
    probab = np.arange(d_p, d_p*(n_sample+1), d_p)
    prob_mid = 0.5 * (probab[0:n_sample-1] + probab[1:n_sample])

    p_light = np.array([lambda_light, 0.002])
    lambda_gr, lambda_mid, d_lambda = gen_sample(p_light, probab)
    L_gr = L - a * np.cos(np.pi * probab)
    L_mid = 0.5 * (L_gr[0:n_sample-1] + L_gr[1:n_sample])
    d_L = L_gr[1:n_sample] - L_gr[0:n_sample-1]

    prob_L = np.gradient(prob_mid, L_mid)
    prob_light = funcGauss(p_light, lambda_mid)

    prob_L_lambda = []
    L_all = []
    lambda_all = []
    for i_L in range(0, n_sample-1):
        prob_L_lambda = np.append(prob_L_lambda, \
                                  prob_L[i_L]*prob_light*d_L[i_L]*d_lambda)
        L_all = np.append(L_all, L_mid[i_L]*np.ones(n_sample-1))
        lambda_all = np.append(lambda_all, lambda_mid)
    n_all = len(prob_L_lambda)
    print(n_all, (n_sample-1)**2)

    np.random.seed(0)
    dslit_esti = np.zeros((n_MC, n_all))
    sigma_dslit = np.zeros((n_MC, n_all))
    funcBeugungsbild = lambda p, r : p[0] * (np.sin(p[1]*r) / (p[1]*r))**2
    resiBeugungsbild = lambda p, r, g : g - funcBeugungsbild(p, r)
    for i_MC in range(0, n_MC):
        x_simu = x_det + np.random.normal(0, sigma_dx, num_det)
        for i_all in range(0, n_all):
            r_simu =  np.sin(np.arctan(x_simu/L_all[i_all])) * np.pi/lambda_all[i_all]
            u_simu = d_slit * r_simu
            intens = intens_0 * (np.sin(u_simu)/u_simu)**2
            detnoise = np.random.normal(0,pix_noise,num_det)
            detsignal = intens * (1+detnoise)
            p_init = np.array([np.max(detsignal), d_nom])
            p_fit, inv_rmatrix, infodict, errmsg, ier = leastsq(resiBeugungsbild, \
                                p_init, args=(r_simu, detsignal), full_output=True)
            residual_variance = np.sum(infodict['fvec']**2) / (num_det - 2)
            covar = inv_rmatrix * residual_variance
            sigma_dslit[i_MC][i_all] = np.sqrt(covar[1][1])
            dslit_esti[i_MC][i_all] = p_fit[1]

    d_spalt = np.linspace(np.min(dslit_esti), np.max(dslit_esti), n_spalt)
    d_d = d_spalt[1] - d_spalt[0]
    prob_d = np.zeros(n_spalt)
    for i_d in range(0, n_spalt-1):
        idx_list = np.where((d_spalt[i_d] <= dslit_esti[i_MC]) & \
                            (dslit_esti[i_MC] < d_spalt[i_d+1]))[0]
        if len(idx_list) > 0:
            prob_d[i_d] += np.sum(prob_L_lambda[idx_list])
    idx_list = np.where(d_spalt[n_spalt-1] <= dslit_esti[i_MC])[0]
    if len(idx_list) > 0:
        prob_d[n_spalt-1] += np.sum(prob_L_lambda[idx_list])

    prob_d_kum = np.cumsum(prob_d)
    prob_d /= (d_d*prob_d_kum[-1])
    prob_d_kum /= prob_d_kum[-1]
    i_li = np.argmin(np.abs(prob_d_kum-0.025))
    i_re = np.argmin(np.abs(prob_d_kum-0.975))
    d_erw = np.sum(prob_d * d_spalt) * d_d
    print('Ergebnis: bar d', d_erw)
    print('Ueberdeckungsintervall: [', d_spalt[i_li],',', d_spalt[i_re],']')
    print(np.mean(sigma_dslit), np.std(sigma_dslit, ddof=1))