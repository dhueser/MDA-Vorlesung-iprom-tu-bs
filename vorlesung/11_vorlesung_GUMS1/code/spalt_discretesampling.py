import numpy as np 
from scipy.special import erfinv

funcinvP = lambda p, x : np.sqrt(2.0) * p[1] * erfinv(2.0 * x - 1.0) + p[0]
funcGauss = lambda p, x : np.exp(-0.5 * ((x - p[0]) / p[1])**2) / \
                          (np.sqrt(2.0*np.pi) * p[1])

def gen_sample(param_in, prob_in):
    n_gr = len(prob_in)
    measurand_gr = funcinvP(param_in, prob_in)
    measurand_mid = 0.5 * (measurand_gr[0:n_gr-1] + measurand_gr[1:n_gr])
    d_measurand = measurand_gr[1:n_gr] - measurand_gr[0:n_gr-1]
    return measurand_gr, measurand_mid, d_measurand

if __name__ == "__main__":
    n_gauss = 87
    n_spalt = 250
    L = 4000 # micron
    a = 20 # micron
    lambda_light = 0.55 # micron
    x_min_0 = 384.3711
    
    d_p = 1 / (n_gauss+1)
    probab = np.arange(d_p, d_p*(n_gauss+1), d_p)
    prob_mid = 0.5 * (probab[0:n_gauss-1] + probab[1:n_gauss])
    print(probab)

    p_light = np.array([lambda_light, 0.002]) 
    p_xmin = np.array([x_min_0, 6.5])

    lambda_gr, lambda_mid, d_lambda = gen_sample(p_light, probab)
    xmin_gr, xmin_mid, d_xmin = gen_sample(p_xmin, probab)
    L_gr = L - a * np.cos(np.pi * probab)
    L_mid = 0.5 * (L_gr[0:n_gauss-1] + L_gr[1:n_gauss])
    d_L = L_gr[1:n_gauss] - L_gr[0:n_gauss-1]

    prob_L = np.gradient(prob_mid, L_mid)
    prob_light = funcGauss(p_light, lambda_mid)
    prob_xmin = funcGauss(p_xmin, xmin_mid)

    nominator = []
    prob_nom = []
    for i_L in range(0, n_gauss-1):
        nominator = np.append(nominator, L_mid[i_L]*lambda_mid)
        prob_nom = np.append(prob_nom, prob_L[i_L]*prob_light*d_L[i_L]*d_lambda)

    d_spalt_all = []
    prob_d_all = []
    for i_x in range(0, n_gauss-1):
        d_spalt_all = np.append(d_spalt_all, nominator/xmin_mid[i_x])
        prob_d_all = np.append(prob_d_all, prob_nom*prob_xmin[i_x]*d_xmin[i_x])

    d_spalt = np.linspace(np.min(d_spalt_all), np.max(d_spalt_all), n_spalt)
    d_d = d_spalt[1] - d_spalt[0]
    prob_d = np.zeros(n_spalt)
    for i_d in range(0, n_spalt-1):
        idx_list = np.where((d_spalt[i_d] <= d_spalt_all) & \
                            (d_spalt_all < d_spalt[i_d+1]))[0]
        if len(idx_list) > 0:
            prob_d[i_d] = np.sum(prob_d_all[idx_list])
    idx_list = np.where(d_spalt[n_spalt-1] <= d_spalt_all)[0]
    if len(idx_list) > 0:
        prob_d[n_spalt-1] = np.sum(prob_d_all[idx_list])

    d_erw = np.sum(prob_d * d_spalt) * d_d
    prob_d_kum = np.cumsum(prob_d)
    prob_d /= (d_d*prob_d_kum[-1])
    prob_d_kum /= prob_d_kum[-1]
    i_li = np.argmin(np.abs(prob_d_kum-0.025))
    i_re = np.argmin(np.abs(prob_d_kum-0.975))
    print('Ergebnis: bar d', d_erw)
    print('Ueberdeckungsintervall: [', d_spalt[i_li],',', d_spalt[i_re],']')