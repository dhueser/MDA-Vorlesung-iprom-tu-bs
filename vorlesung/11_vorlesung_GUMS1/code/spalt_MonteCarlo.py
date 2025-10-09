import numpy as np 
from scipy.optimize import leastsq

if __name__ == "__main__":
    n_rnd = 1500
    n_spalt = 200
    L = 4000 # micron
    a = 20 # micron
    lambda_light = 0.55 # micron
    sigma_lambda = 0.002
    num_det = 64
    intens_0 = 240
    d_slit = 5.75 # micron
    d_nom = 5.5
    pix_noise = 0.1
    x_min = L * np.tan(np.arcsin(lambda_light/ d_slit))
    x_det = np.linspace(-2.5*x_min, 2.5*x_min, num_det) 
    sigma_dx = 0.1 * (x_det[1] - x_det[0])

    np.random.seed(0)
    dslit_esti = np.zeros(n_rnd)
    sigma_dslit = np.zeros(n_rnd)

    funcBeugungsbild = lambda p, r : p[0] * (np.sin(p[1]*r) / (p[1]*r))**2
    resiBeugungsbild = lambda p, r, g : g - funcBeugungsbild(p, r)

    x_simu = x_det + np.random.normal(0, sigma_dx, num_det)
    L_vary = L + a * np.cos(np.pi * np.random.uniform(0, 1, n_rnd))
    lambda_vary = np.random.normal(lambda_light, sigma_lambda, n_rnd)
    for i_rnd in range(0, n_rnd):
        r_simu = (np.pi / lambda_vary[i_rnd]) * np.sin(np.arctan(x_simu/L_vary[i_rnd]))
        u_simu =  d_slit * r_simu
        # Intensitaetsverteilung der Beugungsinterferenzen
        intens = intens_0 * (np.sin(u_simu)/u_simu)**2
        # Pixelrauschen
        detnoise = np.random.normal(0,pix_noise,num_det)
        detsignal = intens * (1+detnoise)
        p_init = np.array([np.max(detsignal), d_nom])
        p_fit, inv_Jac, infodict, errmsg, ier = leastsq(resiBeugungsbild, p_init, \
                                    args=(r_simu, detsignal), full_output=True)
        # Varianz der Residuen
        residual_variance = np.sum(infodict['fvec']**2) / (num_det - 2)
        # Kovarianzmatrix
        covar = inv_Jac * residual_variance
        sigma_dslit[i_rnd] = np.sqrt(covar[1][1])
        dslit_esti[i_rnd] = p_fit[1]

    d_mean = np.mean(dslit_esti)
    d_sort = np.sort(dslit_esti)
    prob_d_kum = np.arange(0, n_rnd) / n_rnd
    i_li0 = np.argmin(np.abs(prob_d_kum-0.025))
    i_re0 = np.argmin(np.abs(prob_d_kum-0.975))
    print('Ergebnis: bar d', d_mean)
    print('Ueberdeckungsintervall: [', d_sort[i_li0],',', d_sort[i_re0],']')
    print(np.mean(sigma_dslit), np.std(sigma_dslit, ddof=1))
