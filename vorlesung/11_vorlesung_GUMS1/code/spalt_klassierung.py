import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

def plot_distributions(
        quantity, 
        cdf, 
        bar_quanti, 
        quanti_interval, 
        str_abcissa):
    plt.figure()
    plt.plot(quantity, cdf, 'k-', lw=1.5)
    plt.plot([bar_quanti, bar_quanti], [0, 1], 'r-')
    plt.plot([quanti_interval[0], quanti_interval[1]], [0.025, 0.975], 'r+', ms=8)
    plt.grid()
    plt.xlabel(str_abcissa,fontsize=14)
    plt.ylabel('kumulative Wahrscheinlichkeit',fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig('cdf_spaltbreite.svg', bbox_inches='tight')
    plt.show()

    mass, edges = np.histogram(quantity, bins=530)
    print(len(quantity), np.sum(mass))
    pdf = mass/np.sum(mass)
    edge_mid = 0.5 * (edges[1:] + edges[0:-1])
    plt.figure()
    plt.plot(edge_mid, pdf)
    plt.plot([bar_quanti, bar_quanti], [0, np.max(pdf)], 'r-')

    i = np.searchsorted(edge_mid, quanti_interval[0], side="left")  
    print('plot', i, edge_mid[0], edge_mid[-1]) 
    print(edge_mid[i], quanti_interval[0]) 
    plt.plot([quanti_interval[0], quanti_interval[0]], [0.0, pdf[i]], 'r--')

    i = np.searchsorted(edge_mid, quanti_interval[1], side="left")   
    print(edge_mid[i], quanti_interval[1]) 
    plt.plot([quanti_interval[1], quanti_interval[1]], [0.0, pdf[i]], 'r--')
    plt.grid()
    plt.xlabel(str_abcissa,fontsize=14)
    plt.ylabel('Wahrscheinlichkeitsdichte',fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig('pdf_spaltbreite.svg', bbox_inches='tight')
    plt.show()    

def klassierung_spalt(K=86,
                     xmin_0=386.3, sigma_xmin = 5.9,
                     lambda_0=0.5500, sigma_lambda=0.0020,
                     L_0=4000.0, a_L=20.0, 
                     alpha=0.05):
    # Quantile in der Mitte der Wahrscheinlichkeitsklassen
    r = (np.arange(K) + 0.5) / K
    z = norm.ppf(r)
    lambda_rep = lambda_0 + sigma_lambda * z
    L_rep = L_0 - a_L * np.cos(np.pi * r)
    xmin_rep = xmin_0 + sigma_xmin * z
    if np.any(xmin_rep <= 0.0):
        raise ValueError("Positive Minimumspositionen erforderlich")

    # Alle K**3 Kombinationen, keine Zufallsziehung
    # berechne alle Produkte von lambda und L
    numerator = lambda_rep[:, None] * L_rep[None, :]
    # und alle Produkte von numerator und xmin_rep
    # diese sortieren fuer die kumulative Verteilung
    d_all = np.sort((numerator[:, :, None]
             / xmin_rep[None, None, :]).ravel())
    
    # erzeuge ein Array von derselben Groesse n wie d_all,
    # in dem alle Arrayelemente denselben Wert haben,
    # hier also 1/n
    # in diesem Fall ist weights schon so gewaehlt,
    # dass die Summe ueber alle weights schon 1 ist
    weights = np.full(d_all.size, 1.0 / d_all.size)

    # Momente direkt aus den gewichteten Modellwerten
    d_est = np.sum(weights * d_all)
    u_d = np.sqrt(np.sum(weights * (d_all - d_est)**2))

    # # Klassenwahrscheinlichkeiten und Dichte getrennt berechnen,
    # mass, edges = np.histogram(d_all, bins=n_bins,
    #                            weights=weights, density=False)
    # mass /= mass.sum()  # Ausgleich numerischer Rundung
    # widths = np.diff(edges)
    # density = mass / widths
    # cdf_edges = np.r_[0.0, np.cumsum(mass)]

    N_d = len(d_all)
    cdf_slit = (np.arange(N_d) + 0.5) / N_d

    def quantile(probability):
        # Finde die Stelle, in deren Naehe das gesuchte Quantil liegt.
        i = np.searchsorted(
            cdf_slit[0:-1], probability, side="left"
        )

        # lineare Interpolation zwischen 2 Stützstellen
        probability_before = cdf_slit[i]
        probability_after = cdf_slit[i+1]
        fraction = (probability - probability_before) / (probability_after - probability_before)
        value = d_all[i] + fraction * (d_all[i+1] - d_all[i])

        return float(value)


    interval = (quantile(alpha / 2.0),
                quantile(1.0 - alpha / 2.0))
    at50percent = quantile(0.5)
    return {"d": d_est, "u": u_d, "interval": interval,
            "d_all": d_all, "d_cdf": cdf_slit,
            "at50percent": at50percent}

if __name__ == "__main__":
    result = klassierung_spalt()
    print(result["d"], 'pm', result["u"])
    print(result["d"]-1.96*result["u"], result["d"]+1.96*result["u"])
    print(result["interval"])
    print('d at 50% probability', result["at50percent"])
    print((result["d"]-result["interval"][0])/result["u"])
    print((result["interval"][1]-result["d"])/result["u"])

    plot_distributions(
        result["d_all"], result["d_cdf"],
        result["d"], result["interval"],
        r'Spaltbreite $d$ / µm'
    )

