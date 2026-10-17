"""Einfache Demonstration des KS- und Lilliefors-Tests.

Benoetigte Pakete:
    numpy
    scipy
    matplotlib

Installation, falls erforderlich:
    python -m pip install numpy scipy matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm


def pruefgroesse_D(x, mu=None, sigma=None):
    """Berechnet D und die Werte der empirischen Verteilungsfunktion.

    Die Daten werden sortiert. Fuer den j-ten sortierten Wert gilt

        H(x_(j)) = j / J.

    Die Funktion verwendet den Mittelwert und die empirische
    Standardabweichung der Stichprobe, um eine Normalverteilung
    an die Daten anzupassen.
    """
    # Stichprobe sortieren: x_1, ..., x_J
    x_sortiert = np.sort(x)
    J = len(x_sortiert)

    # In NumPy muss ddof=1 angegeben werden, damit durch J-1 geteilt wird.
    mittelwert = np.mean(x_sortiert)
    standardabweichung = np.std(x_sortiert, ddof=1)

    # Empirische kumulierte relative Haeufigkeit H(x_j) = j/J
    rang = np.arange(1, J + 1)
    H = rang / J

    # Angepasste Normalverteilung: F(x) = Phi((x - y) / s)
    if mu and sigma:
        print(f"Verwende festgelegte Verteilungsparameter mit  mu2 = {mu}, s2 = {sigma}")
        F = norm.cdf(x_sortiert, loc=mu, scale=sigma)
    else:
        print("Verwende empirisch geschaetzte Parameter der Verteilung")
        F = norm.cdf(x_sortiert, loc=mittelwert, scale=standardabweichung)

    # Beide Seiten jeder Stufe der empirischen Verteilungsfunktion pruefen.
    D_plus = np.max(H - F)
    D_minus = np.max(F - (rang - 1) / J)
    D = max(D_plus, D_minus)

    return D, x_sortiert, H, F, mittelwert, standardabweichung


def lilliefors_grenzwert(J, alpha=0.05):
    """Naeherungswert fuer den Lilliefors-Grenzwert.

    Fuer alpha=0.05 gilt fuer groessere Stichproben naeherungsweise

        c_(alpha,J)^L = 0.886 / sqrt(J).

    Die Formel ist eine Naeherung und sollte hier nur fuer J >= 30
    verwendet werden. Fuer kleine Stichproben sind tabellierte Werte
    oder eine Simulation geeigneter.
    """
    if alpha != 0.05:
        raise ValueError("Diese einfache Funktion ist nur fuer alpha=0.05 programmiert.")
    if J < 30:
        print(
            "Hinweis: Fuer J < 30 ist der Lilliefors-Grenzwert "
            "nur eine grobe Naeherung."
        )
    return 0.886 / np.sqrt(J)


def ks_grenzwert(J, alpha=0.05):
    """Asymptotischer KS-Grenzwert fuer eine vollstaendig festgelegte
    Verteilung F0, d.h. Erwartungswert und Standardabweichung sind
    vordefiniert und nicht geschaetzt, ausserdem fuer J >= 30.
    """
    return np.sqrt(-np.log(alpha / 2) / (2 * J))


def ergebnis_ausgeben(name, D, grenzwert):
    """Gibt das Testergebnis aus."""
    print(f"\n{name}")
    print("-" * len(name))
    print(f"Pruefgroesse D = {D:.4f}")
    print(f"Grenzwert   = {grenzwert:.4f}")

    if D > grenzwert:
        print("Ergebnis: H0 wird verworfen.")
    else:
        print("Ergebnis: H0 wird nicht verworfen.")


# ---------------------------------------------------------------------------
# Teil 1: Eine normalverteilte Stichprobe
# ---------------------------------------------------------------------------
print("TEIL 1: Eine Stichprobe")

alpha = 0.05
rng = np.random.default_rng(2026)  # fester Startwert: reproduzierbare Zahlen

mu1 = 23.0
s1 = 3.0
J1 = 200

# Zufallszahlen aus N(mu1, s1^2) erzeugen
Xarray = rng.normal(loc=mu1, scale=s1, size=J1)

D1, Xsortiert, H1, F1, y1, s_emp1 = pruefgroesse_D(Xarray)
c_lilliefors = lilliefors_grenzwert(J1, alpha)

print(f"Parameter zur Erzeugung der ersten Stichprobe: mu2 = {mu1}, s2 = {s1}")
print(f"Stichprobenumfang J = {J1}")
print(f"Stichprobenmittelwert y = {y1:.3f}")
print(f"Empirische Standardabweichung s = {s_emp1:.3f}")
print(f"Lilliefors-Grenzwert c_L = {c_lilliefors:.4f}")
ergebnis_ausgeben("Lilliefors-Test", D1, c_lilliefors)

# Abbildung der beiden kumulierten Verteilungen
plt.step(Xsortiert, H1, where="post", label="empirische Verteilung H")
plt.plot(Xsortiert, F1, label="angepasste Normalverteilung F")
plt.xlabel("x")
plt.ylabel("kumulierte Hufigkeit bzw. Wahrscheinlichkeit")
plt.title("Teil 1: Empirische und angepasste Verteilungsfunktion")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


# ---------------------------------------------------------------------------
# Teil 2: Zwei erzeugte Stichproben werden zusammengelegt
# ---------------------------------------------------------------------------
print("\nTEIL 2: Vereinigte Stichproben")

mu2 = 35.0
s2 = 5.0
J2 = 31  # deutlich kleiner als J1

Yarray = rng.normal(loc=mu2, scale=s2, size=J2)

# Die beiden Stichproben werden zu einem Array vereinigt.
gemeinsames_array = np.concatenate((Xarray, Yarray))
J_gesamt = len(gemeinsames_array)

# Die folgende Auswertung behandelt den vereinigten Datensatz so,
# als kaeme er aus einer einzigen Grundgesamtheit.
D2, Xgesamt_sortiert, H2, F2, y2, s_emp2 = pruefgroesse_D(gemeinsames_array, mu1, s1)
K = ks_grenzwert(J_gesamt, alpha)

print(f"Stichprobenumfaenge: J1 = {J1}, J2 = {J2}, insgesamt J = {J_gesamt}")
print(f"Parameter zur Erzeugung der zweiten Stichprobe: mu2 = {mu2}, s2 = {s2}")
print(f"Vorgegebene Parameter der ersten Stichprobe fuer den KS-Test: mu1 = {mu1}, s1 = {s1}")

print(f"Mittelwert der vereinigten Stichproben y = {y2:.3f}")
print(f"Standardabweichung der vereinigten Stichproben s = {s_emp2:.3f}")
print(
    "Fuer den geforderten Vergleich wird der "
    "Standard-KS-Grenzwert verwendet."
)
ergebnis_ausgeben("Standard-KS-Test fuer die vereinigte Stichprobe", D2, K)

plt.step(Xgesamt_sortiert, H2, where="post", label="empirische Verteilung H")
plt.plot(Xgesamt_sortiert, F2, label="angepasste Normalverteilung F")
plt.xlabel("x")
plt.ylabel("kumulierte Haeufigkeit bzw. Wahrscheinlichkeit")
plt.title("Teil 2: Vereinigte Stichproben")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
