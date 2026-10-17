"""Loest eine Aufgabe zu zwei Stichproben und erzeugt eine LaTeX-Loesung.

Das Skript berechnet:

    a) F-Test auf gleiche Varianzen,
    b) Mittelwerte und empirische Standardabweichungen,
    c) gepoolten Zweistichproben-t-Test auf gleiche Mittelwerte,
    d) und e) einseitige Chi-Quadrat-Tests auf eine vorgegebene Varianz.

Benötigte Pakete:
    numpy
    scipy

Installation, falls erforderlich:
    python -m pip install numpy scipy

Nach dem Start wird die Datei

    loesung_stichprobenaufgabe.tex

geschrieben.
"""

from pathlib import Path

import numpy as np
from scipy import stats


ALPHA = 0.05


def latex_number(value, digits=4):
    """Formatiert eine Zahl fuer deutsche LaTeX-Schreibweise."""
    text = f"{value:.{digits}f}"
    return text.replace(".", "{,}")


def latex_integer_list(values):
    """Formatiert eine Liste von Messwerten fuer eine mathematische Umgebung."""
    return r",\,\,".join(str(int(value)) for value in values)


def calculate_results(sample1, sample2):
    """Berechnet alle statistischen Groessen der Aufgabe."""
    n1 = len(sample1)
    n2 = len(sample2)

    # Empirische Mittelwerte und Standardabweichungen.
    # ddof=1 bedeutet: Division durch n-1 statt durch n.
    mean1 = np.mean(sample1)
    mean2 = np.mean(sample2)
    s1 = np.std(sample1, ddof=1)
    s2 = np.std(sample2, ddof=1)
    variance1 = s1**2
    variance2 = s2**2

    # F-Test. Die Reihenfolge der Stichproben bleibt fest:
    # F = s1^2 / s2^2, mit nu1=n1-1 und nu2=n2-1.
    f_value = variance1 / variance2
    f_df1 = n1 - 1
    f_df2 = n2 - 1
    f_lower = stats.f.ppf(ALPHA / 2, f_df1, f_df2)
    f_upper = stats.f.ppf(1 - ALPHA / 2, f_df1, f_df2)
    f_p_value = 2 * min(
        stats.f.cdf(f_value, f_df1, f_df2),
        stats.f.sf(f_value, f_df1, f_df2),
    )
    equal_variances = f_lower <= f_value <= f_upper

    # Gepoolter t-Test, weil der F-Test die Annahme gleicher Varianzen
    # in dieser Lehrbuchaufgabe nicht verwirft.
    pooled_variance = (
        (n1 - 1) * variance1 + (n2 - 1) * variance2
    ) / (n1 + n2 - 2)
    pooled_standard_deviation = np.sqrt(pooled_variance)
    t_df = n1 + n2 - 2
    t_value = (mean1 - mean2) / (
        pooled_standard_deviation * np.sqrt(1 / n1 + 1 / n2)
    )
    t_critical = stats.t.ppf(1 - ALPHA / 2, t_df)
    t_p_value = 2 * stats.t.sf(abs(t_value), t_df)
    equal_means = abs(t_value) <= t_critical

    # Einseitiger Chi-Quadrat-Test fuer die Standardabweichung der
    # zweiten Stichprobe: H0: sigma=sigma0 gegen H1: sigma>sigma0.
    chi_df = n2 - 1
    chi_critical = stats.chi2.ppf(1 - ALPHA, chi_df)

    def variance_test(sigma0):
        chi_value = chi_df * variance2 / sigma0**2
        p_value = stats.chi2.sf(chi_value, chi_df)
        reject = chi_value > chi_critical
        return chi_value, p_value, reject

    chi_for_6 = variance_test(6.0)
    chi_for_9 = variance_test(9.0)

    return {
        "n1": n1,
        "n2": n2,
        "mean1": mean1,
        "mean2": mean2,
        "s1": s1,
        "s2": s2,
        "variance1": variance1,
        "variance2": variance2,
        "f_value": f_value,
        "f_df1": f_df1,
        "f_df2": f_df2,
        "f_lower": f_lower,
        "f_upper": f_upper,
        "f_p_value": f_p_value,
        "equal_variances": equal_variances,
        "pooled_variance": pooled_variance,
        "pooled_standard_deviation": pooled_standard_deviation,
        "t_value": t_value,
        "t_df": t_df,
        "t_critical": t_critical,
        "t_p_value": t_p_value,
        "equal_means": equal_means,
        "chi_df": chi_df,
        "chi_critical": chi_critical,
        "chi_for_6": chi_for_6,
        "chi_for_9": chi_for_9,
    }


def decision_text(reject):
    """Liefert eine deutsche Beschreibung der Testentscheidung."""
    if reject:
        return "Die Nullhypothese wird verworfen."
    return "Die Nullhypothese wird nicht verworfen."


def create_latex_solution(sample1, sample2, result):
    """Erzeugt den LaTeX-Quelltext der Musterloesung."""
    f = lambda key, digits=4: latex_number(result[key], digits)

    if result["equal_variances"]:
        f_decision = "Da der F-Wert im Annahmebereich liegt, werden die Varianzen als gleich angesehen."
    else:
        f_decision = "Da der F-Wert ausserhalb des Annahmebereichs liegt, werden die Varianzen als ungleich angesehen."

    if result["equal_means"]:
        t_decision = "Die Nullhypothese gleicher Erwartungswerte wird nicht verworfen."
    else:
        t_decision = "Die Nullhypothese gleicher Erwartungswerte wird verworfen."

    chi6_value, chi6_p, chi6_reject = result["chi_for_6"]
    chi9_value, chi9_p, chi9_reject = result["chi_for_9"]

    # Platzhalter (@@...@@) vermeiden Konflikte zwischen Python-format()
    # und den geschweiften Klammern, die LaTeX benoetigt.
    latex = r"""\documentclass[a4paper,11pt]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{amsmath,amssymb}
\usepackage{booktabs}
\usepackage[margin=25mm]{geometry}

\title{Lösung: Vergleich zweier Stichproben}
\author{}
\date{}

\begin{document}
\maketitle

Gegeben sind die Stichproben
\begin{align*}
X_1 &= (@@SAMPLE1@@),\\
X_2 &= (@@SAMPLE2@@).
\end{align*}
Das Signifikanzniveau beträgt in allen Tests
\[
\alpha=0{,}05.
\]

\begin{enumerate}


\item[a)] \textbf{Mittelwerte und Standardabweichungen}

Für die empirische Standardabweichung wird durch $J_i-1$ geteilt:
\[
s_i=\sqrt{\frac{1}{J_i-1}\sum_{j=1}^{J_i}
(x_{i,j}-\overline{x}_i)^2}.
\]
Damit erhält man
\[
\begin{array}{c|c|c|c}
\text{Stichprobe} & J_i & \overline{x}_i & s_i\\
\midrule
1 & @@N1@@ & @@MEAN1@@ & @@S1@@\\
2 & @@N2@@ & @@MEAN2@@ & @@S2@@
\end{array}
\]

\item[b)] \textbf{F-Test auf gleiche Varianzen}

Wir prüfen
\[
H_0:\ \sigma_1^2=\sigma_2^2
\qquad\text{gegen}\qquad
H_1:\ \sigma_1^2\ne\sigma_2^2.
\]
Die Prüfgröße wird hier als Quotient in der vorgegebenen Reihenfolge
der Stichproben gebildet:
\[
F=\frac{s_1^2}{s_2^2}.
\]
Unter $H_0$ gilt
\[
F\sim F_{\nu_1,\nu_2},
\qquad
\nu_1=J_1-1=@@FDF1@@,
\qquad
\nu_2=J_2-1=@@FDF2@@.
\]
Mit den empirischen Varianzen ergibt sich
\[
s_1^2=@@VAR1@@,
\qquad
s_2^2=@@VAR2@@,
\qquad
F=@@FVALUE@@.
\]
Der zweiseitige Annahmebereich lautet
\[
F_{\alpha/2;\,@@FDF1@@,@@FDF2@@}
\leq F \leq
F_{1-\alpha/2;\,@@FDF1@@,@@FDF2@@}.
\]
Numerisch ist dies
\[
@@FLOWER@@\leq F\leq@@FUPPER@@.
\]
@@FDECISION@@
Damit werden die Stichproben für den folgenden Test gepoolt.

\item[c)] \textbf{Test auf gleichen Erwartungswert}

Da der F-Test die Gleichheit der Varianzen nicht verworfen hat,
verwenden wir den gepoolten Zweistichproben-t-Test. Die Hypothesen lauten
\[
H_0:\ \mu_1=\mu_2
\qquad\text{gegen}\qquad
H_1:\ \mu_1\ne\mu_2.
\]
Die gepoolte Varianz ist
\[
s_p^2=
\frac{(J_1-1)s_1^2+(J_2-1)s_2^2}{J_1+J_2-2}
=@@POOLVAR@@.
\]
Die Prüfgröße lautet
\[
t=\frac{\overline{x}_1-\overline{x}_2}
{s_p\sqrt{1/J_1+1/J_2}}
=@@TVALUE@@.
\]
Unter $H_0$ ist diese Prüfgröße t-verteilt mit
\[
\nu=J_1+J_2-2=@@TDF@@
\]
Freiheitsgraden. Für den zweiseitigen Test wird das Quantil
\[
t_{1-\alpha/2;\,@@TDF@@}=@@TCRIT@@
\]
verwendet. Da
\[
|t|=@@TABS@@>@@TCRIT@@,
\]
gilt:
@@TDECISION@@

\item[d)] \textbf{Test der Standardabweichung mit $\sigma_0=6$}

Für die zweite Stichprobe prüfen wir einseitig
\[
H_0:\ \sigma=\sigma_0=6
\qquad\text{gegen}\qquad
H_1:\ \sigma>6.
\]
Bei normalverteilten Beobachtungen gilt unter $H_0$
\[
Q=\frac{(J_2-1)s_2^2}{\sigma_0^2}
\sim\chi^2_{J_2-1}.
\]
Es handelt sich daher um einen Chi-Quadrat-Test. Für den rechten Rand
ist der kritische Wert
\[
\chi^2_{1-\alpha;\,@@CHIDF@@}=@@CHICRIT@@.
\]
Die Prüfgröße beträgt
\[
Q=\frac{@@CHIDF@@\cdot@@VAR2@@}{6^2}=@@CHI6@@.
\]
Da $Q>\chi^2_{1-\alpha;\,@@CHIDF@@}$, folgt:
@@CHI6DECISION@@
Der zugehörige einseitige p-Wert beträgt ungefähr
\[
p=@@CHI6P@@.
\]

\item[e)] \textbf{Test der Standardabweichung mit $\sigma_0=9$}

Nun wird geprüft:
\[
H_0:\ \sigma=9
\qquad\text{gegen}\qquad
H_1:\ \sigma>9.
\]
Die Prüfgröße ist
\[
Q=\frac{@@CHIDF@@\cdot@@VAR2@@}{9^2}=@@CHI9@@.
\]
Der kritische Wert bleibt unverändert:
\[
\chi^2_{1-\alpha;\,@@CHIDF@@}=@@CHICRIT@@.
\]
Damit gilt:
@@CHI9DECISION@@
Der einseitige p-Wert beträgt ungefähr
\[
p=@@CHI9P@@.
\]

\end{enumerate}

\end{document}
"""

    replacements = {
        "@@SAMPLE1@@": latex_integer_list(sample1),
        "@@SAMPLE2@@": latex_integer_list(sample2),
        "@@N1@@": str(result["n1"]),
        "@@N2@@": str(result["n2"]),
        "@@FDF1@@": str(result["f_df1"]),
        "@@FDF2@@": str(result["f_df2"]),
        "@@FVALUE@@": f("f_value"),
        "@@FLOWER@@": f("f_lower"),
        "@@FUPPER@@": f("f_upper"),
        "@@VAR1@@": f("variance1"),
        "@@VAR2@@": f("variance2"),
        "@@MEAN1@@": f("mean1"),
        "@@MEAN2@@": f("mean2"),
        "@@S1@@": f("s1"),
        "@@S2@@": f("s2"),
        "@@POOLVAR@@": f("pooled_variance"),
        "@@TVALUE@@": f("t_value"),
        "@@TABS@@": latex_number(abs(result["t_value"])),
        "@@TDF@@": str(result["t_df"]),
        "@@TCRIT@@": f("t_critical"),
        "@@CHIDF@@": str(result["chi_df"]),
        "@@CHICRIT@@": f("chi_critical"),
        "@@CHI6@@": latex_number(chi6_value),
        "@@CHI6P@@": latex_number(chi6_p),
        "@@CHI9@@": latex_number(chi9_value),
        "@@CHI9P@@": latex_number(chi9_p),
        "@@FDECISION@@": f_decision,
        "@@TDECISION@@": t_decision,
        "@@CHI6DECISION@@": decision_text(chi6_reject),
        "@@CHI9DECISION@@": decision_text(chi9_reject),
    }

    for placeholder, replacement in replacements.items():
        latex = latex.replace(placeholder, replacement)

    return latex


def main():
    sample1 = np.array([53, 69, 63, 47, 49, 44, 47, 44, 38], dtype=float)
    sample2 = np.array([21, 33, 19, 39, 7], dtype=float)

    result = calculate_results(sample1, sample2)
    latex_solution = create_latex_solution(sample1, sample2, result)

    output_path = Path("loesung_stichprobenaufgabe.tex")
    output_path.write_text(latex_solution, encoding="utf-8")

    print("LaTeX-Loesung geschrieben nach:", output_path)
    print("F-Wert:", f"{result['f_value']:.4f}")
    print("t-Wert:", f"{result['t_value']:.4f}")
    print("Chi-Quadrat-Wert fuer sigma0=6:", f"{result['chi_for_6'][0]:.4f}")
    print("Chi-Quadrat-Wert fuer sigma0=9:", f"{result['chi_for_9'][0]:.4f}")


if __name__ == "__main__":
    main()
