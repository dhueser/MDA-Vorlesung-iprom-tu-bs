#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Veranschauliche die Linearisierung von Y = X**2 bei x = 0.

Voraussetzungen: Python 3, NumPy und Matplotlib (kein LaTeX erforderlich).

Aufruf:
    python verschwindende_sensitivitaet.py
    python verschwindende_sensitivitaet.py --sigma 2 --show

Das Diagramm zeigt die Modellfunktion und ihre horizontale Tangente,
die normalverteilte Eingangsgröße und ein Monte-Carlo-Histogramm der
Ausgangsgröße. Das Histogramm stellt die ursprüngliche Modellfunktion
dar, nicht ihre Linearisierung. Alle linearisierten Ausgangswerte sind null.

Exakte Kenngrößen:
    E(X) = 0, Var(X) = sigma**2
    E(Y) = sigma**2, Var(Y) = 2 * sigma**4
    u(Y) = sqrt(2) * sigma**2

Die sichtbaren Achsen umfassen X in [-3*sigma, 3*sigma] und
Y bis 9*sigma**2. Das Histogramm wird mit der Gesamtzahl aller Ziehungen
normiert; die wenigen Werte außerhalb des sichtbaren Bereichs werden
nicht in diesen Bereich umverteilt.

PNG, PDF und SVG werden standardmäßig neben dem Skript gespeichert.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def create_figure(sigma=1.0, samples=200_000, seed=7):
    """Erzeuge die Abbildung und gib die Monte-Carlo-Kenngrößen zurück."""
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("sigma muss endlich und größer als null sein.")
    if samples < 1:
        raise ValueError("samples muss mindestens 1 sein.")

    rng = np.random.default_rng(seed)
    draws_x = rng.normal(loc=0.0, scale=sigma, size=samples)
    draws_y = draws_x**2
    x = np.linspace(-3 * sigma, 3 * sigma, 1201)
    y = x**2
    density_x = np.exp(-0.5 * (x / sigma)**2) / (sigma * np.sqrt(2 * np.pi))

    blue, orange, green = "#245ba8", "#d97706", "#15803d"
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
    fig = plt.figure(figsize=(11.5, 7.2), layout="constrained")
    grid = fig.add_gridspec(2, 2, width_ratios=[3.0, 1.25], height_ratios=[3.2, 1.05])
    ax_model = fig.add_subplot(grid[0, 0])
    ax_output = fig.add_subplot(grid[0, 1], sharey=ax_model)
    ax_input = fig.add_subplot(grid[1, 0], sharex=ax_model)
    ax_info = fig.add_subplot(grid[1, 1])

    fig.suptitle(
        r"Verschwindende Sensitivität: $Y=X^2$, $X\sim\mathcal{N}(0,\sigma^2)$"
        + rf", $\sigma={sigma:g}$",
        fontsize=15,
    )

    # Modell und Linearisierung besitzen am Ursprung denselben Wert.
    ax_model.plot(x, y, color=blue, lw=2.6, label=r"Modell: $f(X)=X^2$")
    ax_model.axhline(0, color=orange, lw=2.3, label=r"Tangente: $Y_{\mathrm{lin}}=0$")
    ax_model.axhline(sigma**2, color=green, lw=1.6, ls="--",
                     label=r"Erwartungswert: $\mathrm{E}(Y)=\sigma^2$")
    ax_model.axvline(0, color="0.65", lw=0.9, ls=":")
    ax_model.scatter([0], [0], s=55, color="0.2", zorder=6)
    ax_model.annotate(
        "Entwicklungspunkt\n" + r"$x=0$, $f(0)=0$",
        xy=(0, 0), xytext=(0.25 * sigma, 0.22 * sigma**2),
        fontsize=9, va="bottom",
        arrowprops={"arrowstyle": "->", "color": "0.35", "lw": 0.8},
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9},
    )

    # Zwei symmetrische Eingangsabweichungen: echte Modellwerte vs. Tangente.
    for value_x in [-1.6 * sigma, 1.6 * sigma]:
        ax_model.plot([value_x, value_x], [0, value_x**2], color=blue, ls=":", lw=1.2)
        ax_model.scatter([value_x], [value_x**2], color=blue, s=43, zorder=5)
        ax_model.scatter([value_x], [0], facecolors="white", edgecolors=orange,
                         linewidths=1.8, s=48, zorder=5)
    ax_model.annotate(
        r"$X=\pm1{,}6\sigma$" + "\n" + r"$Y=2{,}56\sigma^2$",
        xy=(1.6 * sigma, (1.6 * sigma)**2), xytext=(0.65 * sigma, 4.2 * sigma**2),
        fontsize=10, color=blue,
        arrowprops={"arrowstyle": "->", "color": blue, "lw": 0.9},
    )
    ax_model.set(xlim=(-3 * sigma, 3 * sigma), ylim=(-0.5 * sigma**2, 9.5 * sigma**2),
                 ylabel=r"Ausgangswert $Y$", title="Modellfunktion und lineare Näherung")
    ax_model.tick_params(labelbottom=False)
    ax_model.grid(alpha=0.18)
    ax_model.legend(loc="upper center", fontsize=9, framealpha=0.95)

    # Normaldichte der Eingangsgröße; farbiger Bereich: +/- eine Standardabweichung.
    ax_input.plot(x, density_x, color=blue, lw=2)
    ax_input.fill_between(x, density_x, color=blue, alpha=0.10)
    within_sigma = np.abs(x) <= sigma
    ax_input.fill_between(x, density_x, where=within_sigma, color=blue, alpha=0.25)
    ax_input.axvline(0, color="0.5", ls=":", lw=1)
    ax_input.set(xlabel=r"Eingangswert $X$", ylabel=r"Dichte $p_X(X)$",
                 ylim=(0, 0.46 / sigma), title="Eingangsverteilung: Streuung um null")
    ax_input.grid(alpha=0.15)

    # Horizontales Histogramm mit denselben Y-Koordinaten wie die Modellfunktion.
    bins = np.linspace(0, 9 * sigma**2, 76)
    counts, edges = np.histogram(draws_y, bins=bins)
    heights = np.diff(edges)
    density_y = counts / (samples * heights)
    ax_output.barh(edges[:-1], density_y, height=heights, align="edge",
                   color=blue, alpha=0.60, edgecolor="white", linewidth=0.3)
    density_max = float(density_y.max())
    ax_output.axhline(sigma**2, color=green, lw=1.6, ls="--")
    ax_output.axhline(0, color=orange, lw=2)
    ax_output.scatter([0], [0], facecolors="white", edgecolors=orange,
                      linewidths=1.8, s=55, zorder=5, clip_on=False)
    ax_output.annotate(
        "Linearisierung:\nalle Ausgangswerte\nbei " + r"$Y=0$",
        xy=(0, 0), xytext=(0.33 * density_max, 3.5 * sigma**2),
        fontsize=9, color=orange,
        arrowprops={"arrowstyle": "->", "color": orange, "lw": 1.1},
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.95},
    )
    ax_output.set(xlim=(0, 1.08 * density_max), xlabel=r"Dichte $p_Y(Y)$",
                  title="Ausgangsverteilung\n(Monte-Carlo-Histogramm)")
    ax_output.tick_params(labelleft=False)
    ax_output.grid(axis="y", alpha=0.15)

    # Exakte Momente: das Histogramm dient allein der Veranschaulichung.
    ax_info.axis("off")
    ax_info.text(0, 1, "Exakte Kenngrößen", fontsize=11, weight="bold", va="top")
    ax_info.text(0, 0.80,
                 r"$\mathrm{E}(Y)=\sigma^2$" + "\n" +
                 r"$u_{\mathrm{c}}(y)=\sqrt{2}\,\sigma^2$" + "\n\n" +
                 r"$y_{\mathrm{lin}}=0$, $u_{\mathrm{c}}(y_{\mathrm{lin}})=0$",
                 fontsize=11, va="top", linespacing=1.6)

    statistics = {
        "mean_mc": float(draws_y.mean()),
        "standard_deviation_mc": float(draws_y.std(ddof=0)),
        "mean_exact": sigma**2,
        "standard_deviation_exact": np.sqrt(2) * sigma**2,
        "visible_probability_mc": float(counts.sum() / samples),
    }
    return fig, statistics


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sigma", type=float, default=1.0, help="Standardabweichung von X (Standard: 1)")
    parser.add_argument("--samples", type=int, default=200_000, help="Anzahl der Monte-Carlo-Ziehungen")
    parser.add_argument("--seed", type=int, default=7, help="Zufallsstartwert für reproduzierbare Ergebnisse")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent,
                        help="Verzeichnis für die Abbildungen")
    parser.add_argument("--show", action="store_true", help="Diagramm zusätzlich in einem Fenster anzeigen")
    args = parser.parse_args()

    fig, stats = create_figure(args.sigma, args.samples, args.seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for extension in ["png", "pdf", "svg"]:
        output_path = args.output_dir / ("verschwindende_sensitivitaet." + extension)
        fig.savefig(output_path, dpi=220, facecolor="white")
        print("Gespeichert:", output_path)
    print("Erwartungswert: exakt = {mean_exact:.6g}, Monte Carlo = {mean_mc:.6g}".format(**stats))
    print("Standardunsicherheit: exakt = {standard_deviation_exact:.6g}, "
          "Monte Carlo = {standard_deviation_mc:.6g}".format(**stats))
    print("Wahrscheinlichkeitsanteil im dargestellten Y-Bereich: "
          "{visible_probability_mc:.3%}".format(**stats))
    if args.show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
