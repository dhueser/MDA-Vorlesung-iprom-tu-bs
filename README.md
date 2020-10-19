# Vorlesung _Messdatenauswertung und Messunsicherheit am iprom_ der TU Braunschweig
## Einleitung
Diese Vorlesung ist ein Modul des Studiengangs _Messtechnik und Analytik_ im Fachbereich Maschinenbau der TU Braunschweig.

Der Modulverantwortliche ist Prof. Dr.-Ing Rainer Tutsch, Leiter des Instituts für Produktionsmesstechnik. Die Vorlesung wird gehalten von Mitarbeitern der Physikalisch-Technischen Bundesanstalt, die das nationale Metrologie-Institut zur Wahrung der SI-Einheiten und verantworlich für das gesetzliche Messwesen in Deutschland ist:
- Dr. habil. Dorothee Hüser (Hauptverantwortliche)
- Dr.-Ing. Gerd Ehret
- Dr. rer. nat. Wolfgang Schmid

## [Vorlesungsskript](https://github.com/dhueser/MDA-Vorlesung-iprom-tu-bs/releases/latest/download/MDA_iprom_TUbraunschweig.pdf)
- Das pdf-File des jeweils aktuellen Vorlesungsskriptes steht zum Download [hier](https://github.com/dhueser/MDA-Vorlesung-iprom-tu-bs/releases/latest/download/MDA_iprom_TUbraunschweig.pdf) zur Verfügung.
- Das Skript wird laufend überarbeitet.
- Wer Vorschläge zur Verbesserung hat, ist herzlich eingeladen, die LaTeXquellen des gesamten Skripts zu clonen und Änderungen in einem separaten Branch in das github-Repository zu laden und einen pull Request zu setzen. Tutorial: (https://help.github.com/en/articles/setting-guidelines-for-repository-contributors)[Tutorial: https://help.github.com/en/articles/setting-guidelines-for-repository-contributors]; Beispiel: (https://github.com/tensorflow/tensorflow/blob/master/CONTRIBUTING.md)[https://github.com/tensorflow/tensorflow/blob/master/CONTRIBUTING.md]

## Inhalt
Das Modul gibt einen umfangreichen Überblick über die Grundlagen der Messdatenauswertung und der rechnerischen Ermittlung der zugeordneten Unsicherheiten.

### Einführung: Messtechnik und Statistik
### Lineare Regression (Ausgleichsrechnung)
### Konzepte der Statistik für die Messdatenanalyse
### Optimierung nichtlinearer Modelle
### Wahrscheinlichkeiten und Hypothesentests
### Auswertung von Mess- und Ringvergleichen
### Messunsicherheitsfortpflanzung linearer und linearisierbarer Modelle
### Messunsicherheitsfortpflanzung mit Bayesstatistik
### Monte-Carlo-Methoden
### Infrastrukturen der Metrologie
### Normen und Richtlinien                             

## WS 2019/20

| Nr.  | Datum      | Dozent                                  | Thema                                               |
| ---- | ---------- | --------------------------------------- | --------------------------------------------------- |
| 1    | 2020-10-19 | Dr. habil. Dorothee Hüser               | Einführung: Messtechnik und Statistik               |
| 2    | 2020-10-26 | Dr.-Ing. Gerd Ehret                     | Lineare Regression                                  |
| 3    | 2020-11-02 | Dr.-Ing. Gerd Ehret                     | Optimierung nichtlinearer Modelle                   |
| 4    | 2020-11-09 | Dr. habil. Dorothee Hüser               | Konzepte der Statistik für die Messdatenanalyse     |
| 5    | 2020-11-16 | Dr. habil. Dorothee Hüser               | Wahrscheinlichkeiten und Hypothesentests            |
| 6    | 2020-11-23 | Dr.-Ing. Gerd Ehret                     | Auswertung von Mess- und Ringvergleichen            |
| 7    | 2020-11-30 | Dr. habil. Dorothee Hüser               | Messunsicherheitsfortpflanzung bei linearen Modelle |
| 8    | 2020-12-07 | Dr.-Ing. Gerd Ehret                     | Bayes-Statistik                                     |
| 9    | 2020-12-14 | Dr. Wolfgang Schmid                     | Beispiele zur Messunsicherheitsfortpflanzung nach GUM |
| 10   | 2020-01-11 | Dr.-Ing. Gerd Ehret                     | 1. Konformitätsbewertung (GUM - JCGM 106:2012), 2. Monte-Carlo Methoden I |
| 11   | 2020-01-18 | Dr. Wolfgang Schmid                     | Monte Carlo Methoden II, Beispiele zu Monte-Carlo-Methoden nach GUM-S1  |
| 12   | 2020-01-25 | Dr.-Ing. Gerd Ehret/Dr. Wolfgang Schmid | Infrastrukturen der Metrologie und ISO 17025: (1.) Internationale: BIPM, EURAMET; (2.) nationale: DAkkS, Kalibrierlabore |
| 13   | 2020-02-03 | Dr. Rolf Kumme                          | Normen und Richtlinien                              | 14   | nach Vereinb.| Dr. Rolf Kumme | Exkursion |

# Interaktive Code-Beispiele

Die interaktiven Code-Beispiele sind in folgender Tabelle zusammengefasst. Diese können mit einem Klick auf den Link geöffnet werden. Der Code wird auf einem Server ausgeführt und die Anzeige findet im Browser statt. Es ist also nicht notwendig zusätzliche Programme zu installieren. Es kann mitunter ein paar Minuten dauern bis die passende Programmierumgebung erstellt wird.

| Kapitel | Titel         | Beschreibung                                                 | Link                                                         |
| ------- | ------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| 2       | Kubischer Fit | Matlab/Octave-Code, um einen kubischen Fit zu bestimmen mit anschließender Visualisierung in einem Plot. | [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/dhueser/MDA-Vorlesung-iprom-tu-bs/master?urlpath=/lab/tree/vorlesung/02_vorlesung/code/cubic_fit.ipynb) |
| 5       | Hypothesentests | Python-Code für Hypothesentests. | [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/dhueser/MDA-Vorlesung-iprom-tu-bs/master?urlpath=/lab/tree/vorlesung/05_vorlesung/code/hypothesentests.ipynb) |
