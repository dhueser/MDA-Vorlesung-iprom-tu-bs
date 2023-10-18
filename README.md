# Vorlesung _Messdatenauswertung und Messunsicherheit am iprom_ der TU Braunschweig
## Einleitung
Diese Vorlesung ist ein Modul des Studiengangs _Metrologie und Messtechnik_ im Fachbereich Maschinenbau der TU Braunschweig.

Der Modulverantwortliche ist Prof. Dr.-Ing Rainer Tutsch, Leiter des Instituts für Produktionsmesstechnik. Die Vorlesung wird gehalten von Mitarbeitern der Physikalisch-Technischen Bundesanstalt, die das nationale Metrologie-Institut zur Wahrung der SI-Einheiten und verantworlich für das gesetzliche Messwesen in Deutschland ist:
- Dr. habil. Dorothee Hüser (Hauptverantwortliche)
- Dr.-Ing. Gerd Ehret

## [Vorlesungsskript](https://github.com/dhueser/MDA-Vorlesung-iprom-tu-bs/releases)
- Das pdf-File des jeweils aktuellen Vorlesungsskriptes liegt innerhalb des zip bzw. des tar.gz-Archivs.
- Das Skript wird laufend überarbeitet, zuletzt wurde es am 5. Dez 2021 ge&auml;ndert und liegt unter dem Datenamen `MDA_iprom_TUbraunschweig_v2021-12-05.pdf` hier im Hauptverzeichnis [MDA_iprom_TUbraunschweig_v2021-12-05.pdf](MDA_iprom_TUbraunschweig_v2021-12-05.pdf). Ab dem WS 2023/24 wird das Skript neu gegliedert und die Nummerierung der Kapitel mit den Vorlesungen in Einklang gebracht.
- Wer Vorschläge zur Verbesserung hat, ist herzlich eingeladen, die LaTeXquellen des gesamten Skripts zu clonen und Änderungen in einem separaten Branch in das github-Repository zu laden und einen pull Request zu setzen. [Tutorial](https://help.github.com/en/articles/setting-guidelines-for-repository-contributors); [Beispiel](https://github.com/tensorflow/tensorflow/blob/master/CONTRIBUTING.md)

## Inhalt
Das Modul gibt einen umfangreichen Überblick über die Grundlagen der Messdatenauswertung und der rechnerischen Ermittlung der zugeordneten Unsicherheiten.

1. Einführung: Messtechnik und Statistik
2. Lineare Regression (Ausgleichsrechnung)
3. Konzepte der Statistik für die Messdatenanalyse
4. Optimierung nichtlinearer Modelle
5. Wahrscheinlichkeiten und Hypothesentests
6. Konformitätesbewertung
7. Auswertung von Mess- und Ringvergleichen
8. Infrastrukturen der Metrologie
9. Konzepte und Methoden der Messunsicherheitsfortpflanzung
10. Messunsicherheitsfortpflanzung linearisierbarer Modelle
11. Monte-Carlo-Methoden für die Messunsicherheitsfortpflanzung
12. Messunsicherheitsfortpflanzung mit Bayesstatistik                      

# Interaktive Code-Beispiele

Die interaktiven Code-Beispiele sind in folgender Tabelle zusammengefasst. Diese können mit einem Klick auf den Link geöffnet werden. Der Code wird auf einem Server ausgeführt und die Anzeige findet im Browser statt. Es ist also nicht notwendig zusätzliche Programme zu installieren. Es kann mitunter ein paar Minuten dauern bis die passende Programmierumgebung erstellt wird.

| Kapitel | Titel         | Beschreibung                                                 | Link                                                         |
| ------- | ------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| 2       | Kubischer Fit | Matlab/Octave-Code, um einen kubischen Fit zu bestimmen mit anschließender Visualisierung in einem Plot. | [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/dhueser/MDA-Vorlesung-iprom-tu-bs/master?urlpath=/lab/tree/vorlesung/02_vorlesung/code/cubic_fit.ipynb) |
| 5       | Hypothesentests | Python-Code für Hypothesentests. | [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/dhueser/MDA-Vorlesung-iprom-tu-bs/master?urlpath=/lab/tree/vorlesung/05_vorlesung/code/hypothesentests.ipynb) |
