# AGGP #

#### Générateur de graphes biologiquement plausibles par algorithme génétique.

### Fichiers ###
Le projet est divisé en plusieurs fichiers :
* main.py		*lancement du programme*
* population.py		*algorithme génétique*
* individuals.py	*calcul du score*
* mutations.py		*fonctions de mutations*
* globals.py		*paramètres et fonctions globales*
* progressbar.py	*thread d'affichage barre d'avancement*

### Dépendances ###
Les librairies de base de `Python v2.7` + :
* numpy
* networkx
* graphviz (optionnelle)
* matplotlib

### Lancement du programme ###
	> python main.py -h
	usage: main.py [options]

	Biological Graph Generator

	optional arguments:
	  -h, --help            show this help message and exit
	  -p, --param           Ask for every parameters of the simulation
	  -v, --verbose
	  -q, --quiet
	  --no-progress         Disable the progress bar
	  -f FREQ, --freq FREQ  Frequency of displaying informations
	  -s, --no-save         Do not save and plot individuals
	  -d, --delete          Delete all output (files, graphs and pictures) from
		                previous run

### Crédits ###
Jonas Abernot, Arthur Bailly, Johan Chan & Maxime Sainlot

4BIM - INSA Lyon - AGGP
