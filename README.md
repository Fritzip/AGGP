# AGGP #

#### Générateur de graphes biologiquement plausibles par algorithme génétique.

### Fichiers ###
Le projet est divisé en plusieurs fichiers :
	main.py			*lancement du programme*
	population.py		*algorithme génétique (évaluation, sélection, mutations/reproduction)*
	individuals.py		*calcul du score de l'individu*
	mutations.py		*fonctions de mutations*
	globals.py		*paramètres et fonctions globales*
	progressbar.py		*thread d'affichage barre d'avancement*

### Dépendances ###
Les librairies de base de `Python v2.7` + :
* numpy
* networkx
* graphviz (optionnelle mais conseillée)
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
	  -d, --delete          Delete all output (files, graphs and pictures) from previous run

L'option `-v` peut etre multipliée pour graduer son effet. Sont alors accepté `-vv` ou `-vvv`. Cette option est exclusive avec `-q`.
### Fichiers de Sortie ###
Suivant les paramètres utilisés, le programme peut créer dans le dossier `out/` les réseaux biologiques, la liste d'adjacence, l'évolution des scores pour la totalité ou une proportion des meilleurs individus obtenus en sortie de simulation.

### Crédits ###
Jonas Abernot, Arthur Bailly, Johan Chan & Maxime Sainlot

4BIM - INSA Lyon - AGGP
