# AGGP #

#### Générateur de graphes biologiquement plausibles par algorithme génétique.

### Fichiers & Dossiers ###


Le projet est divisé en plusieurs dossiers :
* src

	Sources du programme

        main.py                 lancement du programme
        population.py           algorithme génétique (évaluation, sélection, mutations/reproduction)
        individuals.py          calcul du score de l'individu
        mutations.py            fonctions de mutations
        globals.py              paramètres et fonctions globales
        progressbar.py          thread d'affichage barre d'avancement
* in

	Fichiers d'entrée nécessaires au programme 

        names			Liste de prénoms pour identifier les individus
        scores.gp		Fichier gnuplot (-t ou --no-stat pour ne pas lancer gnuplot)
* ref_network

	Réseaux de références (non nécessaires au programme)
* img

	Dossier créé si inexistant au lancement, contiendra les visualisations des individus 
* out

	Dossier créé si inexistant au lancement, contiendra des fichiers .sif pour une visualisation Cytoscape et des fichiers evo_* pour une visualisation par type de scores au cours des génerations (gnuplot).

### Dépendances ###
* librairies de base de `Python v2.7`
* numpy
* networkx
* graphviz (optionnelle mais conseillée)
* matplotlib
* gnuplot

### Lancement du programme ###
Une aide est disponible, expliquant les différents paramètres optionels possibles et l'usage général du programme :

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
          -t, --stat            Plot final stats graphs
          -d, --delete          Delete all output (files, graphs and pictures) from
                                previous run


L'option `-v` peut etre multipliée pour graduer son effet. Sont alors acceptés `-vv` ou `-vvv`. Ces options sont exclusives avec `-q`.

### Crédits ###
Jonas Abernot, Arthur Bailly, Johan Chan & Maxime Sainlot

4BIM - INSA Lyon - AGGP
