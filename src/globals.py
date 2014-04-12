#! /usr/bin/python
# -*- coding: utf-8 -*-

import networkx as nx
import sys, os

####################################################################
#			Global Parameters
####################################################################
# Paths
IMG = "../img/"
IN = "../in/"
OUT = "../out/"

for PATH in [IMG,OUT]:
    if not os.path.exists(PATH):
        os.makedirs(PATH)
        
# Plot png in /img
PLOT_PDL = 1001 # plot degree graph every X generation
#PLOT_CF = 1001 # plot clique formation graph every X generation
PLOT_GR = 1001 # plot graph every X generation
PLOT_GEN_ZERO = False # plot initials individuals ?

# Parameters of algogen
NB_GEN = 1000 # genetic algo's iteration number
NB_NODES = 60
NB_INDIV = 30

# Plot information box
INFO_INDIV = False
INFO_BEST = True
INFO_SELECT = False
INFO_GEN = True
INFO_FREQ = 10 # information frequency (every X generation)
PROGRESS_GEN = False if INFO_INDIV or NB_NODES*NB_INDIV < 1000 or PLOT_PDL < NB_GEN or PLOT_GR < NB_GEN else True

# Rates
RATE_ELITISM = 0.2 
RATE_TOURNAMENT = 1-RATE_ELITISM
RATE_CROSS = 0.6
RATE_MUT = 0.3

# Scores Rates
PDL = 1
SW = 1
CF = 1

# Random Reference
#G_RAND = nx.fast_gnp_random_graph(NB_NODES,0.2)
C_RAND = 1 #nx.average_clustering(G_RAND)
L_RAND = 1 #nx.average_shortest_path_length(G_RAND)

# Miscellaneous
NAMES = open(IN+"names").read().splitlines()
ERROR = False
EPS = 0.001 # log(x+EPS) to avoid log(0)

# Colors
HEADER = '\033[1m' # bold
OKBLUE = '\033[94m' # blue
OKGREEN = '\033[92m' # green
WARNING = '\033[93m' # yellow
FAIL = '\033[91m' # red
ENDC = '\033[0m' # back to normal


####################################################################
#			Global Functions
####################################################################
def symetrize(a):
    """ deprecated """
    return -a*a.T + a + a.T

def update_progress(label,progress,bar_length=25): # small 20, medium 25, large 50
    progress = int(progress)
    if progress > 100 : progress = 100
    sys.stdout.write('\r{2:<25} [{0}] {1:3d}%'.format('#'*(progress/int(100./bar_length))+'-'*(bar_length-(progress/int(100./bar_length))), progress,label))
    sys.stdout.flush()

def weighted_sample(items, n):
    """ deprecated (used in wheel selection)"""
    total = float(sum(w for w, v in items))
    i = 0
    w, v = items[0]
    while n:
        x = total * (1 - rd.random() ** (1.0 / n))
        total -= x
        while x > w:
            x -= w
            i += 1
            w, v = items[i]
            w -= x
            yield v
            n -= 1
