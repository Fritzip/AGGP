#! /usr/bin/python
# -*- coding: utf-8 -*-
# Dependancies : networkx, numpy, graphviz, matplotlib

import networkx as nx
import sys, os, shutil
import argparse


####################################################################
#			Global Parameters (default)
####################################################################
# Parameters of algogen
NB_GEN = 100 # genetic algo's iteration number
NB_NODES = 60
NB_INDIV = 30

# Plot png in /img
PLOT_PDL = NB_GEN+1 # plot degree graph every X generation
PLOT_GR = NB_GEN+1 # plot graph every X generation
PLOT_CF = NB_GEN+1
CYTO3D = False
SAVE = True
PLOT_GEN_ZERO = False # plot initials individuals ?

# Plot information box
INFO_INDIV = False
INFO_BEST = True
INFO_SELECT = False
INFO_GEN = True
INFO_FREQ = 5 # information frequency (every X generation)
PROGRESS_GEN = False if INFO_INDIV or NB_NODES*NB_INDIV < 1000 or PLOT_PDL < NB_GEN or PLOT_GR < NB_GEN else True

# Rates
RATE_ELITISM = 0.2 
RATE_TOURNAMENT = 1-RATE_ELITISM
RATE_CROSS = 0.4
# Mutation rates
RATE_SUB = 0.15
RATE_LOCAL_DEL = 0.1
RATE_LOCAL_INS = 0.1
RATE_GLOBAL_INS = 0.05
RATE_GLOBAL_DEL = 0.05

# Scores Rates
PDL = 10
SW = 1
CF = 100



####################################################################
#			Arguments parser
####################################################################
parser = argparse.ArgumentParser(description="Biological Graph Generator",usage='%(prog)s [options]')
group = parser.add_mutually_exclusive_group()

#parser.add_argument("-i", metavar="FILE",
#                    help="Take file of parameters as input")

#parser.add_argument("-p","--param",action="store_true",
#                    help="Ask for every parameters of the simulation")

group.add_argument("-v", "--verbose", action="count", default=0)

group.add_argument("-q", "--quiet", action="store_true", default=0)

parser.add_argument("--no-progress", action="store_true", default=0,
                    help="Disable the progress bar")

parser.add_argument("-f","--freq", default=INFO_FREQ, type=int,
                    help="Frequency of displaying informations")

#parser.add_argument("-g","--graph", metavar="X",
#                    help="Plot graph output every X generation")

#parser.add_argument("-s","--save",
#                    help="Ask at the end if and how you to save and plot individuals")

parser.add_argument("-d","--delete",action="store_true",
                    help="Delete all output (files, graphs and pictures) from previous run")


args = parser.parse_args()

####################################################################
#		Arguments parser consequences and other initialization stuff
####################################################################

INFO_FREQ=args.freq

if args.quiet:
    INFO_GEN = False
    INFO_BEST = False
    INFO_SELECT = False
    INFO_INDIV = False
    PROGRESS_GEN = False
    QUIET = True
    VERBOSE = False
else:
    VERBOSE = True
    QUIET = False
    
if args.verbose >= 4:
    INFO_GEN = True
    INFO_BEST = True
    INFO_SELECT = True
    INFO_INDIV = True
    PROGRESS_GEN = False
elif args.verbose >= 3:
    INFO_GEN = True
    INFO_BEST = True
    INFO_SELECT = True
    INFO_INDIV = False
    PROGRESS_GEN = True
elif args.verbose >= 2:
    INFO_GEN = True
    INFO_BEST = True
    INFO_SELECT = False
    INFO_INDIV = False
    PROGRESS_GEN = True
elif args.verbose >= 1: 
    INFO_GEN = True
    INFO_BEST = False
    INFO_SELECT = False
    INFO_INDIV = False
    PROGRESS_GEN = True

if args.no_progress:
    PROGRESS_GEN=False

# Paths
IMG = "../img/"
IN = "../in/"
OUT = "../out/"

if args.delete:
    if VERBOSE:
        print "Removing previous outputs"
    try:
        shutil.rmtree(IMG)
        shutil.rmtree(OUT)
    except:
        pass
    
for PATH in [IMG,OUT]:
    if not os.path.exists(PATH):
        os.makedirs(PATH)
        

# Random Reference
short_path_list=[]
coeff_clustering_list=[]
for i in range(NB_INDIV) :
    short_path_list.append(nx.average_shortest_path_length(nx.fast_gnp_random_graph(NB_NODES,0.2)))
    coeff_clustering_list.append(nx.average_clustering(nx.fast_gnp_random_graph(NB_NODES,0.2)))
C_RAND = sum(coeff_clustering_list)/len(coeff_clustering_list) 
L_RAND = sum(short_path_list)/len(short_path_list)

# Miscellaneous
NAMES = open(IN+"names").read().splitlines()
ERROR = False
EPS = 0.0001 # log(x+EPS) to avoid log(0)

# Colors
HEADER = '\033[1m' # bold
OKBLUE = '\033[94m' # blue
OKGREEN = '\033[92m' # green
WARNING = '\033[93m' # yellow
FAIL = '\033[91m' # red
UNDERLINE = '\033[4m'
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

def update_gen(gen):
    sys.stdout.write('\r{0:<25}'.format("Generation {0:3d}/{1}".format(gen,NB_GEN)))
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
