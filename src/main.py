#! /usr/bin/python
# -*- coding: utf-8 -*-
# Dependancies : networkx, numpy, graphviz, matplotlib

## python main.py -h
## usage: main.py [options]

## Biological Graph Generator

## optional arguments:
##   -h, --help            show this help message and exit
##   -p, --param           Ask for every parameters of the simulation
##   -v, --verbose
##   -q, --quiet
##   --no-progress         Disable the progress bar
##   -f FREQ, --freq FREQ  Frequency of displaying informations
##   -s, --no-save         Do not save and plot individuals
##   -t, --stat            Plot final stats graphs
##   -d, --delete          Delete all output (files, graphs and pictures) from
##                         previous run

                                       
from population import *

####################################################################
#			Main
####################################################################
if __name__ == "__main__":
    a = Population()
    a.genetic_algo()
