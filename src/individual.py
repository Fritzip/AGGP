#! /usr/bin/python
# -*- coding: utf-8 -*-
# Dependancies : networkx, numpy, graphviz, matplotlib

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import random as rd
import copy, math, sys, time, subprocess, os
from scipy import stats as stats

from globals import *
from mutations import *

        
####################################################################
#			Individual
####################################################################

class Individual():
    """ Metabolic graph """
    def __init__(self,mat=0, nb_nodes=NB_NODES,id=rd.choice(NAMES)):
        self.id = id
        self.score_pdl = 0
        self.score_sw = 0
        self.score_cf = 0
        self.score = 0
        self.penalite = 0
        self.go_on = True

        self.list_degrees=[]
        self.list_count=[]
        self.list_degrees_log=[]
        self.list_count_log=[]

        self.list_degrees_for_clustering_log=[]
        self.list_clustering_coeff_log=[]
        
        try:
            self.graph = self.adj_mat_to_graph(mat)
        except:
            self.graph = nx.fast_gnp_random_graph(nb_nodes,0.2)
            
    def graph_to_adj_mat(self):
        return np.array(nx.adjacency_matrix(self.graph)).astype(int)
    
    def adj_mat_to_graph(self, mat):
        return nx.from_numpy_matrix(mat)
    
    def graphizer(self, label, i):
        update_progress("Plotting {}".format(label),i)
        try:
            nx.draw_graphviz(self.graph)
        except:
            nx.draw(self.graph)
            if not WARNVIZ:
                print WARNING+"Would be best with graphviz"+ENDC
                WARNVIZ = True
        b=plt.savefig(IMG+label+"_"+self.id+".png") # save as png
        plt.clf()

    def degree_graph(self,label,i):
        update_progress("Plotting {}".format(label),i)
        plt.plot(self.list_degrees_log,self.list_count_log)
        plt.savefig(IMG+label+"_"+self.id+".png") # save as png
        plt.clf()
           
    def clique_graph(self,generation,i):
        plt.plot(self.list_degrees_for_clustering_log,self.list_clustering_coeff_log)
        plt.savefig(IMG+"CG_gen"+str(generation)+"_id"+str(i)+"_graph"+str(self.id)+".png") # save as png
        plt.clf()

    
    def apply_mutations(self):
        m = self.graph_to_adj_mat()
        #Substitutions
        for x in range(int(RATE_SUB*NB_NODES)):
            #Attention petit biais : Un locus peut muter 2 fois et donc rester stable, affaiblissement 
            #du taux de mutation
            i = rd.randint(0,NB_NODES-1) #Je m'interroge sur la pertinence de ce -1
            j = rd.randint(i,NB_NODES-1)
            m = substitution(m,i,j)

        #Local (One node concerns) insertions and deletions.
        for x in range(int(RATE_LOCAL_INS*NB_NODES)):
            i = rd.randint(0,NB_NODES-1) #Je m'interroge toujours sur la pertinence de ce -1
            j = rd.randint(i,NB_NODES-1)
            m = insertion(m,i,j,rd.randint(0,1))

        for x in range(int(RATE_LOCAL_DEL*NB_NODES)):
            i = rd.randint(0,NB_NODES-1) #Je m'interroge encore sur la pertinence de ce -1
            j = rd.randint(i,NB_NODES-1)
            m = deletion(m,i,j)

        #Global insertions and deletions.
        if rd.random() < RATE_GLOBAL_INS:
            i = rd.randint(0,np.sum(np.arange(NB_NODES))) 
            m = ins_in_compr(m,i,rd.randint(0,1))

        if rd.random() < RATE_GLOBAL_DEL:
            i = rd.randint(0,np.sum(np.arange(NB_NODES))) 
            m = del_in_compr(m,i)

        return m

    def power_degree_law(self):
        """ power degree law """

        #if len(self.list_degrees_log) <= 0.05*(NB_NODES) :
        #    self.penalite+=200
        #else :
            #if len(self.list_degrees_log) < 0.15*(NB_NODES) :
            #    self.penalite+=50
            #    print len(self.list_degrees_log)
        nb_diff_deg = int(0.8*len(self.list_degrees_log))
        if nb_diff_deg < 4:
            nb_diff_deg = len(self.list_degrees_log)
            self.penalite += 10
            
        slope=stats.linregress(self.list_degrees_log[:nb_diff_deg],self.list_count_log[:nb_diff_deg])

        SCE=(slope[4]**2)*NB_NODES
        
            #print "SCE = " + str(SCE*5) +" pente = " + str(slope[0]) + "   erreur de pente : " + str(abs(-1-slope[0])*15) +"\n"
            

        if slope[0] > 0 : self.penalite += 100
        
        self.score_pdl = abs(-2-slope[0])+SCE*10

        
        
    def clique_formation(self):
        """ Compute clique formation score through log(linear) regression 

        lin_regress_clique=stats.linregress(self.list_degrees_for_clustering_log,self.list_clustering_coeff_log)
        #print "\n liste des degrés en log : " +str(self.list_degrees_for_clustering_log)
        #print "liste des degrés : " + str(self.deg_dict.values())
        #print "\n liste des coeff de clustering : " +str(self.clustl_dict.values()) + "\n"


        SCE_clique=(lin_regress_clique[4]**2)*NB_NODES

        if lin_regress_clique[0] > 0 : self.penalite += 20 
        self.score_cf = abs(-0.5-lin_regress_clique[0])+SCE_clique
        """

        tri = np.mean(nx.triangles(self.graph).values())
        self.score_cf=1/(tri+EPS)
        #self.score_cf=abs((self.list_degrees[-1::])[0]-25)*5

                
    def small_world(self):
        """ Compute small world score of graph """
        L = nx.average_shortest_path_length(self.graph)
        #C = nx.average_clustering(self.graph)

        self.score_sw=abs(L-L_RAND)
        if (self.score_sw*SW >6) :
            self.penalite+=100
        """
        try:
            a= nx.hits(self.graph)[0]
            print a
            print self.deg_dict

            plt.plot(self.deg_dict.values(),a.values() , 'ro')
            plt.axis([0, 40, 0, 0.1])
            plt.savefig(IMG+self.id+"hub")
            plt.clf()

        except:
            print FAIL+"pas réussi"+ENDC
        """   
        #self.score_sw = (1-C)*L # A préciser !
        # self.score_sw=abs(1-S)*50
        
    def reconnect(self,main,sub):
        recon = range(int(round(0.4*len(sub),0))) if len(sub) != 1 else range(1)
        new_edge = []
        for i in recon:
            new_edge.append((rd.sample(main,1)[0],rd.sample(sub,1)[0]))
        return new_edge
    

    def score_toolbox(self):
        # Correction for unconnected_graph
        if not nx.is_connected(self.graph):
            subnods = nx.connected_components(self.graph)
            map(lambda x:self.graph.add_edges_from(self.reconnect(subnods[0],x)) ,subnods[1:len(subnods)])

        # Usefull dict
        self.deg_dict=self.graph.degree() #dictionnaire des degrés : clé = id noeuds ; valeur = degré
        self.clustl_dict = nx.clustering(self.graph) # key = id nodes, values = local clustering 

        self.list_degrees = list(set(self.deg_dict.values())) # [unique dict values]

        # Lists
        values = sorted((self.deg_dict.values()))
        self.list_count = [values.count(x) for x in self.list_degrees]
        # Log
        self.list_degrees_log = [math.log10(x+EPS) for x in self.list_degrees]
        self.list_count_log = [math.log10(x+EPS) for x in self.list_count]
        self.list_degrees_for_clustering_log=[math.log10(x+EPS) for x in self.deg_dict.values()]
        self.list_clustering_coeff_log=[math.log10(x+EPS) for x in self.clustl_dict.values()]


    def calc_score(self):
        """ Fitness function """
        self.score_pdl = 1
        self.score_sw = 1
        self.score_cf = 1
        self.penalite = 0

        self.score_toolbox()
            
        # Score functions
        self.power_degree_law()
        self.small_world()
        self.clique_formation()

        self.score = PDL*self.score_pdl + SW*self.score_sw + CF*self.score_cf + self.penalite

        return self.score
