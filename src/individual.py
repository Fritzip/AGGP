#! /usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import random as rd
import copy, math, sys, time
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
        
        self.list_degrees_log=[]
        self.list_count_log=[]
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
        nx.draw_graphviz(self.graph)
        b=plt.savefig(IMG+label+"_"+self.id+".png") # save as png
        plt.clf()

    def average_short_path(self):
        return nx.average_shortest_path_length(self.graph)

    def average_coeff_clustering(self):
        return nx.average_clustering(self.graph) 

    def degree_graph(self,label,i):
        update_progress("Plotting {}".format(label),i)
        plt.plot(self.list_degrees_log,self.list_count_log)
        plt.savefig(IMG+label+"_"+self.id+".png") # save as png
        plt.clf()

    """        
    def clique_graph(self,generation,i):
        plt.plot(self.list_degrees,self.list_meanclust_log)
        plt.savefig(IMG+"CG_gen"+str(generation)+"_id"+str(i)+"_graph"+str(self.id)+".png") # save as png
        plt.clf()
    """          
    def apply_mutations(self):
        m = self.graph_to_adj_mat()
        for x in range(int(0.3*NB_NODES)):
            i = rd.randint(0,NB_NODES-1)
            j = rd.randint(i,NB_NODES-1)
            m = substitution(m,i,j)
        return m

    def power_degree_law(self,generation,i):
        """ power degree law """
        #if generation%PLOT_PDL==0 :
        #    self.degree_graph(generation,i) # Plot
        
        slope=stats.linregress(self.list_degrees_log,self.list_count_log)
        
        SCE=(slope[4]**2)*NB_NODES

        if slope[0] > 0 : self.penalite += 20
            
        self.score_pdl = abs(-2.5-slope[0])*10+SCE*10

        if generation%1001==0:
            print ("\n" + str(self.id))
            print ("id="+str(i))
            print "liste des degrés :"
            print self.list_degrees
            print "list des count  : "
            print self.list_count
            print slope
            print "erreur de pente : "+str(abs(-1.5-slope[0])*10)
            print "SCE : " + str(SCE)
            print "score :",self.score_pdl
        
        
    def clique_formation(self,generation,i):
        """ Compute clique formation score through log(linear) regression """
        """
        #if generation%PLOT_CF==0 :
        #    self.clique_graph(generation,i) # Plot
        
        slope=stats.linregress(self.list_degrees,self.list_meanclust_log)
        SCE=(slope[4]**2)*NB_NODES

        if slope[0] > 0 : self.penalite += 20 
        self.score_cf = abs(-2-slope[0])*10+SCE
        
        """
        moysup = np.mean(self.clustsup) if len(self.clustsup)>0 else 0#sum(self.clustsup)/len(self.clustsup)
        moyinf = np.mean(self.clustinf) if len(self.clustsup)>0 else 0#sum(self.clustinf)/len(self.clustinf)
        #print "moysup = ",(moysup, len(self.clustsup))
        #print "moyinf = ",(moyinf, len(self.clustinf))
        #if len(self.clustsup)<int(0.2*NB_NODES):self.penalite += 20
        self.score_cf = (1-moysup)+moyinf+ abs(0.3*NB_NODES-len(self.clustsup))*0.5
                
    def small_world(self):
        """ Compute small world score of graph """
        L = nx.average_shortest_path_length(self.graph)
        C = nx.average_clustering(self.graph)
        S=(C/C_RAND)/(L/L_RAND)   
        #self.score_sw = (1-C)*L # A préciser !
        self.score_sw=abs(1-S)*20
        
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
        self.clustl_dict = nx.clustering(self.graph) # key = id nodes, values = local clustering value

        self.list_degrees = list(set(self.deg_dict.values())) # [unique dict values]

        """
        # Build clustfk (clustering~k)
        self.clustfk = {} # {k:mean(local_clust_of_k_deg_nodes)}
        for x in self.list_degrees : self.clustfk[x]=[]
        zipp = zip(self.deg_dict.values(), self.deg_dict.keys()) 
        map(lambda x: self.clustfk[x[0]].append(self.clustl_dict[x[1]]), zipp)
        for x in self.clustfk.keys() : self.clustfk[x]=np.mean(self.clustfk[x])
        """
        #a = np.mean(self.clustl_dict.values)
        self.clustsup = filter(lambda x: x>0.45, self.clustl_dict.values())
        self.clustinf = filter(lambda x: x<=0.45, self.clustl_dict.values())

        # Lists
        values = sorted((self.deg_dict.values()))
        self.list_count = [values.count(x) for x in self.list_degrees]
        #self.list_meanclust = self.clustfk.values()

        # Log
        self.list_degrees_log = [math.log10(x+EPS) for x in self.list_degrees]
        self.list_count_log = [math.log10(x+EPS) for x in self.list_count]
        #self.list_meanclust_log = [math.log10(x+EPS) for x in self.list_meanclust]
        
    def calc_score(self,generation,i):
        """ Fitness function """
        self.score_pdl = 1
        self.score_sw = 1
        self.score_cf = 1
        self.penalite = 0

        self.score_toolbox()
            
        # Score functions
        self.power_degree_law(generation,i)
        self.small_world()
        self.clique_formation(generation,i)

        self.score = PDL*self.score_pdl + SW*self.score_sw + CF*self.score_cf + self.penalite

        return self.score
