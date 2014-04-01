#! /usr/bin/python
# -*- coding: utf-8 -*-

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random as rd
import copy
import math
import sys
from scipy import stats as stats

####################################################################
#			Global Parameters
####################################################################
while True:
    try:
        PLOT_DG = 101 # plot degree graph every X generation
        PLOT_GR = 10 # plot graph every X generation
        NB_GEN = 100
        NB_NODES = 25
        NB_INDIV = 20
        G_RAND = nx.fast_gnp_random_graph(NB_NODES,0.2)
        C_RAND = nx.average_clustering(G_RAND)
        L_RAND = nx.average_shortest_path_length(G_RAND)
        NAMES = open("names.txt").read().splitlines()
        ERROR = False
        break
    except:
        pass


####################################################################
#			Global Functions
####################################################################
def symetrize(a):
    return -a*a.T + a + a.T

def update_progress(progress,bar_length=20):
    sys.stdout.write('\r[{0}] {1}%'.format('#'*(progress/int(100./bar_length))+' '*(bar_length-(progress/int(100./bar_length))), progress))
    sys.stdout.flush()

def weighted_sample(items, n):
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

####################################################################
#			Gaïa
####################################################################
class Gaia():
    """ Run the genetic algorithm and select the individuals """
    def __init__(self):
        pass
    
####################################################################
#			Population
####################################################################
class Population():
    """ Population including all the individuals """
    def __init__(self,size_pop=NB_INDIV,size_indiv=NB_NODES):
        self.indiv = []
        self.score = []
        
        self.selected_indiv = []
        self.selected_score = []
        
        self.best_last_gen_score = []
        self.best_last_gen_indiv = []
        
        self.best_ever_score = []
        self.best_ever_indiv = []
        
        self.next_gen = []
        
        self.generation = 0
        
        self.size_pop = size_pop
        self.size_indiv = size_indiv
        self.nb_best = int(self.size_pop*0.2)+1

        for i in range(self.size_pop):
            self.indiv.append(Individual(nb_nodes=self.size_indiv,id=rd.choice(NAMES)))
            self.score.append(0)
            
    def genetic_algo(self):
        """ Inspired by
        http://fr.wikipedia.org/wiki/Algorithme_g%C3%A9n%C3%A9tique#Sch.C3.A9ma_r.C3.A9capitulatif
        """
        
        while self.generation<NB_GEN and not ERROR:
            print "\n###############################\n##\tGénération = %d\n###############################\n"%self.generation
            if self.generation%PLOT_GR==0:
                i=0
                for indi in self.indiv : 
                    indi.graphizer(self.generation,i)
                    i += int(100./len(self.indiv))
            self.evaluation()
            self.selection()
            self.crossormut()
            self.indiv = []
            self.indiv = copy.deepcopy(self.next_gen)
            self.next_gen = []
            self.generation += 1
            
    def evaluation(self):
        """ Compute the score using fitness function """
        count = 0
        for i in range(len(self.indiv)):
            try:
                self.score[i] = self.indiv[i].calc_score(self.generation,count)
                count+=1
            except nx.NetworkXError:
                self.score[i] = 100
            print "%.2f\t%s"%(self.score[i],self.indiv[i].id) #"\n",indi.graph_to_adj_mat()
            print "\n==================================="
            
    def roulette_wheel_selec(self):
        """ Roulette wheel selection """
        zipper = zip(self.score,self.indiv)
        generator = weighted_sample(zipper,self.size_pop/2)
        for i in generator:
            self.selected_indiv.append(i)
            #	print "selected =",self.selected_indiv
            #	print "len(selected) =",len(self.selected_indiv)
            
    def elitist_selec(self):
        """ Elitist (20%) and tournament selection (80%) """
        #self.best_last_gen_score, self.best_last_gen_indiv = map(lambda x : list(x[0:self.nb_best]),zip(*(sorted(zip(self.score,self.indiv)))))
        #self.best_ever_score, self.best_ever_indiv = copy.deepcopy(map(lambda x : list(x[0:self.nb_best]),zip(*(sorted(zip(self.best_last_gen_score+self.best_ever_score,self.best_last_gen_indiv+self.best_ever_indiv))))))

        self.best_last_gen_score, self.best_last_gen_indiv = map(lambda x : list(x),zip(*(sorted(zip(self.score,self.indiv)))))
        self.best_temp_score, self.best_temp_indiv = copy.deepcopy(map(lambda x : list(x),zip(*(sorted(zip(self.best_last_gen_score+self.best_ever_score,self.best_last_gen_indiv+self.best_ever_indiv))))))

        self.indices_best = [0]
        indice = 1
        while len(self.indices_best) != self.nb_best or indice >= len(self.best_temp_indiv):
            if not np.array_equal(self.best_temp_indiv[indice].graph_to_adj_mat(), self.best_temp_indiv[self.indices_best[-1]].graph_to_adj_mat()):
                self.indices_best.append(indice)
            indice += 1
        if len(self.indices_best) != self.nb_best:
            ERROR = True

        self.best_ever_indiv = map(lambda x: self.best_temp_indiv[x] ,self.indices_best)
        self.best_ever_score = map(lambda x: self.best_temp_score[x] ,self.indices_best)

        self.selected_indiv = copy.deepcopy(self.best_ever_indiv)
        self.selected_score = copy.deepcopy(self.best_ever_score)
        while len(self.selected_indiv) != self.size_pop:
            a,b = rd.sample(range(len(self.indiv)),2)
            indice = self.tournament(a,b)
            self.selected_indiv.append(self.indiv[indice])
            self.selected_score.append(self.score[indice])
        
        ### Print, just print
        ever = map(lambda x: x.id ,self.best_ever_indiv)
        last = map(lambda x: x.id ,self.best_last_gen_indiv)
        select = map(lambda x: x.id ,self.selected_indiv)
        print ever, last, select
        print self.best_ever_score,self.best_last_gen_score

        print "\n####################################################\n## BEST EVER \t\t BEST LAST GENERATION\n####################################################"
        for i in range(len(self.best_ever_score)):
            print "%.2f\t%s\t | \t%.2f\t%s"%(self.best_ever_score[i],ever[i],self.best_last_gen_score[i],last[i])
            
        print "\n##########################\n## SELECTED\n##########################"
        for i in range(len(self.selected_indiv)):
            print "%.2f\t%s"%(self.selected_score[i],select[i])
        
    def tournament(self,a,b):
        if self.score[a]>=self.score[b]:
            higher = a
            lower = b
        else:
            higher = b
            lower = a
        
        rand = rd.random()
        prob = 0.5+((self.score[higher]-self.score[lower])/self.score[higher])*0.5
        print "prob = %.2f, rand  = %.2f"%(prob,rand)
        print self.indiv[higher].id,self.indiv[lower].id
        print self.score[higher]==self.score[lower]
        print self.indiv[higher]==self.indiv[lower]
        print np.array_equal(self.indiv[higher].graph_to_adj_mat(),self.indiv[lower].graph_to_adj_mat())
        if rand > prob: return higher
        else: return lower
        
        
    def selection(self):
        """ Choose type of selection """
        #self.roulette_wheel_selec()
        self.elitist_selec()
        
    def crossormut(self):
        while len(self.selected_indiv) != 0:
            rand = rd.random()
            if rand < 0.7 and len(self.selected_indiv) > 1:
                sample = rd.sample(self.selected_indiv,2)
                map(self.trans_indiv ,sample)
                self.cross(tuple(sample))
            elif rand < 0.9:
                sample = rd.sample(self.selected_indiv,1)[0]
                self.trans_indiv(sample)
                self.mutation(sample)
            else :
                #print "rien"
                sample = rd.sample(self.selected_indiv,1)[0]
                self.trans_indiv(sample)
                self.next_gen.append(sample)
                    
    def trans_indiv(self,x):
        #self.next_gen.append(x)
        self.selected_indiv.remove(x)
        
    def cross(self,(ind1,ind2)):
        #print "croisement"
        """ Simulate the reproduction beetween two individuals : 70% of variability """
        A = ind1.graph_to_adj_mat()
        B = ind2.graph_to_adj_mat()
        #X = symetrize(np.random.binomial(1, 0.2, self.size_indiv**2).reshape(self.size_indiv,self.size_indiv))
        X = np.array(nx.adjacency_matrix(nx.fast_gnp_random_graph(self.size_indiv,0.2))).astype(int)
        self.next_gen.append(Individual(mat=((A+B)%2)*X+(A+B)/2,id=rd.choice(NAMES)))
        self.next_gen.append(Individual(mat=((A+B)%2)*((X+np.ones((self.size_indiv,self.size_indiv),dtype=np.int))%2)+(A+B)/2,id=rd.choice(NAMES)))
        
    def mutation(self,ind):
        #print "mutation"
        """ 30% of variability """
        self.next_gen.append(Individual(mat = ind.apply_mutations(), id=rd.choice(NAMES)))
        
####################################################################
#			Individual
####################################################################
class Individual():
    """ Metabolic graph """
    def __init__(self,mat=0, nb_nodes=NB_NODES,id=rd.choice(NAMES)):
        self.id = id
        self.score_pdl = 0
        self.score_sw = 0
        self.score = 0
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
    
    def graphizer(self, gen,i):
        update_progress(i)
        nx.draw(self.graph)
        b=plt.savefig("img/gen"+str(gen)+"_graph"+str(self.id)+".png") # save as png
        plt.clf()

    def degree_graph(self,list_degrees_log,list_count_log,generation,i):
        if generation%PLOT_DG==0 :
            a=plt.plot(list_degrees_log,list_count_log)
            plt.savefig("img/plot"+"_gen"+str(generation)+"_id"+str(i)+"_graph"+str(self.id)+".png") # save as png
            plt.clf()
        
    def apply_mutations(self):
        return self.graph_to_adj_mat()

    def power_degree_law(self,generation,i):
        """ power degree law """
        dict=self.graph.degree() #dictionnaire des degrés : clé = id noeuds ; valeur = degré
        #print dict
        values = sorted((dict.values()))
	#print values
        list_degrees=[]
        for x in values :
            if x not in list_degrees and x != 0 :
                list_degrees.append(x)
                
        list_count = [values.count(x) for x in list_degrees]
        
        list_degrees_log=[math.log10(x) for x in list_degrees]
        list_count_log=[math.log10(x) for x in list_count]
        #print "liste des degrés en log :"
        self.degree_graph(list_degrees_log,list_count_log,generation,i)
        
        slope=stats.linregress(list_degrees_log,list_count_log)
        #print slope
        SCE=(slope[4]**2)*25

        if nx.is_connected(self.graph) == False:
            sanction = 100
        else :
            sanction=0
        if slope[0] > 0 :
            sanction_pente = 20
        else : 
            sanction_pente=0
        
        self.score_pdl = abs(-3-slope[0])*10+SCE+sanction+sanction_pente
        if generation%100==0:
            print ("\n" + str(self.id))
            print ("id="+str(i))
            #print "liste des degrés :"
            #print list_degrees
            #print "list des count  : "
            #print list_count
            print slope
            print "erreur de pente : "+str(abs(-2.5-slope[0])*10)
            print "SCE : " + str(SCE)
        print "score_pdl = ", self.score_pdl
        
    def small_world(self):
        """ small world """
        L = nx.average_shortest_path_length(self.graph)
        C = nx.average_clustering(self.graph)
        self.score_sw = (1-C)+(abs(L-L_RAND))
        print "score sw = ",self.score_sw
    
    def calc_score(self,generation,i):
        """ Fitness function """
        print "\n## %s ##"%self.id
        self.score_pdl = 0
        self.score_sw = 0

        self.power_degree_law(generation,i)
        self.small_world()

        self.score = self.score_sw + self.score_pdl
        print "score global =",self.score
        return self.score

####################################################################
#			Main
####################################################################
if __name__ == "__main__":
    a = Population()
    a.genetic_algo()
