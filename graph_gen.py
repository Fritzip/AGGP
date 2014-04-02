#! /usr/bin/python
# -*- coding: utf-8 -*-

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random as rd
import copy
import math
import sys
import time
from scipy import stats as stats

####################################################################
#			Global Parameters
####################################################################
while True:
    try:
        # Plot information box
        INFO_INDIV = False
        INFO_BEST = False
        INFO_SELECT = False
        INFO_GEN = True
        INFO_FREQ = 1 # information frequency (every X generation)
        
        # Plot png in /img
        PLOT_DG = 1001 # plot degree graph every X generation
        PLOT_GR = 10 # plot graph every X generation
        PLOT_GEN_ZERO = True

        # Parameters of algogen
        NB_GEN = 1000 # genetic algo's iteration number
        NB_NODES = 25
        NB_INDIV = 20
        
        # Rates
        RATE_ELITISM = 0.2 
        RATE_TOURNAMENT = 1-RATE_ELITISM
        RATE_CROSS = 0.7
        RATE_MUT = 0.2
        
        # Random Reference
        G_RAND = nx.fast_gnp_random_graph(NB_NODES,0.2)
        C_RAND = nx.average_clustering(G_RAND)
        L_RAND = nx.average_shortest_path_length(G_RAND)
        
        # Miscellanous
        NAMES = open("names.txt").read().splitlines()
        ERROR = False
        break
    except:
        pass

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
    return -a*a.T + a + a.T

def update_progress(label,progress,bar_length=25): # small 20, medium 25, large 50
    progress = int(progress)
    if progress > 100 : progress = 100
    sys.stdout.write('\r{2:<25}[{0}]{1:3d}%'.format('#'*(progress/int(100./bar_length))+'-'*(bar_length-(progress/int(100./bar_length))), progress,label))
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
        
        self.generation = 1
        
        self.size_pop = size_pop
        self.size_indiv = size_indiv
        self.nb_best = int(self.size_pop*RATE_ELITISM)+1
        
        for i in range(self.size_pop):
            self.indiv.append(Individual(nb_nodes=self.size_indiv,id=rd.choice(NAMES)))
            self.score.append(0)
        
    
    def genetic_algo(self):
        start_algo = time.time()
        
        while self.generation<NB_GEN and not ERROR:
            try:
                start = time.time()
                self.evaluation()
                self.selection()
                self.crossormut()
                self.indiv = []
                self.indiv = copy.deepcopy(self.next_gen)
                self.next_gen = []
                
                # INFO GENERATION BOX
                if INFO_GEN and self.generation%INFO_FREQ==0:
                    print "\n╔{0}╗\n║ {2}{1}{3} ║".format('═'*39,"GENERATION {0}".format(self.generation).center(37),HEADER,ENDC)
                    print "║ {2}{0:6.2f} ┆ {1:<10} : {4:<15}{3} ║".format(min(self.selected_score),"Best Score",OKGREEN,ENDC,self.best_ever_indiv[0].id)
                    print "║ {2}{0:6.2f} ┆ {1:<28}{3} ║".format(float(sum(self.selected_score))/len(self.selected_score), "Mean Score",OKBLUE,ENDC)
                    print "║ {2}{0:6.2f} ┆ {1:<28}{3} ║".format(max(self.selected_score), "Worst Score",FAIL,ENDC)
                    print "╚{0}╝".format('═'*39)
                    
                    if self.generation%PLOT_GR==0 or (self.generation==1 and PLOT_GEN_ZERO):
                        i=1
                        for indi in self.indiv :
                            i += 100./len(self.indiv)
                            indi.graphizer(self.generation,i)
                
                self.generation += 1
            except KeyboardInterrupt:
                print "\nKeyboard Interrupt"
                break
            
        print "Done in %.3f sec"%(time.time()-start_algo)
        # en sortie de l'algorithme : lancer des plots, des stats, des summary, des feux d'artifices de pop-up…
        
            
    def evaluation(self):
        """ Compute the score using fitness function """
        count = 0
        for i in range(len(self.indiv)):
            try:
                self.score[i] = self.indiv[i].calc_score(self.generation,count)
                count+=1
            except nx.NetworkXError:
                self.score[i] = 100
                self.indiv[i].score = 100
            if INFO_INDIV and self.generation%INFO_FREQ==0:
                indi = self.indiv[i]
                print "+{}+".format('-'*30)
                print "|{}|".format(indi.id.center(30))
                print "| {0:6.2f} | {1:<19} |".format(indi.score_pdl,"Power Degree Law")
                print "| {0:6.2f} | {1:<19} |".format(indi.score_sw,"Small World")
                print "| {0:6.2f} | {1:<19} |".format(indi.score,"Global")
                print "+{}+\n".format('-'*30)

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
        if INFO_BEST and self.generation%INFO_FREQ==0:
            ever = map(lambda x: x.id ,self.best_ever_indiv)
            last = map(lambda x: x.id ,self.best_last_gen_indiv)
                                    
            LONG = 50
            print "\n╔{0}╦{0}╗".format('═'*((LONG-3)/2))
            print "║{0}║{1}║".format("BEST EVER".center((LONG-3)/2),"BEST LAST GEN".center((LONG-3)/2))
            print "╟{0}╫{0}╢".format('─'*((LONG-3)/2))
            
            for i in range(len(self.best_ever_score)):
                print "║{0:6.2f} ┆ {1:<14}║{2:6.2f} ┆ {3:<14}║".format(self.best_ever_score[i],ever[i],self.best_last_gen_score[i],last[i])

            print "╚{0}╩{0}╝\n".format('═'*((LONG-3)/2))
        
        if INFO_SELECT and self.generation%INFO_FREQ==0:
            select = map(lambda x: x.id ,self.selected_indiv)
            
            LONG = 75
            FAC = len(self.selected_indiv)/3

            print "\n╔{0}╗".format('═'*(LONG-4))
            print "║{0}║".format("SELECTED".center(LONG-4))
            print "╟{0}┬{0}┬{0}╢".format('─'*((LONG-4)/3))
            for i in range(FAC):
                print "║{0:6.2f} ┆ {1:<14}│{2:6.2f} ┆ {3:<14}│{4:6.2f} ┆ {5:<14}║".format(self.selected_score[i],select[i],self.selected_score[i+FAC],select[i+FAC],self.selected_score[i+FAC*2],select[i+FAC*2])
            if NB_INDIV%3==1:
                print "║{0:6.2f} ┆ {1:<14}│{2:<6} ┆ {3:<14}│{4:<6} ┆ {5:<14}║".format(self.selected_score[-1],select[-1],'','','','')
            elif NB_INDIV%3==2:
                print "║{0:6.2f} ┆ {1:<14}│{2:6.2f} ┆ {3:<14}│{4:<6} ┆ {5:<14}║".format(self.selected_score[-1],select[-1],self.selected_score[-2],select[-2],'','')
            print "╚{0}╧{0}╧{0}╝".format('═'*((LONG-4)/3))
            
                

    def tournament(self,a,b):
        if self.score[a]>=self.score[b]:
            higher = a
            lower = b
        else:
            higher = b
            lower = a
        
        rand = rd.random()
        prob = 0.5+((self.score[higher]-self.score[lower])/self.score[higher])*0.5
        #print "prob = %.2f, rand  = %.2f"%(prob,rand)
        #print self.indiv[higher].id,self.indiv[lower].id
        if rand > prob: return higher
        else: return lower
        
        
    def selection(self):
        """ Choose type of selection """
        #self.roulette_wheel_selec()
        self.elitist_selec()
        
    def crossormut(self):
        while len(self.selected_indiv) != 0:
            rand = rd.random()
            if rand < RATE_CROSS and len(self.selected_indiv) > 1:
                sample = rd.sample(self.selected_indiv,2)
                map(self.trans_indiv ,sample)
                self.cross(tuple(sample))
            elif rand < RATE_MUT+RATE_CROSS:
                sample = rd.sample(self.selected_indiv,1)[0]
                self.trans_indiv(sample)
                self.mutation(sample)
            else :
                sample = rd.sample(self.selected_indiv,1)[0]
                self.trans_indiv(sample)
                self.next_gen.append(sample)
                    
    def trans_indiv(self,x):
        self.selected_indiv.remove(x)
        
    def cross(self,(ind1,ind2)):
        """ Simulate the reproduction beetween two individuals : ~70% of variability """
        A = ind1.graph_to_adj_mat()
        B = ind2.graph_to_adj_mat()
        X = np.array(nx.adjacency_matrix(nx.fast_gnp_random_graph(self.size_indiv,0.2))).astype(int)
        self.next_gen.append(Individual(mat=((A+B)%2)*X+(A+B)/2,id=rd.choice(NAMES)))
        self.next_gen.append(Individual(mat=((A+B)%2)*((X+np.ones((self.size_indiv,self.size_indiv),dtype=np.int))%2)+(A+B)/2,id=rd.choice(NAMES)))
        
    def mutation(self,ind):
        """~20% of variability """
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
    
    def graphizer(self, gen, i):
        update_progress("Plotting Generation {}".format(gen),i)
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
        if generation%10000==0:
            print ("\n" + str(self.id))
            print ("id="+str(i))
            #print "liste des degrés :"
            #print list_degrees
            #print "list des count  : "
            #print list_count
            print slope
            print "erreur de pente : "+str(abs(-2.5-slope[0])*10)
            print "SCE : " + str(SCE)
        
        
    def small_world(self):
        """ small world """
        L = nx.average_shortest_path_length(self.graph)
        C = nx.average_clustering(self.graph)
        self.score_sw = 1-C+L
        
    
    def calc_score(self,generation,i):
        """ Fitness function """
        self.score_pdl = 0
        self.score_sw = 0

        #self.power_degree_law(generation,i)
        self.small_world()

        self.score = self.score_sw + self.score_pdl

        return self.score

####################################################################
#			Main
####################################################################
if __name__ == "__main__":
    a = Population()
    a.genetic_algo()
