#! /usr/bin/python
# -*- coding: utf-8 -*-

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random as rd
import copy, math, sys, time, threading
from scipy import stats as stats

import mutations
####################################################################
#			Global Parameters
####################################################################
while True:
    try:
        # Plot png in /img
        PLOT_PDL = 1001 # plot degree graph every X generation
        PLOT_CF = 1001 # plot clique formation graph every X generation
        PLOT_GR = 1001 # plot graph every X generation
        PLOT_GEN_ZERO = False # plot initials individuals ?

        # Parameters of algogen
        NB_GEN = 1000 # genetic algo's iteration number
        NB_NODES = 60
        NB_INDIV = 20

        # Plot information box
        INFO_INDIV = False
        INFO_BEST = True
        INFO_SELECT = False
        INFO_GEN = True
        INFO_FREQ = 10 # information frequency (every X generation)
        PROGRESS_GEN = False if INFO_INDIV or NB_NODES*NB_INDIV < 1000 else True

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
        G_RAND = nx.fast_gnp_random_graph(NB_NODES,0.2)
        C_RAND = nx.average_clustering(G_RAND)
        L_RAND = nx.average_shortest_path_length(G_RAND)
        
        # Miscellaneous
        NAMES = open("names.txt").read().splitlines()
        ERROR = False
        EPS = 0.001 # log(x+EPS) to avoid log(0)

        # Colors
        HEADER = '\033[1m' # bold
        OKBLUE = '\033[94m' # blue
        OKGREEN = '\033[92m' # green
        WARNING = '\033[93m' # yellow
        FAIL = '\033[91m' # red
        ENDC = '\033[0m' # back to normal
        break
    except:
        pass


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

class Progressbar(threading.Thread):
    def __init__(self,gen,time_last):
        threading.Thread.__init__(self)
        self.gen = gen
        self.time_last = time_last
        self.remaining = time_last*0.95
        self.zero = time.time()
        self.stopped = False
        
    def stop(self):
        self.stopped = True
           
    def run(self):
        while self.remaining > 0:
            if not self.stopped :#1self.stopped():
                self.remaining = self.time_last-(time.time()-self.zero)
                i = int(100-(100.*self.remaining/self.time_last))
                update_progress("Generation "+str(self.gen),i)
                time.sleep(0.05)
            else:
                break
   

class Population():
    """ Population including all the individuals """
    def __init__(self,size_pop=NB_INDIV,size_indiv=NB_NODES):
        self.indiv = []
        self.score = []
        self.score_pdl = []
        self.score_sw = []
        self.score_cf = []
        
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
            self.score_pdl.append(0)
            self.score_sw.append(0)
            self.score_cf.append(0)
        
    
    def genetic_algo(self):
        start_algo = time.time()
        time_laps = 0

        # Open all files
        self.fscore = open('evo_score','w')
        self.fpdl = open('evo_pdl','w')
        self.fsw = open('evo_sw','w')
        self.fcf = open('evo_cf','w')

        # Genetic algorithm
        while self.generation<NB_GEN and not ERROR:
            try:
                if PROGRESS_GEN and self.generation != 0:
                    bar = Progressbar(self.generation,time_laps)
                    bar.start()
                start = time.time()
                self.evaluation()
                self.save() # in files
                self.selection()
                self.crossormut()
                time_laps = time.time()-start
                self.indiv = []
                self.indiv = copy.deepcopy(self.next_gen)
                self.next_gen = []
                
                self.prints()
                
                self.generation += 1

            except KeyboardInterrupt:
                bar.stop()
                print "\nKeyboard Interrupt"
                break
        
        # Close all files
        self.fscore.close()
        self.fpdl.close()
        self.fsw.close()
        self.fcf.close()

        self.best_ever_indiv[0].graphizer("Best Indiv",100)
        self.save2sif(self.best_ever_indiv[0].graph_to_adj_mat())
        print "\nDone in %.3f sec"%(time.time()-start_algo)
        # en sortie de l'algorithme : lancer des plots, des stats, des summary, des feux d'artifices de pop-up…
        
            
    def evaluation(self):
        """ Compute the score using fitness function """
        count = 0
        for i in range(len(self.indiv)):
            self.score[i] = self.indiv[i].calc_score(self.generation,count)
            self.score_pdl[i] = self.indiv[i].score_pdl
            self.score_sw[i] = self.indiv[i].score_sw
            self.score_cf[i] = self.indiv[i].score_cf
            count+=1
            
            # PRINT
            if INFO_INDIV and self.generation%INFO_FREQ==0:
                indi = self.indiv[i]
                print "+{}+".format('-'*30)
                print "|{}|".format(indi.id.center(30))
                print "| {0:7.2f} | {1:<18} |".format(indi.score_pdl,"Power Degree Law")
                print "| {0:7.2f} | {1:<18} |".format(indi.score_sw,"Small World")
                print "| {0:7.2f} | {1:<18} |".format(indi.score_cf,"Clique Formation")
                print "| {0:7.2f} | {1:<18} |".format(indi.score,"Global")
                print "+{}+\n".format('-'*30)

    def roulette_wheel_selec(self):
        """ Roulette wheel selection """
        zipper = zip(self.score,self.indiv)
        generator = weighted_sample(zipper,self.size_pop/2)
        for i in generator:
            self.selected_indiv.append(i)
            
            
    def elitist_selec(self):
        """ Elitist (20%) and tournament selection (80%) """
        #self.best_last_gen_score, self.best_last_gen_indiv = map(lambda x : list(x[0:self.nb_best]),zip(*(sorted(zip(self.score,self.indiv)))))
        #self.best_ever_score, self.best_ever_indiv = copy.deepcopy(map(lambda x : list(x[0:self.nb_best]),zip(*(sorted(zip(self.best_last_gen_score+self.best_ever_score,self.best_last_gen_indiv+self.best_ever_indiv))))))

        self.best_last_gen_score, self.best_last_gen_indiv = map(lambda x : list(x),zip(*(sorted(zip(self.score,self.indiv)))))
        self.best_temp_score, self.best_temp_indiv = copy.deepcopy(map(lambda x : list(x),zip(*(sorted(zip(self.best_ever_score+self.best_last_gen_score,self.best_ever_indiv+self.best_last_gen_indiv))))))

        self.indices_best = [0]
        indice = 1
        while len(self.indices_best) != self.nb_best:
            if indice < len(self.best_temp_indiv):
                if not np.array_equal(self.best_temp_indiv[indice].graph_to_adj_mat(), self.best_temp_indiv[self.indices_best[-1]].graph_to_adj_mat()):
                    self.indices_best.append(indice)
                indice += 1
            else :
                ERROR = True # pour le moment. À éventuellement modifier (petite population !)
                print FAIL+"Pas d'individus suffisamment différents"+ENDC
                break

        self.best_ever_indiv = map(lambda x: self.best_temp_indiv[x] ,self.indices_best)
        self.best_ever_score = map(lambda x: self.best_temp_score[x] ,self.indices_best)

        self.selected_indiv = copy.deepcopy(self.best_ever_indiv)
        self.selected_score = copy.deepcopy(self.best_ever_score)
        while len(self.selected_indiv) != self.size_pop:
            a,b = rd.sample(range(len(self.indiv)),2)
            indice = self.tournament(a,b)
            self.selected_indiv.append(self.indiv[indice])
            self.selected_score.append(self.score[indice])
        
       
    def tournament(self,a,b):
        if self.score[a]>=self.score[b]:
            higher = a
            lower = b
        else:
            higher = b
            lower = a
        
        rand = rd.random()
        try:
            prob = 0.5+((self.score[higher]-self.score[lower])/self.score[higher])*0.5
        except:
            prob = 0.5
        #print "prob = %.2f, rand  = %.2f"%(prob,rand)
        if rand > prob: return higher
        else: return lower
        
        
    def selection(self):
        """ Choose type of selection """
        #self.roulette_wheel_selec()
        self.elitist_selec()
        
    def crossormut(self):
        self.pick_indiv = copy.deepcopy(self.selected_indiv)
        while len(self.pick_indiv) != 0:
            rand = rd.random()
            if rand < RATE_CROSS and len(self.pick_indiv) > 1: # Crossing
                #print OKGREEN+"croise"+ENDC
                sample = rd.sample(self.pick_indiv,2)
                map(self.trans_indiv ,sample)
                self.cross(tuple(sample))
            elif rand < RATE_MUT+RATE_CROSS: # Mutation
                #print FAIL+"mute"+ENDC
                sample = rd.sample(self.pick_indiv,1)[0]
                self.trans_indiv(sample)
                self.mutation(sample)
            else : # Nothing
                #print WARNING+"nothing"+ENDC
                sample = rd.sample(self.pick_indiv,1)[0]
                self.trans_indiv(sample)
                self.next_gen.append(sample)
                    
    def trans_indiv(self,x):
        self.pick_indiv.remove(x)
        
    def cross(self,(ind1,ind2)):
        """ Simulate the reproduction beetween two individuals : ~70% of variability """
        A = ind1.graph_to_adj_mat()
        B = ind2.graph_to_adj_mat()
        X = np.array(nx.adjacency_matrix(nx.fast_gnp_random_graph(self.size_indiv,0.2))).astype(int)
        self.next_gen.append(Individual(mat=((A+B)%2)*X+(A+B)/2,id=rd.choice(NAMES)))
        self.next_gen.append(Individual(mat=((A+B)%2)*((X+np.ones((self.size_indiv,self.size_indiv),dtype=np.int))%2)+(A+B)/2,id=rd.choice(NAMES)))
        
    def mutation(self,ind):
        """~20% of variability """
        self.next_gen.append(Individual(mat = ind.apply_mutations(), id=ind.id))#rd.choice(NAMES)))

    def save(self):
        self.fscore.write(str(self.generation)+'\t'+str(max(self.score))+'\t'+str(sum(self.score)/len(self.score))+'\t'+str(min(self.score))+'\n')
        self.fpdl.write(str(self.generation)+'\t'+str(max(self.score_pdl))+'\t'+str(sum(self.score_pdl)/len(self.score_pdl))+'\t'+str(min(self.score_pdl))+'\n')
        self.fsw.write(str(self.generation)+'\t'+str(max(self.score_sw))+'\t'+str(sum(self.score_sw)/len(self.score_sw))+'\t'+str(min(self.score_sw))+'\n')
        self.fcf.write(str(self.generation)+'\t'+str(max(self.score_cf))+'\t'+str(sum(self.score_cf)/len(self.score_cf))+'\t'+str(min(self.score_cf))+'\n')


    def save2sif(self,m):
        sif = open('cyto.sif','w')
        for j in range(NB_NODES):
            for i in range(j,NB_NODES):
                if m[i][j]==1:
                    sif.write(str(i)+'\tpp\t'+str(j)+'\n')
        sif.close()
        
    def prints(self):
        """ Print, just print """
        # INFO BEST INDIV
        if INFO_BEST and self.generation%INFO_FREQ==0:
            ever = map(lambda x: x.id ,self.best_ever_indiv)
            last = map(lambda x: x.id ,self.best_last_gen_indiv)
            
            LONG = 50
            print "\n╔{0}╦{0}╗".format('═'*((LONG-3)/2))
            print "║{0}║{1}║".format("BEST EVER".center((LONG-3)/2),"BEST LAST GEN".center((LONG-3)/2))
            print "╟{0}╫{0}╢".format('─'*((LONG-3)/2))
            
            for i in range(len(self.best_ever_score)):
                print "║{0:7.2f} ┆ {1:<13}║{2:7.2f} ┆ {3:<13}║".format(self.best_ever_score[i],ever[i],self.best_last_gen_score[i],last[i])

            print "╚{0}╩{0}╝\n".format('═'*((LONG-3)/2))
        
        # INFO SELECTED INDIV
        if INFO_SELECT and self.generation%INFO_FREQ==0:
            select = map(lambda x: x.id ,self.selected_indiv)
            
            LONG = 75
            FAC = len(self.selected_indiv)/3

            print "\n╔{0}╗".format('═'*(LONG-4))
            print "║{0}║".format("SELECTED".center(LONG-4))
            print "╟{0}┬{0}┬{0}╢".format('─'*((LONG-4)/3))
            for i in range(FAC):
                print "║{0:7.2f} ┆ {1:<13}│{2:7.2f} ┆ {3:<13}│{4:7.2f} ┆ {5:<13}║".format(self.selected_score[i],select[i],self.selected_score[i+FAC],select[i+FAC],self.selected_score[i+FAC*2],select[i+FAC*2])
            if NB_INDIV%3==1:
                print "║{0:7.2f} ┆ {1:<13}│{2:<7} ┆ {3:<13}│{4:<7} ┆ {5:<13}║".format(self.selected_score[-1],select[-1],'','','','')
            elif NB_INDIV%3==2:
                print "║{0:7.2f} ┆ {1:<13}│{2:7.2f} ┆ {3:<13}│{4:<7} ┆ {5:<13}║".format(self.selected_score[-1],select[-1],self.selected_score[-2],select[-2],'','')
            print "╚{0}╧{0}╧{0}╝".format('═'*((LONG-4)/3))

        # INFO GENERATION BOX
        if INFO_GEN and self.generation%INFO_FREQ==0:
            print "\n╔{0}╗\n║ {2}{1}{3} ║".format('═'*39,"GENERATION {0}".format(self.generation).center(37),HEADER,ENDC)
            print "║ {2}{0:7.2f} ┆ {1:<10} : {4:<14}{3} ║".format(min(self.selected_score),"Best Score",OKGREEN,ENDC,self.best_ever_indiv[0].id)
            print "║ {2}{0:7.2f} ┆ {1:<27}{3} ║".format(float(sum(self.selected_score))/len(self.selected_score), "Mean Score",OKBLUE,ENDC)
            print "║ {2}{0:7.2f} ┆ {1:<27}{3} ║".format(max(self.selected_score), "Worst Score",FAIL,ENDC)
            print "╚{0}╝".format('═'*39)
            
            if self.generation%PLOT_GR==0 or (self.generation==1 and PLOT_GEN_ZERO):
                i=1
                for indi in self.indiv :
                    i += 100./len(self.indiv)
                    indi.graphizer("Generation {}".format(self.generation),i)
        
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
        nx.draw(self.graph)
        b=plt.savefig("img/"+str(label)+"_"+str(self.id)+".png") # save as png
        plt.clf()

    def degree_graph(self,generation,i):
        plt.plot(self.list_degrees,self.list_count_log)
        plt.savefig("img/DG_gen"+str(generation)+"_id"+str(i)+"_graph"+str(self.id)+".png") # save as png
        plt.clf()
            
    def clique_graph(self,generation,i):
        plt.plot(self.list_degrees,self.list_meanclust_log)
        plt.savefig("img/CG_gen"+str(generation)+"_id"+str(i)+"_graph"+str(self.id)+".png") # save as png
        plt.clf()
                
    def apply_mutations(self):
        m = self.graph_to_adj_mat()
        for x in range(int(0.3*NB_NODES)):
            i = rd.randint(0,NB_NODES-1)
            j = rd.randint(i,NB_NODES-1)
            m = mutations.substitution(m,i,j)
        return m

    def power_degree_law(self,generation,i):
        """ power degree law """
        if generation%PLOT_PDL==0 :
            self.degree_graph(generation,i) # Plot
        
        slope=stats.linregress(self.list_degrees,self.list_count_log)
        SCE=(slope[4]**2)*NB_NODES ## 25 ? pour NB_NODES ?

        if slope[0] > 0 : self.penalite += 20
            
        self.score_pdl = abs(-1.5-slope[0])*10+SCE

        
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
        if generation%PLOT_CF==0 :
            self.clique_graph(generation,i) # Plot
        
        slope=stats.linregress(self.list_degrees,self.list_meanclust_log)
        SCE=(slope[4]**2)*NB_NODES

        if slope[0] > 0 : self.penalite += 20 
        self.score_cf = abs(-2-slope[0])*10+SCE
        
        if generation%1001==0:
            #print ("\n" + str(self.id))
            #print ("id="+str(i))
            print "liste des degrés :"
            print self.list_degrees
            print "list des meanclust  : "
            print self.list_meanclust
            #print slope
            #print "erreur de pente : "+str(abs(-9-slope[0])*10)
            #print "SCE : " + str(SCE)
            #print "score :",self.score_pdl
        
    def small_world(self):
        """ Compute small world score of graph """
        L = nx.average_shortest_path_length(self.graph)
        C = nx.average_clustering(self.graph)
        self.score_sw = (1-C)*L # A préciser !

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

        # Build clustfk (clustering~k)
        self.clustfk = {} # {k:mean(local_clust_of_k_deg_nodes)}
        for x in self.list_degrees : self.clustfk[x]=[]
        zipp = zip(self.deg_dict.values(), self.deg_dict.keys()) 
        map(lambda x: self.clustfk[x[0]].append(self.clustl_dict[x[1]]), zipp)
        for x in self.clustfk.keys() : self.clustfk[x]=np.mean(self.clustfk[x])

        # Lists
        values = sorted((self.deg_dict.values()))
        self.list_count = [values.count(x) for x in self.list_degrees]
        self.list_meanclust = self.clustfk.values()

        # Delete all 0 (for log) (deprecated)
        """
        self.list_degrees_cf = copy.deepcopy(self.list_degrees)
        self.list_degrees_dg = copy.deepcopy(self.list_degrees)
        
        indices_zeros = list(np.where(np.array(self.list_meanclust) == 0)[0])
        if len(indices_zeros) != 0:
            [self.list_meanclust.pop(i) for i in sorted(indices_zeros,reverse=True)]
            [self.list_degrees_cf.pop(i) for i in sorted(indices_zeros,reverse=True)] 

        indices_zeros = list(np.where(np.array(self.list_degrees) == 0)[0])
        if len(indices_zeros) != 0:
            [self.list_meanclust.pop(i) for i in sorted(indices_zeros,reverse=True)]
            [self.list_degrees_cf.pop(i) for i in sorted(indices_zeros,reverse=True)]
            [self.list_degrees_dg.pop(i) for i in sorted(indices_zeros,reverse=True)]
            [self.list_count.pop(i) for i in sorted(indices_zeros,reverse=True)]
        """
        # Log
        self.list_degrees_log = [math.log10(x+EPS) for x in self.list_degrees]
        self.list_count_log = [math.log10(x+EPS) for x in self.list_count]
        self.list_meanclust_log = [math.log10(x+EPS) for x in self.list_meanclust]
        
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

####################################################################
#			Main
####################################################################
if __name__ == "__main__":
    a = Population()
    a.genetic_algo()
