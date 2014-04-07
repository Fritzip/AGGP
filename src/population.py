#! /usr/bin/python
# -*- coding: utf-8 -*-

from individual import *
from progressbar import *

####################################################################
#			Population
####################################################################

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
        self.fscore = open(OUT+'evo_score','w')
        self.fpdl = open(OUT+'evo_pdl','w')
        self.fsw = open(OUT+'evo_sw','w')
        self.fcf = open(OUT+'evo_cf','w')

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
                if PROGRESS_GEN and self.generation != 0:
                    bar.stop()
                print "\nKeyboard Interrupt"
                break
        
        # Close all files
        self.fscore.close()
        self.fpdl.close()
        self.fsw.close()
        self.fcf.close()
        
        print "\nDone in %.3f sec"%(time.time()-start_algo)
        
        while True:
            n = input("Sauvegarde des n meilleurs individus (defaut n=1).\nn [0:{}]= ".format(NB_INDIV))
            if n == "":
                n = 1
                break
            try:
                n = int(n)
                break
            except:
                print "int demandé"
                pass
        for i in range(n):
            self.print_info_indiv(self.selected_indiv[i])
        for i in range(n):
            self.selected_indiv[i].graphizer("Best",(i+1)*100./n)
            self.save2sif(self.selected_indiv[i])
        print "\n"
        

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
                self.print_info_indiv(self.indiv[i])
            

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

            # Crossing
            if rand < RATE_CROSS and len(self.pick_indiv) > 1: 
                sample = rd.sample(self.pick_indiv,2)
                map(self.trans_indiv ,sample)
                self.cross(tuple(sample))

            # Mutation 
            elif rand < RATE_MUT+RATE_CROSS:
                sample = rd.sample(self.pick_indiv,1)[0]
                self.trans_indiv(sample)
                self.mutation(sample)

            # Nothing 
            else :
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


    def save2sif(self,indiv):
        m = indiv.graph_to_adj_mat()
        sif = open(OUT+indiv.id+'.sif','w')
        for j in range(NB_NODES):
            for i in range(j,NB_NODES):
                if m[i][j]==1:
                    sif.write(str(i)+'\t'+str(j)+'\n')
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

    def print_info_indiv(self,indi):
        #indi = self.indiv[i]
        print "+{}+".format('-'*30)
        print "|{}|".format(indi.id.center(30))
        print "| {0:7.2f} | {1:<18} |".format(indi.score_pdl,"Power Degree Law")
        print "| {0:7.2f} | {1:<18} |".format(indi.score_sw,"Small World")
        print "| {0:7.2f} | {1:<18} |".format(indi.score_cf,"Clique Formation")
        print "| {0:7.2f} | {1:<18} |".format(indi.score,"Global")
        print "+{}+\n".format('-'*30)
