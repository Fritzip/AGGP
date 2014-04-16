#! /usr/bin/python
# -*- coding: utf-8 -*-
# Dependancies : networkx, numpy, graphviz, matplotlib

from individual import *
from progressbar import *

####################################################################
#			Population
####################################################################

class Population():
    """ Population including all the individuals """
    def __init__(self,size_pop=NB_INDIV,size_indiv=NB_NODES):
        if VERBOSE : print "Initialization ..." 
        self.indiv = []
        self.score = []
        self.score_pdl = []
        self.score_sw = []
        self.score_cf = []

        self.best_pdl = []
        self.best_sw = []
        self.best_cf = []
        
        self.selected_indiv = []
        self.selected_score = []
        
        self.best_last_gen_score = []
        self.best_last_gen_indiv = []
        
        self.best_ever_score = []
        self.best_ever_indiv = []
        
        self.next_gen = []

        self.matrice_3D = []

        self.generation = 1
        
        self.size_pop = size_pop
        self.size_indiv = size_indiv
        self.nb_best = int(self.size_pop*RATE_ELITISM)+1
        
        mean_shortest_path=[]
        mean_coefficient_clustering=[]
        
        for i in range(self.size_pop):
            self.indiv.append(Individual(nb_nodes=self.size_indiv,id=rd.choice(NAMES)))
            self.score.append(0)
            self.score_pdl.append(0)
            self.score_sw.append(0)
            self.score_cf.append(0)
            self.indiv[i].score_toolbox()

    
    def genetic_algo(self):
        start_algo = time.time()
        time_laps = 0
        update_gen(1)
        
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
                elif INFO_FREQ > 1 or QUIET:
                    update_gen(self.generation)
                start = time.time()
                self.evaluation()
                self.selection()
                self.crossormut()
                self.save() # in files
                time_laps = time.time()-start

                if PROGRESS_GEN and self.generation != 0:
                    bar.stop()
                    update_progress("Generation {0:3d}/{1}".format(self.generation,NB_GEN),100)

                if CYTO3D:
                    matrix_adjacency = self.best_ever_indiv[0].graph_to_adj_mat()
                    matrix_adjacency = np.array(matrix_adjacency)				
                    self.matrice_3D.append(matrix_adjacency)
                
                self.prints()
                self.plots()

                self.indiv = []
                self.indiv = copy.deepcopy(self.next_gen)
                self.next_gen = []
                
                self.generation += 1

            except KeyboardInterrupt:
                if PROGRESS_GEN and self.generation != 0:
                    bar.stop()
                print "\n"+WARNING+"Keyboard Interrupt"+ENDC
                break
        
        # Close all files
        self.fscore.close()
        self.fpdl.close()
        self.fsw.close()
        self.fcf.close()

        if CYTO3D:
            self.convert_xgmml()

        print "\n{1}Done in {0:.2f} sec {2}".format((time.time()-start_algo),OKGREEN,ENDC)

        if self.generation > 1 and SAVE:
            while True:
                try:
                    n = input("Sauvegarde des n [0:{}] meilleurs individus : n = ".format(self.nb_best))
                    if n == "":
                        n = 1
                    elif n > self.nb_best:
                        raise IOError
                    break 
                except:
                    print FAIL+"Wrong answer"+ENDC
                    pass
            
            for i in range(n):
                self.print_info_indiv(self.best_ever_indiv[i])
            for i in range(n):
                self.best_ever_indiv[i].graphizer("Best"+str(i),(i+1)*100./n)
                self.save2sif(self.best_ever_indiv[i])
            self.fitness3D(self.best_pdl,self.best_sw,self.best_cf)
            os.system('gnuplot ../out/scores.gp')
   
                #self.selected_indiv[i].degree_graph("PDL Graphs Generation {}".format(self.generation),i)
                #self.selected_indiv[i].clique_graph("Clique Graphs Generation {}".format(self.generation),i)
        
        print ""
            
    def evaluation(self):
        """ Compute the score using fitness function """
        count = 0
        for i in range(len(self.indiv)):
            self.score[i] = self.indiv[i].calc_score()
            self.score_pdl[i] = self.indiv[i].score_pdl*PDL
            self.score_sw[i] = self.indiv[i].score_sw*SW
            self.score_cf[i] = self.indiv[i].score_cf*CF
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
                ERROR = True
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
                self.cross(tuple(sample),mutations=True)

            #Mutations
            else :
                sample = rd.sample(self.pick_indiv,1)[0]
                self.trans_indiv(sample)
                self.mutation(sample)

        
# Le modèle de mutation a changé. Tous les individus subissent des mutations,
# appliquées aléatoirement avec les probabilités choisies. Ceci a pour but de répartir
# la variabilité dans la population et de ne pas fausser notre image des taux
# en utilisant un produit de probabilités (RATE_MUT * 0.3 dans la version précédente)

# Moins d'étapes => Plus de contrôle

# Du coup, p'tite modif de la fonction cross pour qu'elle inclue des mutations.


#             # Mutation 
#             elif rand < RATE_MUT+RATE_CROSS:
#                 sample = rd.sample(self.pick_indiv,1)[0]
#                 self.trans_indiv(sample)
#                 self.mutation(sample)

#             # Nothing 
#             else :
#                 sample = rd.sample(self.pick_indiv,1)[0]
#                 self.trans_indiv(sample)
#                 self.next_gen.append(sample)
                   
    def trans_indiv(self,x):
        self.pick_indiv.remove(x)
        
    def cross(self,(ind1,ind2),mutations=True):
        """ Simulate the reproduction beetween two individuals : ~70% of variability """
        if mutations :
            A = ind1.apply_mutations()
            B = ind2.apply_mutations()
        else :
            A = ind1.graph_to_adj_mat()
            B = ind2.graph_to_adj_mat()

        X = np.array(nx.adjacency_matrix(nx.fast_gnp_random_graph(self.size_indiv,0.2))).astype(int)
        self.next_gen.append(Individual(mat=((A+B)%2)*X+(A+B)/2,id=rd.choice(NAMES)))
        self.next_gen.append(Individual(mat=((A+B)%2)*((X+np.ones((self.size_indiv,self.size_indiv),dtype=np.int))%2)+(A+B)/2,id=rd.choice(NAMES)))
        
    def mutation(self,ind):
        """~20% of variability """
        self.next_gen.append(Individual(mat = ind.apply_mutations(), id=ind.id))#rd.choice(NAMES)))

    def save(self):
        best = self.best_ever_indiv[0]
        self.fscore.write(str(self.generation)+'\t'+str(max(self.score))+'\t'+str(sum(self.score)/len(self.score))+'\t'+str(min(self.score))+'\t'+str(best.score)+'\n')
        self.fpdl.write(str(self.generation)+'\t'+str(max(self.score_pdl))+'\t'+str(sum(self.score_pdl)/len(self.score_pdl))+'\t'+str(min(self.score_pdl))+'\t'+str(best.score_pdl)+'\n')
        self.fsw.write(str(self.generation)+'\t'+str(max(self.score_sw))+'\t'+str(sum(self.score_sw)/len(self.score_sw))+'\t'+str(min(self.score_sw))+'\t'+str(best.score_sw)+'\n')
        self.fcf.write(str(self.generation)+'\t'+str(max(self.score_cf))+'\t'+str(sum(self.score_cf)/len(self.score_cf))+'\t'+str(min(self.score_cf))+'\t'+str(best.score_cf)+'\n')
        self.best_pdl.append(min(self.score_pdl))
        self.best_sw.append(min(self.score_sw))
        self.best_cf.append(min(self.score_cf))


    def save2sif(self,indiv):
        m = indiv.graph_to_adj_mat()
        sif = open(OUT+indiv.id+'.sif','w')
        for j in range(NB_NODES):
            for i in range(j,NB_NODES):
                if m[i][j]==1:
                    sif.write(str(i)+'\t'+str(j)+'\n')
        sif.close()

    def convert_xgmml(self):            # conversion de la matrice 3D d'adjacence en format XGMML
        
        f = open(OUT+"AGGP.xgmml","w")
        
        # En tête :
        
        f.write('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n')
        f.write('<!-- Created by 4BiM -->\n',)
        f.write('<graph label=" BANANA "\n')
        f.write('    directed="1">\n')
        f.write('    <graphics fill="#000000"/>\n')
        
        # Parcours de la matrice 3D :
        
        # Nodes :
        
        for nodes in range(len(self.matrice_3D[0][0])) :
            sentence = '  <node label="node_' + str(nodes) + '" id="' + str(nodes) + '" start="0" end="' + str(len(self.matrice_3D)+1) + '">\n'
            f.write(sentence)
            sentence = '    <graphics type="CIRCLE" size="8" fill="#0000FF"/>\n'
            f.write(sentence)
            f.write('  </node>\n');
                
        # Edges :

        for j in range(len(self.matrice_3D[0])):
            for i in range(len(self.matrice_3D[0][0])):
                start = 1
                sentence = ''
                sentence_end =  False
                for k in range(len(self.matrice_3D)):
                    if self.matrice_3D[k][j][i] != 0.0 and k < NB_GEN-1 :
                        source = j
                        target = i
                        weight = self.matrice_3D[k][j][i]
                        distance  = 2#1/weight
                        sentence = '  <edge label="edge_' + str(source) + '_' + str(target) + '_' + str(start) + '" source="' + str(source) + '" target="' + str(target) + '" start="' + str(start) + '" end="' + str(k+2) + '">\n'				
                        sentence += '    <graphics width="1" fill="#FFFFFF"/>\n'
                    else : 
                        if sentence != '' :
                            sentence += '  </edge>\n'
                        f.write(sentence)
                        sentence = ''
                        start = k+2						
        f.write('</graph>\n')
        f.close()

    def fitness3D(self,x,y,z):
        mpl.rcParams['legend.fontsize'] = 10

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.invert_yaxis()

        ax.plot(x, y, z, '-o', label='fitness3D')
        ax.legend()
    
        plt.show()
        
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
                print "║{0:8.2f} ┆ {1:<12}║{2:8.2f} ┆ {3:<12}║".format(self.best_ever_score[i],ever[i],self.best_last_gen_score[i],last[i])

            print "╚{0}╩{0}╝".format('═'*((LONG-3)/2))
        
        # INFO SELECTED INDIV
        if INFO_SELECT and self.generation%INFO_FREQ==0:
            select = map(lambda x: x.id ,self.selected_indiv)
            
            LONG = 75
            FAC = len(self.selected_indiv)/3

            print "\n╔{0}╗".format('═'*(LONG-4))
            print "║{0}║".format("SELECTED".center(LONG-4))
            print "╟{0}┬{0}┬{0}╢".format('─'*((LONG-4)/3))
            for i in range(FAC):
                print "║{0:8.2f} ┆ {1:<12}│{2:8.2f} ┆ {3:<12}│{4:8.2f} ┆ {5:<12}║".format(self.selected_score[i],select[i],self.selected_score[i+FAC],select[i+FAC],self.selected_score[i+FAC*2],select[i+FAC*2])
            if NB_INDIV%3==1:
                print "║{0:8.2f} ┆ {1:<12}│{2:<8} ┆ {3:<12}│{4:<8} ┆ {5:<12}║".format(self.selected_score[-1],select[-1],'','','','')
            elif NB_INDIV%3==2:
                print "║{0:8.2f} ┆ {1:<12}│{2:8.2f} ┆ {3:<12}│{4:<8} ┆ {5:<12}║".format(self.selected_score[-1],select[-1],self.selected_score[-2],select[-2],'','')
            print "╚{0}╧{0}╧{0}╝".format('═'*((LONG-4)/3))

        # INFO GENERATION BOX
        if INFO_GEN and self.generation%INFO_FREQ==0:
            indi = self.best_ever_indiv[0]
            print "\n╔{0}╗\n║ {2}{1}{3} ║".format('═'*39,"GENERATION {0}".format(self.generation).center(37),HEADER,ENDC)
            print "║ {2}{0:8.2f} ┆ {1:<9} : {4:<13}{3} ║".format(min(self.selected_score),"Best Score",OKGREEN,ENDC,self.best_ever_indiv[0].id)
            print "║ {0:8.2f} ┆ {1:<26} ║".format(indi.score_pdl*PDL,"Power Degree Law")
            print "║ {0:8.2f} ┆ {1:<26} ║".format(indi.score_sw*SW,"Small World")
            print "║ {0:8.2f} ┆ {1:<26} ║".format(indi.score_cf*CF,"Clique Formation")
            print "║ {2}{0:8.2f} ┆ {1:<26}{3} ║".format(float(sum(self.selected_score))/len(self.selected_score), "Mean Score",OKBLUE,ENDC)
            print "║ {2}{0:8.2f} ┆ {1:<26}{3} ║".format(max(self.selected_score), "Worst Score",FAIL,ENDC)
            print "╚{0}╝".format('═'*39)
            
    def print_info_indiv(self,indi):
        #indi = self.indiv[i]
        print "\n+{}+".format('-'*30)
        print "|{}|".format(indi.id.center(30))
        print "| {0:8.2f} | {1:<17} |".format(indi.score_pdl*PDL,"Power Degree Law")
        print "| {0:8.2f} | {1:<17} |".format(indi.score_sw*SW,"Small World")
        print "| {0:8.2f} | {1:<17} |".format(indi.score_cf*CF,"Clique Formation")
        print "| {0:8.2f} | {1:<17} |".format(indi.score,"Global")
        print "+{}+".format('-'*30)

    def plots(self):
        if self.generation%PLOT_GR==0 or (self.generation==1 and PLOT_GEN_ZERO):
            i=1
            for indi in self.indiv :
                i += 100./len(self.indiv)
                indi.graphizer("Indiv Generation {}".format(self.generation),i)
            print ""
            
        if self.generation%PLOT_PDL==0:
            i=1
            for indi in self.indiv :
                i += 100./len(self.indiv)
                indi.degree_graph("PDL Graphs Generation {}".format(self.generation),i)
            print ""

        if self.generation%PLOT_CF==0:
            i=1
            for indi in self.indiv :
                i += 100./len(self.indiv)
                indi.clique_graph("CF Graphs Generation {}".format(self.generation),i)
            print ""
        

