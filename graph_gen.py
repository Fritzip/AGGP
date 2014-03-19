#! /usr/bin/python
# -*- coding: utf-8 -*-

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random as rd
import copy

####################################################################
#			Global Functions
####################################################################
def symetrize(a):
		return -a*a.T + a + a.T

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
	def __init__(self,size_pop=20,size_indiv=25):
		self.indiv = []
		self.score = []
		
		self.selected_indiv = []
		self.selected_score = []
		
		self.best_last_gen_score = []
		self.best_last_gen_indiv = []
		
		self.best_ever_score = []
		self.best_ever_indiv = []
		
		self.next_gen = []

		self.size_pop = size_pop
		self.size_indiv = size_indiv
		self.nb_best = int(self.size_pop*0.2)+1

		for i in range(self.size_pop):
			self.indiv.append(Individual(nb_nodes=self.size_indiv,id=str(i)))
			self.score.append(0)

	def genetic_algo(self):
		""" Inspired by
		http://fr.wikipedia.org/wiki/Algorithme_g%C3%A9n%C3%A9tique#Sch.C3.A9ma_r.C3.A9capitulatif
		"""
		i=0
		while True:
			self.evaluation()
			self.selection()
			self.crossormut()
			i += 1
			print "génération = ",i
			self.indiv = []
			self.indiv = copy.deepcopy(self.next_gen)
			self.next_gen = []

	def evaluation(self):
		""" Compute the score using fitness function """
		for i in range(len(self.indiv)):
			self.score[i] = self.indiv[i].calc_score()
		print "score =",self.score		

	def roulette_wheel_selec(self):
		""" Roulette wheel selection """
		zipper = zip(self.score,self.indiv)
		generator = weighted_sample(zipper,self.size_pop/2)
		for i in generator:
			self.selected_indiv.append(i)
		print "selected =",self.selected_indiv
		print "len(selected) =",len(self.selected_indiv)
		
	def elitist_selec(self):
		""" Elitist (20%) and tournament selection (80%) """
		self.best_last_gen_score, self.best_last_gen_indiv = map(lambda x : list(x[0:self.nb_best]),zip(*(sorted(zip(self.score,self.indiv)))[::-1]))
		print "best last",self.best_last_gen_score,self.best_last_gen_indiv
		self.best_ever_score, self.best_ever_indiv = map(lambda x : list(x[0:self.nb_best]),zip(*(sorted(zip(self.best_last_gen_score+self.best_ever_score,self.best_last_gen_indiv+self.best_ever_indiv)))[::-1]))
		print "best ever",self.best_ever_score,map(lambda x: x.id ,self.best_ever_indiv)
		self.selected_indiv += self.best_ever_indiv
		self.selected_score += self.best_ever_score
		while len(self.selected_indiv) != self.size_pop:
			a,b = np.random.choice(range(len(self.indiv)),2)
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
		if rand < 0.7: return higher
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
				print "rien"
				sample = rd.sample(self.selected_indiv,1)[0]
				self.trans_indiv(sample)
				self.next_gen.append(sample)

	def trans_indiv(self,x):
		#self.next_gen.append(x)
		self.selected_indiv.remove(x)
		
	def cross(self,(ind1,ind2)):
		print "croisement"
		""" Simulate the reproduction beetween two individuals : 70% of variability """
		A = ind1.graph_to_adj_mat()
		B = ind2.graph_to_adj_mat()
		X = symetrize(np.random.binomial(1, 0.2, self.size_indiv**2).reshape(self.size_indiv,self.size_indiv))
		self.next_gen.append(Individual(mat=((A+B)%2)*X+(A+B)/2,id=ind1.id))
		self.next_gen.append(Individual(mat=((A+B)%2)*((X+np.ones((self.size_indiv,self.size_indiv),dtype=np.int))%2)+(A+B)/2,id=ind2.id))

	def mutation(self,ind):
		print "mutation"
		""" 30% of variability """
		self.next_gen.append(Individual(mat = ind.apply_mutations(), id=ind.id))

####################################################################
#			Individual
####################################################################
class Individual():
	""" Metabolic graph """
	def __init__(self,mat=0, nb_nodes=25,id=0):
		#self.matrix = symetrize(np.random.binomial(1, 0.2, nb_nodes**2).reshape(nb_nodes,nb_nodes)) # random
		self.id = id
		self.generation = 0
		self.score = 0
		self.graph = self.adj_mat_to_graph(mat) if isinstance(mat,np.ndarray) else nx.fast_gnp_random_graph(nb_nodes,0.2)# random graph
		
	def graph_to_adj_mat(self):
		return nx.adjacency_matrix(self.graph)
	
	def adj_mat_to_graph(self, mat):
		return nx.from_numpy_matrix(mat)
		
	def graphizer(self):
		nx.draw(self.graph)
		plt.savefig("graph"+str(self.id)+"_gen"+str(self.generation)+".png") # save as png
		plt.clf()
		#plt.show() # display
		
	def apply_mutations(self):
		return self.graph_to_adj_mat()
		
	def calc_score(self):
		""" Fitness function """
		self.score += 1 
		self.generation += 1
		return self.score

####################################################################
#			Main
####################################################################
if __name__ == "__main__":
	a = Population()
	a.genetic_algo()
