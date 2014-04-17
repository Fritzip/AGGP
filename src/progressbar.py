#! /usr/bin/python
# -*- coding: utf-8 -*-
# Dependancies : networkx, numpy, graphviz, matplotlib

import threading, time

from globals import *

####################################################################
#			Progress Bar
####################################################################

class Progressbar(threading.Thread):
""" Thread for displaying progress bar"""
    def __init__(self,gen,time_last):
        threading.Thread.__init__(self)
        self.gen = gen
        self.time_last = time_last*0.85
        self.remaining = self.time_last
        self.zero = time.time()
        self.stopped = False
        
    def stop(self):
        self.stopped = True
        
    def run(self):
        while self.remaining > 0:
            if not self.stopped :
                self.remaining = self.time_last-(time.time()-self.zero)
                i = int(100-(100.*self.remaining/self.time_last))
                update_progress("Generation {0:3d}/{1}".format(self.gen,NB_GEN),i)
                time.sleep(0.1)
            else:
                break
   
