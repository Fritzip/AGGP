#! /usr/bin/python
# -*- coding: utf-8 -*-

import threading, time

from globals import *

####################################################################
#			Progress Bar
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
            if not self.stopped :
                self.remaining = self.time_last-(time.time()-self.zero)
                i = int(100-(100.*self.remaining/self.time_last))
                update_progress("Generation "+str(self.gen),i)
                time.sleep(0.05)
            else:
                break
   
