# -*- coding: utf-8 -*-
"""
Created on Wed May 26 12:31:35 2021

@author: Thomas Maynadié
"""

from threading import Thread, Lock
from EVS.evs_debug import debug, EVSErrCode

class EVSTask(object):
    def __init__(self, taskID, method, args):
        self.id = taskID
        self.method = method
        self.args = args
        
    def run(self):
        self.method(*self.args)

class EVSWorker(object):
    def __init__(self, _id):
        self.id = _id

        self.task_list = []
        self.exitflag = False
        self.active = False
        self.mutex = Lock()
        
    def thread_routine(self):
        self.active = True # waiting for task
        
        while(len(self.task_list) > 0 and self.exitflag == False):
            self.task_list[0].run()
            self.task_list.pop(0)
                            
        self.active = False # idle
        self.exitflag = False
                
    def activate(self, _activate=True):
        self.mutex.acquire()
        self.active = _activate

        if _activate == True: 
            self.exitflag = False
            
            self.thread = Thread(target=self.thread_routine)
            self.thread.start()
                
        else:
            self.exitflag = True
                
        self.mutex.release()
    
    def addTaskToList(self, task):
        self.task_list.append(task)
        
    def isActive(self):
        self.mutex.acquire()
        ret = self.active
        self.mutex.release()
        
        return ret
    
class EVSThreadManager(object):
    def __init__(self):
      self.workers = []
      self.active = False
      
      self.lastID = 0
      self.taskIterator = 0
      
    def __init__(self, n_workers):
      self.workers = []
      self.active = False

      self.lastID = 0
      self.taskIterator = 0
      
      self.createWorkers(n_workers)

    def createTask(self, threadID, method, args):
        task = EVSTask(self.lastID, method, args)
        self.lastID = self.lastID + 1
        
        self.workers[self.taskIterator].addTaskToList(task)
        
        self.taskIterator = self.taskIterator + 1
        if self.taskIterator == len(self.workers): self.taskIterator = 0
        
    def createWorkers(self, n_workers):
        for i in range(n_workers):
            worker = EVSWorker(i)
            self.workers.append(worker)

    def startSession(self, wait=True):
        if len(self.workers) == 0: debug("no workers in list, Exiting...", EVSErrCode.EVS_ERROR)
        
        else :
            self.active = True
            
            for worker in self.workers:
                if (worker.isActive() == False):
                    worker.activate(True)
            
            if wait == True:
                self.waitForWorkers()
                
            else: 
                debug("not waiting for threads to finish, Exiting...", EVSErrCode.EVS_STATUS)
            
            self.active = False
    
    def waitForWorkers(self):
        flag = True
        while(flag): 
            flag = False
            
            for i in range(len(self.workers)):
                if (self.workers[i].active == True):
                    flag = True
                    break
        
        debug("All threads finished. Exiting...", EVSErrCode.EVS_STATUS)