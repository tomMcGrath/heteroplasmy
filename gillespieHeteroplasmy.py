# -*- coding: utf-8 -*-
# Uses the Gillespie algorithm to simulate mitochondrial heteroplasmy dynamics in a population of cells
# Author: Tom McGrath
# Date started: 21/01/14

import numpy as np
import matplotlib.pyplot as plt
import collections as coll

# define constants
mitoGenRate = 1 # needs changing!
divisionRate = 1
targetNumA = 100
cellcycleTime = 3000

# initialise population dictionary (typeA, typeB):(number, [cell cycle timers as a deque])
population = {}

# setup initial population
population[(1,1)] = (1, coll.deque([0]))

# transition functions
def probAPlusOne((numA, numB)):
    heteroplasmy = float(numA)/(float(numA + numB))
    transitionProb = mitoGenRate * heteroplasmy * population[(numA, numB)][0]
    
    return transitionProb
    
def probBPlusOne((numA, numB)_:
    heteroplasmy = float(numA)/(float(numA + numB))
    transitionProb = mitoGenRate * (1 - heteroplasmy) * population[(numA, numB)][0]
    
    return transitionProb

# need to implement division
    
# calculate a0
def calcA0(population):
    # loop over each element of the population dictionary and sum the total probabilities
    
def loop(runTime):
    for t in range(0, runtime):
        
            