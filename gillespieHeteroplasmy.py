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
cellCycleTime = 3000

# initialise population dictionary (typeA, typeB):(number, [cell cycle timers as a deque])
population = {}

# setup initial population
population[(1,1)] = (1, coll.deque([0]))

# return shortest time to deterministic division event
def getShortestClock(population):
    shortestClock = cellCycleTime # this is the worst it can be    
    for x in population.keys():
        if population[x][1][0] < shortestClock:
            shortestClock = population[x][1][0]
        else:
            continue
    return shortestClock

# transition functions
def probAPlusOne((numA, numB)):
    heteroplasmy = float(numA)/(float(numA + numB))
    transitionProb = mitoGenRate * heteroplasmy * population[(numA, numB)][0]    
    return transitionProb
    
def probBPlusOne((numA, numB)):
    heteroplasmy = float(numA)/(float(numA + numB))
    transitionProb = mitoGenRate * (1 - heteroplasmy) * population[(numA, numB)][0]   
    return transitionProb
    
def chooseDivisiom((numA, numB)):
    numALost = np.random.binomial(numA, 0.5)
    numBLost = np.random.binomial(numB, 0.5)
    
    return ((numA-numALost, numB-numBLost),(numALost, numBLost))
    
# calculate a0
def calcA0(population):
    # loop over each element of the population dictionary and sum the total probabilities
    a0 = 0.0    
    for x in population.keys():
        #print x #useful to see the workings here
        a0 += probAPlusOne(x)
        #print a0
        a0 += probBPlusOne(x)
        #print a0
    return a0
        
            