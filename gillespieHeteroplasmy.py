# -*- coding: utf-8 -*-
# Uses the Gillespie algorithm to simulate mitochondrial heteroplasmy dynamics in a population of cells
# Author: Tom McGrath
# Date started: 21/01/14

import numpy as np
import scipy.misc
import matplotlib.pyplot as plt
import collections as coll

# define constants
mitoGenRate = 1 # needs changing!
divisionRate = 1
targetNumA = 100
cellCycleTime = 3000
divisionHelper = {} # stores results of factorial computation ((numA, numB), (numALost, numBLost))

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
    
def choose(a,b):
    returnVar = scipy.misc.factorial(a)/(scipy.misc.factorial(b)*scipy.misc.factorial((a-b)))
    return int(returnVar)
    
def totalChoices(a):
    returnVar = 0
    for i in range(0, a+1):
        returnVar += choose(a,i)
    return int(returnVar)

# transition functions
def probAPlusOne((numA, numB)):
    heteroplasmy = float(numA)/(float(numA + numB))
    transitionProb = mitoGenRate * heteroplasmy * population[(numA, numB)][0]    
    return transitionProb
    
def probBPlusOne((numA, numB)):
    heteroplasmy = float(numA)/(float(numA + numB))
    transitionProb = mitoGenRate * (1 - heteroplasmy) * population[(numA, numB)][0]   
    return transitionProb

def splittingProb((numA, numB), (numALost, numBLost), divisionHelper):
    # can this be done with dynamic programming to speed up calculations of factorials?
    if ((numA, numB), (numALost, numBLost)) in divisionHelper.keys():
        return divisionHelper[((numA, numB), (numALost, numBLost))]        
    else:
        returnVar = float(choose(numA, numALost)*choose(numB, numBLost))/float(totalChoices(numA)*totalChoices(numB))
        divisionHelper[((numA, numB), (numALost, numBLost))] = returnVar
        return returnVar
    
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
        
            