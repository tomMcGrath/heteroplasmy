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
population[(1,1)] = [1, coll.deque([0])]
population[(2,1)] = [1, coll.deque([1])]
population[(1,2)] = [1, coll.deque([2])]

# return shortest time to deterministic division event
def getShortestClock(population):
    shortestClock = cellCycleTime # this is the worst it can be    
    for x in population.keys():
        if population[x][1][0] < shortestClock:
            shortestClock = population[x][1][0]
            populationMember = x
        else:
            continue
    return (shortestClock, populationMember)

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
    
def createCell((numA, numB)):
    try:    
        population[(numA, numB)][0] += 1
        population[(numA, numB)][1].append(cellCycleTime + np.random.uniform()) # init with some random delta to prevent collisions
    except KeyError:
        population[(numA, numB)] = [1]
        population[(numA, numB)].append(coll.deque([cellCycleTime + np.random.uniform()]))
    
# the Gillespie loop:
def run(runTime):
    t = 0
    while t < runTime:
        a0 = calcA0(population)
        tau = np.random.exponential(float(1/a0))
        
        if tau < getShortestClock(population)[0]:
            # partition cell from shortest clock cell group
            # create new cells with appropriate mito numbers (initialise cell clocks w/some randomness to prevent collision)
            # delete old cell
            # advance time by shortest clock time
            print 'not here yet!'
            
        elif tau >= getShortestClock(population)[0]:
            # partition [0, 1) appropriately by a_i
            # draw from [0, 1)
            # do the appropriate action
            print 'not here yet either!'
        
    
# calculate a0
def calcA0(population):
    # loop over each element of the population dictionary and sum the total probabilities
    returnVar = 0.0    
    for x in population.keys():
        #print x #useful to see the workings here
        returnVar += probAPlusOne(x)
        #print a0
        returnVar += probBPlusOne(x)
        #print a0
    return returnvar
        
            