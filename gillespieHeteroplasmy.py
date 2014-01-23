# -*- coding: utf-8 -*-
# Uses the Gillespie algorithm to simulate mitochondrial heteroplasmy dynamics in a population of cells
# Author: Tom McGrath
# Date started: 21/01/14

import numpy as np
import matplotlib.pyplot as plt
import collections as coll

# define constants
mitoGenRate = 1e-3 # needs changing!
divisionRate = 1
targetNumA = 100
cellCycleTime = 3000

# initialise population dictionary (typeA, typeB):(number, [cell cycle timers as a deque])
population = {}

# setup initial population
population[(100,100)] = [1, coll.deque([2500])]
population[(150,100)] = [1, coll.deque([1000])]
population[(100,200)] = [1, coll.deque([2000])]

# return shortest time to deterministic division event
def getShortestClock(population):
    shortestClock = cellCycleTime # this is the worst it can be    
    for x in population.keys():
        try:        
            if population[x][1][0] < shortestClock:
                shortestClock = population[x][1][0]
                populationMember = x
            else:
                continue
        except IndexError:
            continue # the population dictionary entry is empty
    return (shortestClock, populationMember)

# transition functions
def probAPlusOne((numA, numB)):
    if numA > targetNumA:
        return 0
    else:
        heteroplasmy = float(numA)/(float(numA + numB))
        transitionProb = mitoGenRate * heteroplasmy * population[(numA, numB)][0]
        if population[(numA, numB)][0] == 0:
            transitionProb = 0
        return transitionProb
    
def probBPlusOne((numA, numB)):
    if numA > targetNumA:
        return 0
    else:
        heteroplasmy = float(numA)/(float(numA + numB))
        transitionProb = mitoGenRate * (1 - heteroplasmy) * population[(numA, numB)][0]
        if population[(numA, numB)][0] == 0:
            transitionProb = 0
        return transitionProb
        
def chooseDivision((numA, numB)):
    numALost = np.random.binomial(numA, 0.5)
    numBLost = np.random.binomial(numB, 0.5)    
    if numALost == numA:
        numALost -= 1
    if numBLost == numB:
        numBLost -= 1
    return ((numA-numALost, numB-numBLost),(numALost, numBLost))
    
def createCell((numA, numB)):
    # creates a new entry in population[(numA, numB)] with timer cellCycleLength plus some random number in [0,1)
    try:    
        population[(numA, numB)][0] += 1
        population[(numA, numB)][1].append(cellCycleTime + np.random.uniform()) # init with some random delta to prevent collisions
    except KeyError:
        population[(numA, numB)] = [1]
        population[(numA, numB)].append(coll.deque([cellCycleTime + np.random.uniform()]))    
        
def advanceAllClocks(tau):
    # steps all clocks in the population clock deques forward by tau
    for key in population.keys():
        for i in range(0, len(population[key][1])):
            population[key][1][i] -= tau
    
# the Gillespie loop:
def run(runTime):
    t = 0
    while t < runTime:
        a0 = calcA0(population)
        tau = np.random.exponential(float(1/a0))
        print 'tau ', tau
        
        if tau > getShortestClock(population)[0]:
            # partition cell from shortest clock cell group
            (numA, numB) = getShortestClock(population)[1]
            ((newA1, newB1), (newA2, newB2)) = chooseDivision((numA, numB))
            # create new cells with appropriate mito numbers (initialise cell clocks w/some randomness to prevent collision)
            createCell((newA1, newB1))
            createCell((newA2, newB2))
            print 'time ', t
            print 'dividing ', (numA, numB)
            print 'creating ', (newA1, newB1)
            print 'creating ', (newA2, newB2)
            # delete old cell
            population[(numA, numB)][0] -= 1
            timeAdvance = population[(numA, numB)][1].popleft()
            if population[(numA, numB)][0] == 0:
                del population[(numA, numB)]
            # advance time by shortest clock time
            t += timeAdvance
            advanceAllClocks(timeAdvance)
            
        elif tau <= getShortestClock(population)[0]: # THIS IS THE BIT THAT NEEDS WORK
            # partition [0, 1) appropriately by a_i
            # draw from [0, 1)
            # do the appropriate action
            target = np.random.uniform()
            ai = 0
            for x in population.keys():
                if (ai + probAPlusOne(x) > target):
                    print 'time ', t                    
                    print 'adding type A mitochondria to ', x
                    createCell((x[0]+1, x[1])) # PROBLEM - needs to inherit cycle timer!
                    population[x][0] -= 1
                    population[x][1].popleft()
                    if population[x][0] == 0:
                        del population[x]
                    continue
                else:
                    ai += probAPlusOne(x)/a0
                    
                if (ai + probBPlusOne(x) > target):
                    print 'time ', t
                    print 'adding type B mitochondria to ', x
                    createCell((x[0], x[1]+1))
                    population[x][0] -= 1
                    population[x][1].popleft()
                    if population[x][0] == 0:
                        del population[x]
                    continue
                else:
                    ai += probBPlusOne(x)/a0
                    
            t += tau
            advanceAllClocks(tau)
        
    
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
    return returnVar
    
def getPop(popDict = population):
    for x in popDict:
        print x, popDict[x]
        
            