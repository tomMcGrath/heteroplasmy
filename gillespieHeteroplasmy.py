# -*- coding: utf-8 -*-
# Uses the Gillespie algorithm to simulate mitochondrial heteroplasmy dynamics in a population of cells
# Author: Tom McGrath
# Date started: 21/01/14

import numpy as np
import matplotlib.pyplot as plt

# define constants
mitoGenRate = 0.4 # needs changing!
divisionRate = 1
targetNumA = 100
cellCycleTime = 3000

# initialise population dictionary (typeA, typeB):(number, [cell cycle timers as a deque])
population = {}
# setup initial population
population[(50,1)] = [1, [2500]]

# the Gillespie loop:
def run(runTime):
    resultHolder = []
    timeSeries = []
    popHolder = []
    t = 0
    while t < runTime:
        getPop()
        print 'time: ', t
        a0 = calcA0(population)
        if a0 == 0.0:
            tau = getShortestClock(population)[0]
        else:
            tau = np.random.exponential(float(1/a0))
        print 'tau: ', tau
        
        if tau >= getShortestClock(population)[0]:
            # partition cell from shortest clock cell group
            (numA, numB) = getShortestClock(population)[1]
            ((newA1, newB1), (newA2, newB2)) = chooseDivision((numA, numB))
            # create new cells with appropriate mito numbers (initialise cell clocks w/some randomness to prevent collision)
            createCell((newA1, newB1))
            createCell((newA2, newB2))
            #print 'dividing ', (numA, numB)
            #print 'creating ', (newA1, newB1)
            #print 'creating ', (newA2, newB2)
            # delete old cell
            population[(numA, numB)][0] -= 1
            population[(numA, numB)][1].reverse()
            timeAdvance = population[(numA, numB)][1].pop()
            population[(numA, numB)][1].sort()
            if population[(numA, numB)][0] == 0:
                del population[(numA, numB)]
            # advance time by shortest clock time
            t += timeAdvance
            #print 'time ', t
            advanceAllClocks(timeAdvance)
            
        elif tau < getShortestClock(population)[0]: # THIS IS THE BIT THAT NEEDS WORK
            # partition [0, 1) appropriately by a_i
            # draw from [0, 1)
            # do the appropriate action
            target = np.random.uniform()
            ai = 0
            for x in population.keys():
                if (ai + probAPlusOne(x) > target):                   
                    #print 'adding type A mitochondria to ', x
                    createCell((x[0]+1, x[1]))
                    population[(x[0]+1, x[1])][1].pop() # remove the 'wrong' timer - this is all a bit messy & would be easier without deque
                    population[x][0] -= 1
                    population[x][1].reverse()
                    timer = population[x][1].pop()
                    population[x][1].sort()
                    population[(x[0]+1, x[1])][1].append(timer) # transfers cell clock timer            
                    population[(x[0]+1,x[1])][1].sort()
                    if population[x][0] == 0:
                        del population[x]
                    break
                else:
                    ai += probAPlusOne(x)/a0
                    
                if (ai + probBPlusOne(x) > target):
                    #print 'adding type B mitochondria to ', x
                    createCell((x[0], x[1]+1))
                    population[(x[0], x[1]+1)][1].pop() # remove the 'wrong' timer
                    population[x][0] -= 1
                    population[x][1].reverse()
                    timer = population[x][1].pop()
                    population[x][1].sort()
                    population[(x[0], x[1]+1)][1].append(timer) # transfers cell clock timer            
                    population[(x[0],x[1]+1)][1].sort()
                    if population[x][0] == 0:
                        del population[x]
                    break
                else:
                    ai += probBPlusOne(x)/a0
                    
            t += tau
            advanceAllClocks(tau)
            print t, calcHeteroplasmy(population)
            resultHolder.append(calcHeteroplasmy(population)[0])
            popHolder.append(calcHeteroplasmy(population)[1])
            timeSeries.append(t)
    return timeSeries, resultHolder, popHolder

# return shortest time to deterministic division event
def getShortestClock(population):
    shortestClock = cellCycleTime + 2 # this is the worst it can be    
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
    if numA >= targetNumA:
        return 0        
    elif numA == 0 and numB == 0: # these cells are 'idle' (biologically shouldn't exist!)
        return 0
    else:
        heteroplasmy = float(numA)/(float(numA + numB))
        transitionProb = mitoGenRate * heteroplasmy * population[(numA, numB)][0]
        if population[(numA, numB)][0] == 0:
            transitionProb = 0
        return transitionProb
    
def probBPlusOne((numA, numB)):
    if numA >= targetNumA:
        return 0
    elif numA == 0 and numB == 0:
        return 0
    else:
        heteroplasmy = float(numA)/(float(numA + numB))
        transitionProb = mitoGenRate * (1 - heteroplasmy) * population[(numA, numB)][0]
        if population[(numA, numB)][0] == 0:
            transitionProb = 0
        return transitionProb
        
def chooseDivision((numA, numB)):   
    if numA == 0:
        numALost = 0
    else:
        numALost = np.random.binomial(numA, 0.5)
    if numB == 0:
        numBLost = 0
    else:
        try:
            numBLost = np.random.binomial(numB, 0.5)    
        except ValueError:
            print (numA, numB)
    return ((numA-numALost, numB-numBLost),(numALost, numBLost))
    
def createCell((numA, numB)):
    # creates a new entry in population[(numA, numB)] with timer cellCycleLength plus some random number in [0,1)
    try:    
        population[(numA, numB)][0] += 1
        population[(numA, numB)][1].append(cellCycleTime + np.random.uniform()) # init with some random delta to prevent collisions
    except KeyError:
        population[(numA, numB)] = [1]
        population[(numA, numB)].append([cellCycleTime + np.random.uniform()])
        
def advanceAllClocks(tau):
    # steps all clocks in the population clock deques forward by tau - bugged when advance by float
    for key in population.keys():
        for i in range(0, len(population[key][1])):
            #print population[key][1][i]            
            population[key][1][i] = population[key][1][i] - tau
            
def sortedInsert(dq, val):
    # takes a sorted deque and inserts val
    for i in range(0, len(dq)):
        if val < dq[i]:
            dq.rotate(-(i))
            dq.appendleft(val)
            dq.rotate(i)
            break
    dq.append(val)
        
    
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
    
def calcHeteroplasmy(population):
    totalA = 0
    totalB = 0
    numCells = 0
    for x in population.keys():
        totalA += x[0]*population[x][0]
        totalB += x[1]*population[x][0]
        numCells += population[x][0]
    heteroplasmy = float(totalA)/float(totalA + totalB)
    return (heteroplasmy, numCells)

def graphResults(t, h, p):
    fig, ax1 = plt.subplots()
    ax1.plot(t,h)
    ax1.set_ylim([0,1])
    ax1.set_xlabel('time')
    ax1.set_ylabel('heteroplasmy')
    
    ax2 = ax1.twinx()
    ax2.plot(t, p, color = 'green')
    ax2.set_ylabel('population')
    plt.show()
    

# debug function    
def getPop(popDict = population):
    for x in popDict:
        print x, popDict[x]
