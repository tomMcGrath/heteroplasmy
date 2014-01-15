# Stochastic simulation of heteroplasmy dynamics in mitochondrial populations
# Tom McGrath, 15/01/14
#
# Assume a population of N cells, each of which can have mitochondria of varying type
# Cells can produce new mitochondria, or divide and split mitochondria by some rule
# Cells 'want' to achieve a certain level of 'good' mitochondria
# However they can only influence total mitochondrial population

# This version:
# - hard-code selection rules
# - assume only types A (good) and B (bad)
# - cells don't die even if fully B-heteroplasmic

# Libraries
import numpy as np

# Base cell class
class cell(object):
    nextID = 0
    def __init__(self, probTypeA, initialMitoNumber):
        # Initialise unique cell ID
        self.ID = cell.nextID
        cell.nextID += 1
        # Store probability of type A
        self.probTypeA = probTypeA
        # Initialise type A & B populations
        self.numTypeA = 0
        self.numTypeB = 0
        # Grow type A & B mitochondria up to initialMitoNumber
        while (self.numTypeA + self.numTypeB) < initialMitoNumber:
            self.addMito()

    # Mitochondrial replication code
    def addMito(self):
        randomNumber = np.random.random_sample()        
        if randomNumber <= self.probTypeA:
            self.numTypeA += 1
        elif randomNumber > self.probTypeA:
            self.numTypeB += 1
        else:
            print "Unexpected random value"
            raise ValueError

    # Cell division code
    def cellDivision(self):
        print 'not implemented yet!'

    # Interfaces
    def getNumTypeA(self):
        return self.numTypeA

    def getNumTypeB(self):
        return self.numTypeB
    
    def getID(self):
        return self.ID

def test():
    cellList = []
    for i in range(0, 10):
        cellList.append(cell(0.5, 5))

    print 'testing initialisation:'

    for x in cellList:
        print x.getID(), x.getNumTypeA(), x.getNumTypeB(), x.getNumTypeA() + x.getNumTypeB()

    print 'testing adding mitochondria, 1 per cell'
    for x in cellList:
        x.addMito()

    print 'ID, typeA, typeB, total'
    for x in cellList:
        print x.getID(), x.getNumTypeA(), x.getNumTypeB(), x.getNumTypeA() + x.getNumTypeB()




















	
