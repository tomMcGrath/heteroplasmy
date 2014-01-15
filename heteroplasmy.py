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
    def __init__(self, probTypeA, mitosToAdd, numTypeA, numTypeB):
        # Initialise unique cell ID
        self.ID = cell.nextID
        cell.nextID += 1
        # Store probability of type A
        self.probTypeA = probTypeA
        # Store initialised variables (in case we need some sort of history)
        self.mitosToAdd = mitosToAdd
        self.numTypeA = numTypeA
        self.numTypeB = numTypeB
        # Initialise type A & B populations
        self.numTypeA = numTypeA
        self.numTypeB = numTypeB
        # Grow type A & B mitochondria up to mitosToAdd
        while (self.numTypeA + self.numTypeB) < mitosToAdd:
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
        if self.numTypeA > 0:
            numToLoseA = np.random.binomial(self.numTypeA, 0.5)
        else:
            numToLoseA = 0

        if self.numTypeB > 0:
            numToLoseB = np.random.binomial(self.numTypeB, 0.5)
        else:
            numToLoseB = 0

        return [numToLoseA, numToLoseB] # tell the petri dish what the new cell has

    # Interfaces
    def getNumTypeA(self):
        return self.numTypeA

    def getNumTypeB(self):
        return self.numTypeB
    
    def getID(self):
        return self.ID

class petriDish(object):
    def __init__(self, numCells, probTypeA, cellMitosToAdd):
        # Store cells in this list
        self.cellList = []

        # Store initial variables
        self.numCells = numCells
        self.probTypeA = probTypeA
        self.cellMitosToAdd = cellMitosToAdd

        # Add the initial cell population
        for i in range(0, numCells):
            self.cellList.append(cell(self.probTypeA, self.cellMitosToAdd, 0, 0))

    def newCell(self, probTypeA, mitosToAdd, numTypeA, numTypeB):
        cellList.append(cell(probTypeA, mitosToAdd, numTypeA, numTypeB))

    def makeCellDivide(self, chosenID):
        newCellInfo = cellList[0].cellDivision()
        print newCellInfo

    # Interfaces
    def getNumCells(self):
        return len(self.cellList)

    def printCellList(self): # really just for testing
        print 'ID, typeA, typeB, total'
        for x in self.cellList:
            print x.getID(), x.getNumTypeA(), x.getNumTypeB(), x.getNumTypeA() + x.getNumTypeB()



    
def test():
    petri = petriDish(10, 0.5, 10)
    print petri.getNumCells()
    petri.printCellList()
    print 'testing division'
    petri.makeCellDivide(0)
    petri.printCellList()
    
test()













	
