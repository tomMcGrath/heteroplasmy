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
import matplotlib as mpl

# Base cell class
class cell(object):
    nextID = 0
    def __init__(self, probTypeA, mitosToAdd, numTypeA, numTypeB, targetMitos):
        # Initialise unique cell ID
        self.ID = cell.nextID
        cell.nextID += 1
        # Store probability of type A
        self.probTypeA = probTypeA
        # Store initialised variables (in case we need some sort of history)
        self.mitosToAdd = mitosToAdd
        self.numTypeA = numTypeA
        self.numTypeB = numTypeB
        self.targetMitos = targetMitos
        # Initialise type A & B populations
        self.numTypeA = numTypeA
        self.numTypeB = numTypeB
        # Grow type A & B mitochondria up to mitosToAdd
        while (self.numTypeA + self.numTypeB) < mitosToAdd:
            self.addMito()

    # Mitochondrial replication code
    def addMito(self):
        randomNumber = np.random.random_sample()
        heteroplasmy = self.getCellHeteroplasmy()
        
        assert heteroplasmy >= 0 # check I've written heteroplasmy right!
        assert heteroplasmy <= 1

        if randomNumber < heteroplasmy:
            self.numTypeA += 1
        elif randomNumber >= heteroplasmy:
            self.numTypeB += 1
#        if randomNumber <= self.probTypeA:
#            self.numTypeA += 1
#        elif randomNumber > self.probTypeA:
#            self.numTypeB += 1
        else:
            print "Unexpected random value"
            raise ValueError

    # Cell division code
    def cellDivision(self):
        initialTypeA = self.numTypeA
        initialTypeB = self.numTypeB
       
        if self.numTypeA > 0:
            numToLoseA = np.random.binomial(self.numTypeA, 0.5)
        else:
            numToLoseA = 0

        if self.numTypeB > 0:
            numToLoseB = np.random.binomial(self.numTypeB, 0.5)
        else:
            numToLoseB = 0
        
        # take mitochondria out
        self.numTypeA -= numToLoseA
        self.numTypeB -= numToLoseB
        
        assert initialTypeA + initialTypeB == self.numTypeA + self.numTypeB + numToLoseA + numToLoseB

        return [numToLoseA, numToLoseB] # tell the petri dish what the new cell has
        
    def chooseAction(self):
        if self.numTypeA < self.targetMitos:
            return 'addMito'
            
        if self.numTypeA >= self.targetMitos:
            return 'divide'
        
    # Interfaces
    def getNumTypeA(self):
        return self.numTypeA

    def getNumTypeB(self):
        return self.numTypeB
    
    def getID(self):
        return self.ID
        
    def getCellHeteroplasmy(self):
        heteroplasmy = 0.0        
        typeA = float(self.numTypeA)
        typeB = float(self.numTypeB)
        
        if typeA == 0:
            return 0
        
        elif typeB == 0:
            return 1
            
        else:
            heteroplasmy = typeA/(typeA + typeB)
            return heteroplasmy

class petriDish(object):
    def __init__(self, numCells, probTypeA, cellMitosToAdd, initialA, initialB, targetMitos):
        # Store cells in this list
        self.cellList = []

        # Store initial variables
        self.numCells = numCells
        self.probTypeA = probTypeA
        self.cellMitosToAdd = cellMitosToAdd
        self.targetMitos = targetMitos
        self.initialA = initialA
        self.initialB = initialB

        # Add the initial cell population
        for i in range(0, numCells):
            self.newCell(self.probTypeA, self.cellMitosToAdd, self.initialA, self.initialB, self.targetMitos)

    def newCell(self, probTypeA, mitosToAdd, numTypeA, numTypeB, targetMitos):
        self.cellList.append(cell(probTypeA, mitosToAdd, numTypeA, numTypeB, targetMitos))

    def makeCellDivide(self, chosenCell):
        newCellInfo = chosenCell.cellDivision()
        self.newCell(self.probTypeA, 0, newCellInfo[0], newCellInfo[1], self.targetMitos)
        
    def makeCellAddMito(self, chosenCell):
        chosenCell.addMito()
        
    def chooseCell(self):
        return(np.random.choice(self.cellList))
        
    # Experiment loop:
    def loop(self, time):
        heteroplasmyTimeSeries = []        
        for t in range(0, time):
            chosen = self.chooseCell()
            if chosen.chooseAction() == 'addMito':
                chosen.addMito()
                
            elif chosen.chooseAction() == 'divide':
                self.makeCellDivide(chosen)
                
            else:
                print 'undetermined action'
                
            heteroplasmyTimeSeries.append(self.getHeteroplasmy())
            
        return heteroplasmyTimeSeries

    # Interfaces
    def getNumCells(self):
        return len(self.cellList)

    def printCellList(self): # really just for testing
        print 'ID, typeA, typeB, total'
        for x in self.cellList:
            print x.getID(), x.getNumTypeA(), x.getNumTypeB(), x.getNumTypeA() + x.getNumTypeB()
            
    def getHeteroplasmy(self):
        totalA = 0.0
        totalB = 0.0
        
        for x in self.cellList:
            totalA += x.getNumTypeA()
            totalB += x.getNumTypeB()
            
        heteroplasmy = totalA/(totalA + totalB)
        
        return heteroplasmy

def doBasicExperiment(numRuns, runTime, numCells, probA, mitosToAdd, initialA, initialB, targetMitos):
    for run in range(0, numRuns):
        petri = petriDish(numCells, probA, mitosToAdd, initialA, initialB, targetMitos)    
        result = petri.loop(runTime)
        mpl.pyplot.plot(result)
    
def test():
    petri = petriDish(10, 0.5, 0, 9, 1, 10)
    print petri.getNumCells()
    petri.printCellList()
    print 'testing division'
    chosen = petri.chooseCell()
    petri.makeCellDivide(chosen)
    petri.printCellList()
    print petri.getHeteroplasmy()
    result = petri.loop(10)
    for x in result:
        print x
        
    mpl.pyplot.plot(result)
    
doBasicExperiment(100, 1000, 10, 0.7, 0, 9, 1, 10)













	
