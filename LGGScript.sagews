#LGGScript file for Local Graph Generator v1.5 (Fully commented.)
#Last updated 08/18/2023.
#Files in our package: "LGGDataGen.sage", LGGClass.sage", "LGGUtilNew.sage", "LGGIsog.sage", "VTools.sage"
#Our files for data production are "LGGScript.sagews", "LGGDataCollector.sagews"
#This file generates a single neighborhood for a single prime from a starting curve.

#Load necessary assets
try: attach("VTools.sage","LGGIsog.sage","LGGClass.sage","LGGUtilNew.sage","LGGDataGen.sage")
except: load("VTools.sage","LGGIsog.sage","LGGClass.sage","LGGUtilNew.sage",'LGGDataGen.sage')
import time
import datetime
from pickle import dump, load
import sys
datetime = datetime.datetime.now()
date = datetime.strftime("%y")+datetime.strftime("%m")+datetime.strftime("%d")
programStart = time.time()

#Change batch number if you calculate more than one thing in a day.
batch = 0

#Initialize fields. (Edit me as needed.)
p = 101  #Field characteristic.   Default to 1021.
approxCount = floor(p^3/2880)
F.<z> = GF(p^2,modulus=x^2-simplestFieldModulus(p))
BF = GF(p)
XF.<x> = F[]

#Define inputs.
R0 = polyToRoots(x^6-1)
radius = 2 #Number of steps from central node.   Default is 2
nodeLimit = 100000 #A limit to the number of nodes the graph can count, just in case.  Default is 100.
    
#Prepare to save or load data.
def saveToFile(SaveMe,filename):
    file = open(filename, 'wb')
    dump(SaveMe, file)
    file.close()

def loadFromFile(filename, Output=0):
    file = open(filename, 'rb')
    data = load(file)
    file.close()
    if Output == 1:
        print(data)
    return(data)

#The main algorithm we run on the data.
def mainAlgorithm(R0,radius,nodeLimit):
    algorithmStart = time.time()
    timeCheck = algorithmStart
    #Initialize program on first node.
    if len(R0) == 6:
        C0 = HyperellipticClass(R0)
    elif len(R0) == 2:
        C0 = EllipticClass(R0)
    C0.degree2Isogenies()
    C0.neighbors()
    #Conditions to keep the while loop operating smoothly.
    LoopRunner = True
    nodeCount = 1
    while LoopRunner == True:
        #Run procedure on a new node.
        newNode = storedNodes[nodeCount]
        newNode.degree2Isogenies()
        newNode.neighbors()
        nodeCount += 1

        #Provides an update on progress of loop.
        if nodeCount%50 == 0:
            newTime = time.time()
            print(str(newTime-timeCheck)+" seconds have passed since the last time check. "+str(newTime-algorithmStart)+" seconds have passed in total. "+str(nodeCount)+" nodes have been calculated. This is approximately, "+str(numerical_approx(nodeCount*100/approxCount,digits=5))+"% of the nodes in the graph.")
            timeCheck = newTime

        #Loop exit instructions.
        if nodeCount >= nodeLimit:
            LoopRunner = False #Here we have calculated the maximum number of nodes we specified.
        if nodeCount >= len(storedNodes):
            LoopRunner = False #We calculated the entire graph.  There's nothing left to calculate.
        elif storedNodes[nodeCount].layer > radius:
            LoopRunner = False #We've exceeded the maximum distance from the starting node.
    
    #We now limit our nodal list to just those which will be used.
    saveTimeStart = time.time()
    graphNodes = storedNodes[0:nodeCount]
    #Compress data for saving.
    for i in range(len(graphNodes)):
        graphNodes[i].collapseDataToString()
    allData = [graphNodes[i].compressedData for i in range(len(graphNodes))]
    saveToFile(allData,"nbhdData_d"+str(date)+"_b"+str(batch))
    saveTimeEnd = time.time()
    print("Data save took "+str(saveTimeEnd - saveTimeStart)+" seconds. Data saved successfully to file")
    print("Total time: "+str(saveTimeEnd-algorithmStart)+" seconds.")
    return(graphNodes)

#Run the main method.
graphNodes = mainAlgorithm(R0,radius,nodeLimit)









