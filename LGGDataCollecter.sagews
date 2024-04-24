#LGGDataCollector file for Local Graph Generator v1.5 (Fully commented.)
#Last updated 08/18/2023.
#Files in our package: "LGGDataGen.sage", LGGClass.sage", "LGGUtilNew.sage", "LGGIsog.sage", "VTools.sage"
#Our files for data production are "LGGScript.sagews", "LGGDataCollector.sagews"
#This file takes neighborhood data from a saved pickle file and returns graph data.

#Load necessary assets.
try: attach("VTools.sage","LGGIsog.sage","LGGClass.sage","LGGUtilNew.sage","LGGDataGen.sage")
except: load("VTools.sage","LGGIsog.sage","LGGClass.sage","LGGUtilNew.sage","LGGDataGen.sage")
import time
import datetime
from pickle import dump, load
import sys
import numpy as np

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

#Initialize fields. (Edit me as needed.)
p = 101 #Field characteristic.   Default to 1021.
approxCount = floor(p^3/2880)
F.<z> = GF(p^2,modulus=x^2-simplestFieldModulus(p))
BF = GF(p)
XF.<x> = F[]

#Load Data file.
string = "nbhdData_d240424_b0"
LoadedData = loadFromFile(string)

#Method that runs submethods for data collection.
def outputMe(LoadedData):
    #Seperate out data for the Spine.
    spinalData = splitSetToolSpinality(LoadedData)[0]
    #Generate dictionaries for the data.
    invarToLabelDict = invarToLabelDictCreator(LoadedData)
    invarToLabelDictSpine = invarToLabelDictCreator(spinalData)
    isogenyList,isogenyDict = isogenyDictCreator(LoadedData,invarToLabelDict)
    isogenyListSpinal,isogenyDictSpinal = isogenyDictCreator(spinalData,invarToLabelDictSpine)
    #Calculates the nodes of the graphs:
    nodes = [n[0][0] for n in LoadedData]
    print("The following is a list of nodes in the main neighborhood graph: ")
    print(nodes)
    spinalNodes = [n[0][0] for n in spinalData]; spinalNodes
    print("The following is a list of nodes in the spinal subgraph:")
    print(spinalNodes)
    #Calculate appropriate matrices.
    stdMatrices = calculateAdjacencyMatrices(nodes,isogenyDict)
    spnMatrices = calculateAdjacencyMatrices(spinalNodes,isogenyDictSpinal)
    #Calculate the number of spinal components.
    print("The spinal subgraph has "+str(spineComponents(spnMatrices[5]))+" connected components.")
    #Calculating four cycles data.
    fourCycles = fourCyclesData(stdMatrices[3],nodes);
    print("The following is a complete list of 4-cycles in the neighborhood graph:")
    print(fourCycles[1])
    fourCyclesSpine = fourCyclesData(spnMatrices[3],spinalNodes)
    print("The following is a complete list of 4-cycles in the spinal subgraph:")
    print(fourCyclesSpine[1])
    print("The following is the number of undirected, unweighted 4-cycles and the number of directed, weighted 4-cycles in the neighborhood graph:")
    print(fourCycleCounter(fourCycles[0],stdMatrices[2]))
    print("The following is the number of undirected, unweighted 4-cycles and the number of directed, weighted 4-cycles in the spinal subgraph:")
    print(fourCycleCounter(fourCyclesSpine[0],spnMatrices[2]))
    return(nodes,spinalNodes,isogenyList,isogenyListSpinal)

#Run the main method.
nodes,spinalNodes,isogenyList,isogenyListSpinal = outputMe(LoadedData)

#Plot the network graphs.
print("A plot of the neighborhood graph:")
plotNbhd(nodes,isogenyList,"graph"+str(p),1,1)
print("A plot of the spinal subgraph:")
plotNbhd(spinalNodes,isogenyListSpinal,"spinalGraph"+str(p),2,1)









