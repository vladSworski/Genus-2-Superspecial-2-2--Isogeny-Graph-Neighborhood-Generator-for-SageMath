︠769f3e5a-13fc-4d10-9040-ca3522626482s︠
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
︡c4b52bcb-11c8-448d-8c44-2e7b2e076beb︡{"stdout":"The following is a list of nodes in the main neighborhood graph: \n['H4_1', 'E22_1', 'E00_1', 'H1_1', 'H2_1', 'H1_2', 'H3_1', 'H3_2', 'E02_1', 'E02_2', 'E22_2', 'E23_1', 'E22_3', 'H1_3', 'H0_1', 'H1_4', 'H0_2', 'H1_5', 'H3_3', 'H1_6', 'H5_1', 'H2_2']\nThe following is a list of nodes in the spinal subgraph:\n['H4_1', 'E22_1', 'E00_1', 'H1_1', 'H2_1', 'H1_2', 'H1_5', 'H3_3', 'H1_6', 'H5_1', 'H2_2']\nThe spinal subgraph has 1 connected components."}︡{"stdout":"\nThe following is a complete list of 4-cycles in the neighborhood graph:"}︡{"stdout":"\n[['H4_1', 'E22_1', 'H1_2', 'H1_1'], ['H4_1', 'E22_1', 'E23_1', 'H1_1'], ['H4_1', 'E22_1', 'H1_2', 'H2_1'], ['H4_1', 'E22_1', 'E23_1', 'H2_1'], ['H4_1', 'H1_1', 'H1_2', 'H2_1'], ['H4_1', 'H1_1', 'E23_1', 'H2_1']]\nThe following is a complete list of 4-cycles in the spinal subgraph:\n[['H4_1', 'E22_1', 'H1_2', 'H1_1'], ['H4_1', 'E22_1', 'H1_2', 'H2_1'], ['H4_1', 'H1_1', 'H1_2', 'H2_1']]\nThe following is the number of undirected, unweighted 4-cycles and the number of directed, weighted 4-cycles in the neighborhood graph:\n(6, 72.0)\nThe following is the number of undirected, unweighted 4-cycles and the number of directed, weighted 4-cycles in the spinal subgraph:\n(3, 36.0)\n"}︡{"stdout":"A plot of the neighborhood graph:\n"}︡{"file":{"filename":"11f3f38a-bbb7-4be7-b46d-8b4883ab59d5.svg","show":true,"text":null,"uuid":"9d468131-a536-4eda-8457-1d1be29d8de4"},"once":false}︡{"stdout":"A plot of the spinal subgraph:\n"}︡{"file":{"filename":"562ad4a6-ad57-4398-8081-4cc2677e337e.svg","show":true,"text":null,"uuid":"a7253bed-5243-4efb-98c8-de8b85658b20"},"once":false}︡{"done":true}
︠feb6edb0-31fb-4b8a-b235-6c7e5fc971d8︠









