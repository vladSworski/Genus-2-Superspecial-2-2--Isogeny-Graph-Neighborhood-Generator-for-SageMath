︠b0957185-0a93-4c6e-8d61-e9b611358a77︠
#Add your preferred directory.
directory = ""

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
︡16584b1a-eef3-4714-ac6d-b190f5490382︡{"done":true}
︠a8b1f478-3654-41d8-b0c3-cf7e6d4d2b58︠
def outputMe(LoadedData):
    #Seperate out data for the Spine.
    spinalData = splitSetToolSpinality(LoadedData)[0]
    #Generate dictionaries for the data.
    invarToLabelDict = invarToLabelDictCreator(LoadedData)
    invarToLabelDictSpine = invarToLabelDictCreator(spinalData)
    isogenyList,isogenyDict = isogenyDictCreator(LoadedData,invarToLabelDict)
    isogenyListSpinal,isogenyDictSpinal = isogenyDictCreator(spinalData,invarToLabelDictSpine)
    #Calculates the nodes of the graphs:
    nodes = [[n[0][0] for n in LoadedData],[n[0][1] for n in LoadedData]]
    spinalNodes = [[n[0][0] for n in spinalData],[n[0][1] for n in spinalData]]
    #Calculate appropriate matrices.
    stdMatrices = calculateAdjacencyMatrices(nodes[0],isogenyDict)
    spnMatrices = calculateAdjacencyMatrices(spinalNodes[0],isogenyDictSpinal)
    #Calculate the number of spinal components.
    connectedComps = spineComponents(spnMatrices[5])
    #Calculating four cycles data.
    #fourCycles = fourCyclesData(stdMatrices[3],nodes[0])
    #fourCyclesSpine = fourCyclesData(spnMatrices[3],spinalNodes[0])
    #fCCnbhd = fourCycleCounter(fourCycles[0],stdMatrices[2])
    #fCCspine = fourCycleCounter(fourCyclesSpine[0],spnMatrices[2])
    return(nodes,spinalNodes,isogenyList,isogenyListSpinal,connectedComps)#,fourCycles,fourCyclesSpine,fCCnbhd,fCCspine)

P = Primes()
piter = 0
p = P[piter]
while p < 599:
    t = 4
    while t < 7:
        skipLoop = False
        try:
            data = loadFromFile("nbhdData_d231002_b0/nbhdData_d231002_b0_r1_t"+str(t)+"_p"+str(p))
        except:
            try:
                data = loadFromFile("nbhdData_d231002_b0/nbhdData_d231002_b0_r1_t"+str(t)+"_p"+str(p))
            except:
                skipLoop = True
        if skipLoop == False:
            #Main body of calculations.
            global BF, YF, F, x, z, XF, z5
            BF = GF(p)
            YF.<y> = BF[]
            F.<z> = GF(p^2,modulus=y^2-simplestFieldModulus(p))
            XF.<x> = F[]
            
            outData = outputMe(data)
            saveToFile(outData,directory+"data_r2_t"+str(t)+"_p"+str(p))
            print(directory+"data_r2_t"+str(t)+"_p"+str(p))

            nodes = outData[0][0]
            spinalNodes = outData[1][0]
            isogenyList = outData[2]
            isogenyListSpinal = outData[3]

            plotNbhd(nodes,isogenyList,directory+"nGraph_r2_t"+str(t)+"_p"+str(p),1,0)
            plotNbhd(spinalNodes,isogenyListSpinal,directory+"sGraph_r2_t"+str(t)+"_p"+str(p),1,0)

        t += 1
    piter += 1
    p = P[piter]









