try: attach("LGGMultiScript.sage")
except: load("LGGMultiScript.sage")
#Prepare to save or load data.
from pickle import dump, load
import time
def loadFromFile(filename, Output=0):
    file = open(filename, 'rb')
    data = load(file)
    file.close()
    if Output == 1:
        print(data)
    return(data)

curvesDatabase = loadFromFile("type456PrimeCurvesUpTo10000")

#19:152 will give range (100,1100)
simplifiedDatabase = curvesDatabase[102:152]
radius = 1
batch = 0
nodeLimit = 1000
holdTime = time.time()
for i in range(len(simplifiedDatabase)):
    p = simplifiedDatabase[i][0]
    for j in range(len(simplifiedDatabase[i])-1):
        curveType = simplifiedDatabase[i][j+1][0]
        rootsInit = simplifiedDatabase[i][j+1][1]
        MultiScript(p,rootsInit,batch,radius,nodeLimit,curveType)
        batch += 0
        newTime = time.time()
        print(newTime-holdTime)
        holdTime = newTime









