def MultiScript(p,rootsInit,batch,radius,nodeLimit,curveType):
    #Load necessary assets
    try: attach("VTools.sage","LGGIsog.sage","LGGClass.sage","LGGUtilNew.sage","LGGDataGen.sage")
    except: load("VTools.sage","LGGIsog.sage","LGGClass.sage","LGGUtilNew.sage",'LGGDataGen.sage')
    try: attach("LGGMultiScript.sage")
    except: load("LGGMultiScript.sage")
    import time
    import datetime
    from pickle import dump, load
    import sys
    datetime = datetime.datetime.now()
    date = datetime.strftime("%y")+datetime.strftime("%m")+datetime.strftime("%d")
    programStart = time.time()

    #Setting up global arrays and dictionaries.

    #Field Setup
    global BF, YF, F, x, z, XF, z5
    BF = GF(p)
    YF.<y> = BF[]
    F.<z> = GF(p^2,modulus=y^2-simplestFieldModulus(p))
    XF.<x> = F[]
    try: z5 = find5thRU(p)
    except: pass
    
    try: attach("VTools.sage","LGGIsog.sage","LGGClass.sage","LGGUtilNew.sage","LGGDataGen.sage")
    except: load("VTools.sage","LGGIsog.sage","LGGClass.sage","LGGUtilNew.sage",'LGGDataGen.sage')
    try: attach("LGGMultiScript.sage")
    except: load("LGGMultiScript.sage")
    roots = reverseRosenhainLFT(rosenhainLFT([F(rootsInit[0]),F(rootsInit[1])],[F(rootsInit[2]),F(rootsInit[3]),F(rootsInit[4])]))
    
    #Prepare to save or load data.
    def saveToFile(SaveMe,filename):
        file = open(filename, 'wb')
        dump(SaveMe, file)
        file.close()

    #The main algorithm we run on the data.
    def mainAlgorithm(R0,radius,nodeLimit):
        #Initialize program on first node.
        print(F)
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
        saveToFile(allData,"nbhdData_d"+str(date)+"_b"+str(batch)+"_r"+str(radius)+"_t"+str(curveType)+"_p"+str(p))
        return(graphNodes)
    
    mainAlgorithm(roots,radius,nodeLimit)