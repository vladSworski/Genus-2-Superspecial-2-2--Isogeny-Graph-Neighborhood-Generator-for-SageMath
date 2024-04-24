#LGGDataGen file for Local Graph Generator v1.5 (Fully commented.)
#Last updated 08/18/2023.
#Other files in our package: "LGGClass.sage", "LGGUtilNew.sage", "LGGIsog.sage", "VTools.sage"
#Our files for data production are "LGGScript.sagews", "LGGDataCollector.sagews"
#This file primarily provides methods used to calculate data for neighborhood graphs.

#Sometimes we don't care about the direction of an edge, but we store both directions separately.
#This forgets everything but that an edge exists.
def graphCompress(data):
    heldData = []
    newList = []
    for i in range(len(data)):
        #New edge?
        if (data[i][2]+"%"+data[i][1]) not in heldData:
            heldData.append(data[i][0])
            newList.append(data[i])
    return(newList)

#This tool splits a set of data into spinal nodes and non-spinal nodes.
def splitSetToolSpinality(data):
    IsSpinal = []
    NotSpinal = []
    for i in range(len(data)):
        #Check spinality.
        if data[i][0][4] == True:
            IsSpinal.append(data[i])
        else:
            NotSpinal.append(data[i])
    return(IsSpinal,NotSpinal)

#Very similar to above. Splits data into included and non-included nodes.
def splitSetToolInclusion(data):
    IsIncluded = []
    NotIncluded = []
    for i in range(len(data)):
        #Checks inclusion.
        if data[i][0][5] == True:
            IsIncluded.append(data[i])
        else:
            NotIncluded.append(data[i])
    return(IsIncluded,NotIncluded)

#Again similar.  Splits data by distance to origin node.
def splitSetToolLayer(data):
    #How many layers does the graph have?
    layerMax = max([n[0][7] for n in data])
    layerArray = [[] for i in range(layerMax+1)]
    #What layer is the node in?
    for i in range(len(data)):
        correctLayer = data[i][0][7]
        layerArray[correctLayer].append(data[i])
    return(layerArray)

#Creates a dictionary that takes in invariants and returns labels.
def invarToLabelDictCreator(data):
    invarToLabelDict = {}
    for i in range(len(data)):
        invarToLabelDict[tuple(data[i][0][1])] = data[i][0][0]
    return(invarToLabelDict)

#Creates a dictionary for edges by their labels.
def isogenyDictCreator(data,labelDict):
    #Produces both a list version and a dictionary version of the data.
    isogenyList = []
    isogenyDict = {}
    for i in range(len(data)):
        for j in range(len(data[i][1][0])):
            #We may want to exclude some edges.   This verifies we actually want to include this one.
            if tuple(data[i][1][0][j]) in labelDict.keys():
                #Build the label for the edge.
                firstLabel = data[i][0][0]
                secondLabel = labelDict[tuple(data[i][1][0][j])]
                NewHandle = firstLabel+"%"+secondLabel
                #0 Edge Label, 1. First node label. 2 Second node label. 3. Weight from first to second. (Second to first is its own edge.)
                newIsogeny = [NewHandle,firstLabel,secondLabel,data[i][1][1][j]]
                isogenyList.append(newIsogeny)
                isogenyDict[NewHandle] = newIsogeny
    return(isogenyList,isogenyDict)

#Calculates all variations of adjacency matrices for a given dataset.
def calculateAdjacencyMatrices(nodeSet,edgeSet):
    #Adjacency matrices are square.  Side length?
    dimension = len(nodeSet)
    import numpy as np
    #Prepare empty numpy arrays.
    adjBase,adjBaseNoWeight,diagonalMat = [np.zeros([dimension,dimension]) for i in range(3)]
    for i in range(dimension):
        for j in range(dimension):
            try:
                #Fill in edge weights.
                adjBase[i][j] = edgeSet[nodeSet[i]+"%"+nodeSet[j]][3]
                #If the above edge exists, fill in a 1. (Unweighted.)
                adjBaseNoWeight[i][j] = 1
            except:
                #Edge doesn't exist.
                adjBase[i][j] = 0
                adjBaseNoWeight[i][j] = 0
    #Prepare new matrices.
    adjNoLoop = adjBase
    adjNoLoopNoWeight = adjBaseNoWeight
    for i in range(dimension):
        #No loop matrices are the standard definition of adjacency matrices.  The diagonal is 0 here (no self-loops).
        adjNoLoop[i][i], adjNoLoopNoWeight[i][i] = [0,0]
    for i in range(dimension):
        #The diagonal matrix only has non-zero entries on the diagonal.  They are the sum of the weights of all outgoing edges to that node.
        diagonalMat[i][i] = sum(adjNoLoopNoWeight[i])
    #0 Basic adjacency matrix for the graph, 1 Ajacency matrix with all non-zero entires set to 1.  2 Basic adjacency matrix with all
    #diagonal entries set to 0.  3 The same as 1 and 2 at the same time.  4. The diagonal matrix of the graph. 5. Laplacian matrix of the graph.
    return(adjBase,adjBaseNoWeight,adjNoLoop,adjNoLoopNoWeight,diagonalMat,diagonalMat-adjNoLoopNoWeight)

#The rank of the null-space of the Laplacian is equal to the number of connected components of a graph.
#The standard neighborhood graph is always connected, but the Spine isn't necessarily.
def spineComponents(Laplace):
    dimension = len(Laplace)
    matrixLaplace = matrix(ZZ,dimension,dimension,Laplace)
    return(kernel(matrixLaplace).rank())

#Calculate all four-cycles in the neighborhood using the adjacency matrices.
#One of my older methods, could be more optimal.
def fourCyclesFromMatrix(M):
    import numpy as np
    #Find the dimension of the matrix.
    n = M.shape[0]
    #Matrix form of the array.
    A = matrix(M.shape[0],M.shape[1],M)
    #Construct a new matrix with n rows and one column for every potential pair of elements that could be part of a cycle one step from the origin.
    width = (n-1)*(n-2)/2
    Trial = np.zeros([n,width])
    Trial.shape
    #This loop prepares the above trial matrix with a pair of ones in every column to represent a pair of elements.
    uptick = 0
    shortTermDict = {}
    for i in range(n-2):
        I = i+1
        for j in range(n-I-1):
            J = j+I+1
            Trial[I][uptick] = 1
            Trial[J][uptick] = 1
            shortTermDict[uptick] = [I,J]
            uptick += 1
    #Matrix form of trial.
    P = matrix(n, width, Trial)
    #The product of the adjacency matrix with the pairs matrix.
    B = A*P
    #These pairs will be those that actually represent a 4-cycle.
    goodpairs = []
    for i in range(width):
        if B[0][i] == 2:
            for j in range(n-1):
                if B[j+1][i] == 2:
                    #If two positions in a given column are equal to two, this implies the existence of a cycle.
                    #(Two nodes have two common neighbors.)
                    goodpairs.append((j+1,i))
    #I don't think this does anything, but I'd rather not risk removing it.
    loopStart = 1
    #Store four-cycles.
    loopOutput = []
    for i in range(len(goodpairs)):
            loopOutput.append([0,shortTermDict[goodpairs[i][1]][0],goodpairs[i][0],shortTermDict[goodpairs[i][1]][1]])
            #Output will be a four-tuple with numerical indexing.
    return(loopOutput)

#Returns four cycles in a useful format.
def fourCyclesData(adjMat,nodes):
    #Runs the above method, returning four-tuples with numerical indexing.
    fourCycleListing = fourCyclesFromMatrix(adjMat)
    fourCycleData = []
    #This loop uses the same data and node list, to rewrite the four-tuples in terms of their labels.
    for i in range(len(fourCycleListing)):
        vect = []
        for j in range(4):
            vect.append(nodes[fourCycleListing[i][j]])
        fourCycleData.append(vect)
    return(fourCycleListing,fourCycleData)

#Returns the number of fourcycles in a four-cycle listing.
def fourCycleCounter(cycleList,weightMatrix):
    #The pure number of cycles.
    baseLength = len(cycleList)
    #Keeps track of all cycles choices.
    finalSum = 0
    for i in range(baseLength):
        loopProdForward = 1
        loopProdBackward = 1
        for j in range(4):
            #The forward product goes through the edges in the order 0,1,2,3,0.
            loopProdForward = loopProdForward*weightMatrix[cycleList[i][j%4]][cycleList[i][(j+1)%4]]
            #The backward product goes through the edges in the order 0,3,2,1,0.
            loopProdBackward = loopProdBackward*weightMatrix[cycleList[i][(j+1)%4]][cycleList[i][j%4]]
            #Either way, we multiply by the edge weight.
        finalSum += (loopProdForward + loopProdBackward)
    #0 The number of unique four cycles in the graph. 1. Returns the number of unique quadruples of isogenies that begin and end at the same base point.
    #i.e., 0. No weight, no direction. 1. Weight, direction.
    return(baseLength,finalSum)

#Plots a rudimentary graph of the neighborhood.  More as a proof of concept.  Cannot do weights or directions.  Can do loops.
def plotNbhd(nodes,isogenyList, title="default", layout = 0, show = 1):
    import matplotlib.pyplot as plt
    import networkx as nx
    #Build edges.
    newList = graphCompress(isogenyList)
    #Setup figure.
    fig = plt.figure(figsize=(12,12))
    H = nx.Graph()
    #Build nodes.
    H.add_nodes_from(nodes)
    edgeWeights = {}
    for i in range(len(newList)):
        H.add_edge(newList[i][1],newList[i][2])
    #Picks a layout for the graph.  Expirement!  Different ones may sometimes be better.
    if layout == 0:
        pos = nx.circular_layout(H)
    if layout == 1:
        pos = nx.kamada_kawai_layout(H)
    if layout == 2:
        pos = nx.spring_layout(H)
    nx.draw(H, pos, with_labels = True)
    #Show the graph in the sage worksheet?
    if show == 1:
        plt.show()
    #Saves the plot using the given title.
    plt.savefig(title+".png", format="PNG")