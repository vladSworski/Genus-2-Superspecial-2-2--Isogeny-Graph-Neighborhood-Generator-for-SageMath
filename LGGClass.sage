#LGGClass file for Local Graph Generator v1.5 (Fully commented.)
#Last updated 08/18/2023.
#Other files in our package: "LGGIsog.sage", "LGGUtilNew.sage", "LGGDataGen.sage", "VTools.sage"
#Our files for data production are "LGGScript.sagews", "LGGDataCollector.sagews"
#This file primarily constructs the nodes of our graph, represented by isomorphism classes of
#products of elliptic curves and Jacobians of hyperelliptic curves of genus 2.

#Setting up global arrays and dictionaries.
storedNodes = []
storedNodesByLabels = {}
storedNodesByInvars = {}
storedInvarsByLabels = {}
storedLabelsByInvars = {}

#Labels are determined by how many curves of a previous type are present.
globalTypeCounter = [0 for i in range(14)]

#Relation between curve types and their automorphism groups. 
typeRelate = [["H0","H1","H2","H3","H4","H5","H6","E23","E22","E12","E11","E02","E01","E00"],\
              ["C1","C2","S3","K4","D12","S4","C5orS5","C2","K4","C4","K4sdpC4","C6","C12","C6dpS3"]]

#Establishes an ordering to the labels for storage.
def labelDictOrder(label1,label2):
    #Finds the point of the label which separates the header from the footer.
    pos1 = label1.find("_")
    pos2 = label2.find("_")
    #Associates a number to each type so the header of the label has a value.
    value1 = typeRelate[0].index(label1[0:pos1])
    value2 = typeRelate[0].index(label2[0:pos2])
    #If header one is less than header two, label1 is less than label2.
    if value1 < value2:
        return(True)
    #If header two is less than header one, label2 is less than label1
    elif value2 < value1:
        return(False)
    #The headers are the same, look at the footers.
    else:
        newValue1 = str(label1[pos1:])
        newValue2 = str(label2[pos2:])
        #If footer one is the same or less than footer 2, label1 is less than label2
        if value1 <= value2:
            return(True)
        #Otherwise, label2 is less than label1
        else:
            return(False)

#A method that stores nodes in a way that we can keep track of.
def nodeStoreProcedure(neighbors):
    #For a given node, look at all of its neighbors.
    for i in range(len(neighbors)):
        #Have we already saved this node? No? Proceed.
        if (tuple(neighbors[i][1]) in storedNodesByInvars.keys()) == False:
            #Store the new node by its invariants.
            storedNodesByInvars[tuple(neighbors[i][1])] = neighbors[i][0]
            #Correctly assign the new node's label.
            labelFront = neighbors[i][0].curveType
            Index = multiIf(labelFront,["H0","H1","H2","H3","H4","H5","H6","E23","E22","E12","E11","E02","E01","E00"],[i for i in range(14)])
            globalTypeCounter[Index] += 1
            labelBack = globalTypeCounter[Index]
            label = labelFront+"_"+str(labelBack)
            #Assign the label to the node.
            neighbors[i][0].defineLabel(label)
            #Store the new node in an array.
            storedNodes.append(neighbors[i][0])
            #Store the new node by its label.
            storedNodesByLabels[label] = neighbors[i][0]
            #Relate the label to its invariants.
            storedLabelsByInvars[tuple(neighbors[i][1])] = label
            storedInvarsByLabels[label] = neighbors[i][1]

#This nodes represents an isomorphism class of the product of two supersingular elliptic curves.
class EllipticClass():
    #Builds an Elliptic Isomorphism Class to represent a node in our graph.
    def __init__(self, roots, layer=0):
        global storedNodes
        global storedNodesByLabels
        global storedNodesByInvars
        global storedInvarsByLabels
        global storedLabelsByInvars
        #How many steps from the origin of the graph are we?
        self.layer = layer
        #Marked true by default.
        self.markedForInclusion = True
        #Takes in roots, standardizes them to Rosenhain normal form (RNF) and also reverse Rosenhain normal form (RRNF).
        holdRoots = []
        adjustRoots = []
        for i in range(2):
            currentRoots = roots[i]
            if len(currentRoots) == 3:
                holdRoots.append(rosenhainLFT(currentRoots[0:2],[currentRoots[2]]))
            if len(currentRoots) == 4:
                holdRoots.append(rosenhainLFT(currentRoots[0:3],[currentRoots[3]]))
            adjustRoots.append(reverseRosenhainLFT(holdRoots[i],elliptic=True))
        self.rootsRNF = holdRoots
        self.rootsRRNF = adjustRoots
        #Elliptic curve methods tend to require Elliptic curve objects.  So we store those.
        coeffs = [self.function()[0].coefficients(sparse=False),self.function()[1].coefficients(sparse=False)]
        self.Elliptics = [EllipticCurve([0,coeffs[0][2],0,coeffs[0][1],coeffs[0][0]]),EllipticCurve([0,coeffs[1][2],0,coeffs[1][1],coeffs[1][0]])]
        #Calculate and store invariants.
        self.jInvariants = [self.Elliptics[0].j_invariant(),self.Elliptics[1].j_invariant()]
        #Invariants stored as a tuple.
        self.stdInvars = tuple(self.jInvariants)
        #Order the elliptic curves so we don't differentiate E1xE2 and E2xE1.
        self.orderElliptic()
        #Determine and save the type of the curve.
        self.curveType = self.typeDetector()
        #Is our curve in the Kohel spine?
        self.KohelSpinality = True
        for i in range(2):
            if self.jInvariants[i] not in BF:
                self.KohelSpinality = False
        #Is our curve superspecial?
        self.superspecial = True
        if self.Elliptics[0].is_supersingular() == False or self.Elliptics[1].is_supersingular() == False:
            self.superspecial = False
        #Store our node to the graph.
        nodeStoreProcedure([[self,self.stdInvars]])
        #Marked as False if the neighborhood hasn't been calculated yet.  Provides a good way to check without try/except.
        self.nCalc = False

    #Usually we store our curves as (R)RNF roots.  We can retreive the associated functions here.
    def function(self,degree=3):
        if degree == 3:
            return([rootsToPoly(self.rootsRNF[0]),rootsToPoly(self.rootsRNF[1])])
        if degree == 4:
            return([rootsToPoly(self.rootsRRNF[0]),rootsToPoly(self.rootsRRNF[1])])

    #Ensures E1xE2 and E2xE1 aren't stored as different curves.   Somewhat overkill to deal with a previous bug.
    def orderElliptic(self):
        if orderInFF(self.jInvariants[0],self.jInvariants[1],z) == False:
            self.rootsRNF = [self.rootsRNF[i-1] for i in range(len(self.rootsRNF))]
            self.rootsRRNF = [self.rootsRRNF[i-1] for i in range(len(self.rootsRRNF))]
            self.Elliptics = [self.Elliptics[i-1] for i in range(len(self.Elliptics))]
            self.jInvariants = [self.jInvariants[i-1] for i in range(len(self.jInvariants))]
            self.stdInvars = tuple([self.stdInvars[i-1] for i in range(len(self.stdInvars))])
        else:
            self.rootsRNF = [self.rootsRNF[i] for i in range(len(self.rootsRNF))]
            self.rootsRRNF = [self.rootsRRNF[i] for i in range(len(self.rootsRRNF))]
            self.Elliptics = [self.Elliptics[i] for i in range(len(self.Elliptics))]
            self.jInvariants = [self.jInvariants[i] for i in range(len(self.jInvariants))]
            self.stdInvars = tuple([self.stdInvars[i] for i in range(len(self.stdInvars))])

    #Determines the type of an elliptic product.  These types are entirely determined by j-invariants.
    #Also stores the automorphism group of the curve and its size.
    def typeDetector(self):
        #Rewrite j-invariants for simplicity.
        j1 = self.jInvariants[0]; j2 = self.jInvariants[1]
        if j1 == j2:
            if j1 == 0:
                #Type E00 curve.  Both have j-invariants 0.
                self.automorphismGroupSize = 36
                self.automorphismGroup = typeRelate[1][13]
                return("E00")
            elif j1 == 1728:
                #Type E11 curve.  Both have j-invariants 1.
                self.automorphismGroupSize = 16
                self.automorphismGroup = typeRelate[1][10]
                return("E11")
            else:
                #Type E22 curve.  Both are the same generic elliptic curve.
                self.automorphismGroupSize = 4
                self.automorphismGroup = typeRelate[1][8]
                return("E22")
        else:
            if j1 == 0:
                if j2 == 1728:
                    #Type E01 curve.  One has j-invariant 0, and the other j-invariant 1.
                    self.automorphismGroupSize = 12
                    self.automorphismGroup = typeRelate[1][12]
                    return("E01")
                else:
                    #Type E02 curve.  One has j-invariant 0, and the other is generic.
                    self.automorphismGroupSize = 6
                    self.automorphismGroup = typeRelate[1][11]
                    return("E02")
            elif j1 == 1728 or j2 == 1728:
                #Type E12 curve.  One has j-invariant 1728, and the other is generic.
                self.automorphismGroupSize = 4
                self.automorphismGroup = typeRelate[1][9]
                return("E12")
            else:
                #Type E23 curve.  Both are different generic elliptic curves.
                self.automorphismGroupSize = 2
                self.automorphismGroup = typeRelate[1][7]
                return("E23")

    #An exterior method to assign a label to the curve.
    def defineLabel(self,label):
        self.label = label

    
    #Generate the degree 2 isogenies of the node.
    def degree2Isogenies(self):
        #Box to store isogenies.
        Isogenies = []
        #Double Elliptic to Hyperelliptic
        for i in range(6):
            outCurve = isogenyEH(self.rootsRNF,self.function(),i)
            if outCurve[0] == True:
                Isogenies.append(HyperellipticClass(outCurve[1],layer = self.layer + 1))
        #Double Elliptic to Double Elliptic
        outCurve1 = isogenyEE(self.rootsRNF[0],self.Elliptics[0])
        outCurve2 = isogenyEE(self.rootsRNF[1],self.Elliptics[1])
        for i in range(3):
            for j in range(3):
                Isogenies.append(EllipticClass([outCurve1[i],outCurve2[j]],layer = self.layer + 1))
        self.degree2Isogenies = Isogenies
        
    #Turns the list of isogenies into detailed neighborhood data.
    def neighbors(self):
        #Stores invariants and curves.
        invarsBox = []
        curveBox = []
        #Takes each isogeny, plucks out the curve and invars, and then stores them.
        for i in range(len(self.degree2Isogenies)):
            curve = self.degree2Isogenies[i]
            invars = curve.stdInvars
            invarsBox.append(invars)
            curveBox.append(curve)
        #Pulls unique neighbors and counts their occurence.
        prepNeighbors = uniqueSetCounter(invarsBox)
        #Stores neighbor information to the curve.
        self.neighbors = [[self.degree2Isogenies[prepNeighbors[i][2][0]],prepNeighbors[i][0],prepNeighbors[i][1],\
                           [self.degree2Isogenies[prepNeighbors[i][2][j]].rootsRRNF for j in range(len(prepNeighbors[i][2]))]]\
                          for i in range(len(prepNeighbors))]
        #Adds these nodes to the global node lists.
        nodeStoreProcedure(self.neighbors)
        #Generate a dictionary version of the neighbors by invariant for ease of use.
        self.neighborsDict = {}
        for i in range(len(self.neighbors)):
            self.neighborsDict[tuple(self.neighbors[i][1])] = self.neighbors[i]
        #We've calculated the neighborhood now.
        self.nCalc = True
        
    #An exterior method to unmark the curve.
    def unMark(self):
        self.markedForInclusion = False

    #Takes all the pivotal information about this class object and reduces it to an easily pickle-able array.
    def collapseDataToString(self):
        #Note nodes use up around 72 bytes of data.
        #0 Label, 1 Invariants, 2 AutGroup, 3 roots for reconstruction, 4 spinality, 5 include me?,
        #6 type, 7 Layer
        DataStringUpper = [self.label,self.stdInvars,self.automorphismGroup,self.rootsRRNF, self.KohelSpinality,\
                           self.markedForInclusion, self.curveType, self.layer]
        #0 Invars of neighbors, 1 Weights to Neighbors, 2 Root Images of Neighbors
        DataStringLower = [[self.neighbors[i][1] for i in range(len(self.neighbors))],\
                           [self.neighbors[i][2] for i in range(len(self.neighbors))],\
                           [self.neighbors[i][3] for i in range(len(self.neighbors))]]
        self.compressedData = [DataStringUpper,DataStringLower]

#This node represents an isomorphism class of the Jacobian of a hyperelliptic curve.
#Note repeat methods will not be commented here.  Only information specific to this class.
#As of the writing of this comment, I have learned how to use inheritance, but it is not
#worth my time to rewrite this right now.
class HyperellipticClass():
    #Builds a Hyperelliptic Isomorphism class to represent a node in our graph.
    def __init__(self, roots, layer=0, manualInvars = False):
        global storedNodes
        global storedNodesByLabels
        global storedNodesByInvars
        global storedInvarsByLabels
        global storedLabelsByInvars
        self.layer = layer
        self.markedForInclusion = True
        if len(roots) == 5:
            self.rootsRNF = rosenhainLFT(roots[0:2],roots[2:])
        elif len(roots) == 6:
            self.rootsRNF = rosenhainLFT(roots[0:3],roots[3:])
        else:
            raise Exception("You must provide 5 or 6 roots.")
        self.rootsRRNF = reverseRosenhainLFT(self.rootsRNF)
        self.rootsRRNFSorted = self.rootsRRNF[3:].sort()
        #Checks if the hyperelliptic curve is superspecial using the Hasse-Witt matrix.
        self.superspecial = HasseWittCheck(self.function())
        #Calculates the automorphisms, group and group size for our curve.
        automorphismPrep = automorphismDetector(self.rootsRRNF)
        self.automorphismGroupSize = automorphismPrep[0]
        self.automorphisms = automorphismPrep[1]
        self.automorphismGroup = automorphismPrep[2][0]
        #Stores the curve type, which is determined by automorphisms.
        self.curveType = "H"+str(automorphismPrep[2][1])
        #I chose to manually construct the invariants.
        if manualInvars == False:
            self.Kohel = HyperellipticCurve(self.function()).absolute_igusa_invariants_kohel()
            self.Wamelen = HyperellipticCurve(self.function()).absolute_igusa_invariants_wamelen()
        else:
            self.Igusa = [symmetric2_6root(self.rootsRRNF),symmetric4_6root(self.rootsRRNF),\
                          symmetric6_6root(self.rootsRRNF),symmetric10_6root(self.rootsRRNF)]
            self.Kohel = IgusaToKohel(self.Igusa)
            self.Wamelen = IgusaToWamelen(self.Igusa)
        #We choose Kohel's invariants as our variants of choice.
        self.stdInvars = tuple(self.Kohel)
        self.KohelSpinality = True
        self.nCalc = False
        #Check that our curve has Kohel spinality.
        for i in range(3):
            if self.Kohel[i] not in BF:
                self.KohelSpinality = False
        #Isomorphism classes may have a "dependent representative," a single isomorphic curve that has a quadratic splitting with determinant 0,
        #despite the fact that all other curves in the class do not.  The author is uncertain of why this happens, but it is clear that this is not
        #representative of the class as a whole.  As such, we check for such a representative, and if our reverse rosenhain normal form coincides with
        #one, we find an alternative representative.  It is worth noting that this is done under the assumption that finding such an alternative
        #representative should be probabilistically very easy.  It is possible that some scenario may refute this conjecture.  If so, this part of the
        #code should be updated and a better effort to understand these "dependent representatives" made.
        self.dependentRepresentative = dependentRepresentativeCheck(self.rootsRRNF)
        if True in self.dependentRepresentative:
            dependencyFixer = dependentRepresentativeFix(roots,self.dependentRepresentative[0],self.dependentRepresentative[1])
            self.dependencyRoots = dependencyFixer[0]
            self.dependencyLFT = dependencyFixer[1]
        nodeStoreProcedure([[self,self.stdInvars]])

    def function(self,degree=5):
        if degree == 5:
            return(rootsToPoly(self.rootsRNF))
        if degree == 6:
            return(rootsToPoly(self.rootsRRNF))

    def defineLabel(self,label):
        self.label = label

    #Calculate the degree 2 isogenies of the curve.
    def degree2Isogenies(self):
        Isogenies = []
        #Roots are choosen based on the above dependency check.
        if True in self.dependentRepresentative:
            roots = self.dependencyRoots
        else:
            roots = self.rootsRRNF
        #Generate all possible quadratic pairings.
        quadraticPairings = quadraticPairingGenerator(roots)

        for i in range(15):
            quadraticPairing = quadraticPairings[i]
            #Construct the matrix D, which checks what type of image we will have.
            determinant = F(determinantDelta(quadraticPairing))
            if determinant != 0:
                #Our output will be the Jacobian of a genus 2 Hyperelliptic Curve.
                OutCurve = isogenyHH(quadraticPairing,determinant)
                Isogenies.append(HyperellipticClass(OutCurve,layer = self.layer + 1))
            else:
                #Our output will be a product of two Elliptic Curves.
                OutCurve = isogenyHE(quadraticPairing)
                Isogenies.append(EllipticClass(OutCurve,layer = self.layer + 1))
        self.degree2Isogenies = Isogenies
    
    def neighbors(self):
        invarsBox = []
        curveBox = []
        for i in range(len(self.degree2Isogenies)):
            curve = self.degree2Isogenies[i]
            invars = curve.stdInvars
            invarsBox.append(invars)
            curveBox.append(curve)
        prepNeighbors = uniqueSetCounter(invarsBox)
        self.neighbors = [[self.degree2Isogenies[prepNeighbors[i][2][0]],prepNeighbors[i][0],prepNeighbors[i][1],\
                           [self.degree2Isogenies[prepNeighbors[i][2][j]].rootsRRNF for j in range(len(prepNeighbors[i][2]))]]\
                          for i in range(len(prepNeighbors))]
        nodeStoreProcedure(self.neighbors)
        self.neighborsDict = {}
        for i in range(len(self.neighbors)):
            self.neighborsDict[tuple(self.neighbors[i][1])] = self.neighbors[i]
        self.nCalc = True
        
    def unMark(self):
        self.markedForInclusion = False
        
    def collapseDataToString(self):
        #0 Label, 1 Invariants, 2 AutGroup, 3 roots for reconstruction, 4 spinality, 5 include me?,
        #6 type, 7 Layer
        DataStringUpper = [self.label,self.stdInvars,self.automorphismGroup,self.rootsRRNF, self.KohelSpinality,\
                           self.markedForInclusion, self.curveType, self.layer]
        #0 Invars of neighbors, 1 Weights to Neighbors, 2 Root Images of Neighbors
        DataStringLower = [[self.neighbors[i][1] for i in range(len(self.neighbors))],\
                           [self.neighbors[i][2] for i in range(len(self.neighbors))],\
                           [self.neighbors[i][3] for i in range(len(self.neighbors))]]
        self.compressedData = [DataStringUpper,DataStringLower]