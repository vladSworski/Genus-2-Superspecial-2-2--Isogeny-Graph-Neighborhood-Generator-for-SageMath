#LGGUtilNew file for Local Graph Generator v1.5 (Fully commented.)
#Last updated 08/18/2023.
#Other files in our package: "LGGClass.sage", "LGGDataGen.sage", "LGGIsog.sage", "VTools.sage"
#Our files for data production are "LGGScript.sagews", "LGGDataCollector.sagews"
#This file contains broadly general methods for this project.  There also exists the file "LGGUtilOld.sage"
#which contains methods from past versions.

#Checks a hyperelliptic curve for superspeciality.
def HasseWittCheck(f):
    g = f^(int((p-1)/2))
    h = g.coefficients(sparse=False)
    #The following is the condition for a curve to be Superspecial.
    if h[p-1] != 0 or h[p-2] != 0 or h[2*p-1] != 0 or h[2*p-2] != 0:
        return(False)
    else:
        return(True)

#This method is an LFT that returns curves with roots infinity,0,1, the requirement to be in "Rosenhain Normal Form"
#(RNF). The elliptic toggle determines whether we are doing this for a Hyperelliptic Curve or an Elliptic Curve.
def rosenhainLFT(InRoot,OutRoot):
    if len(InRoot) == 2:
        def LFT(rootIn):
            Out = (F(rootIn) - F(InRoot[0]))/(F(InRoot[1])-F(InRoot[0]))
            return(Out)
        finalOut = [F(0),F(1)]
    elif len(InRoot) == 3:
        def LFT(rootIn):
            S = (InRoot[2]-InRoot[0])/(InRoot[2]-InRoot[1])
            Out = F(S)*(F(rootIn)-F(InRoot[1]))/(F(rootIn) - F(InRoot[0]))
            return(Out)
        finalOut = [F(0),F(1)]
    else:
        raise exception("Please enter a root list of length 2 or length 3.")
    
    for i in range(len(OutRoot)):
        finalOut.append(LFT(OutRoot[i]))
    return(finalOut)

#This method is an LFT that returns curves with roots 0,1,-1, the requirement to be in "Reverse Rosenhain Normal Form"
#(RRNF). The elliptic toggle determines whether we are doing this for a Hyperelliptic Curve or an Elliptic Curve.
def reverseRosenhainLFT(roots,elliptic=False):
    if elliptic == False:
        s1 = F(roots[2]);s2 = F(roots[3]);s3 = F(roots[4])
        if F(1/2) not in roots:
            outRoot = [F(1/(1-2*s1)),F(1/(1-2*s2)),F(1/(1-2*s3))]
        elif F(-1) not in roots:
            outRoot = [F((s1-1)/(s1+1)),F((s2-1)/(s2+1)),F((s3-1)/(s3+1))]
        elif F(2) not in roots:
            outRoot = [F(s1/(2-s1)),F(s2/(2-s2)),F(s3/(2-s3))]
        else:
            outRoot = [F(1/3),F(1/5),F(1/2)]
        out = [F(0),F(1),F(-1),outRoot[0],outRoot[1],outRoot[2]]
        return(out)
    else:
        s = F(roots[-1])
        if F(1/2) not in roots:
            return([F(0),F(1),F(-1),F(1/(1-2*s))])
        else:
            return([F(0),F(1),F(-1),F(-1/3)])

#A more specific version of the above LFT.
def specifiedReverseRosenhain6LFT(s1,s2,s3,s4,s5,s6,auto=0):
    C1 = F(s2-s3)
    C2 = F(s3+s2-2*s1)
    C3 = F(s1*s2+s1*s3-2*s2*s3)
    B4 = F((s4-s1)*C1/(C2*s4+C3))
    B5 = F((s5-s1)*C1/(C2*s5+C3))
    B6 = F((s6-s1)*C1/(C2*s6+C3))
    if auto == 0:
        return(F(0),F(1),F(-1),B4,B5,B6)
    if auto > 0:
        morph = [F(1),-s1,C2/C1,C3/C1]
    if auto == 1:
        return(morph)
    if auto == 2:
        return(F(0),F(1),F(-1),B4,B5,B6,morph)

#Determines when a representative is Dependent. This occurs when the determinant of the Wronskian of one of the representative's
#quadratic splittings is 0, but no other representative is.  There are two ways this can occur.  If the sum of every pair is the same,
#or when the product of every pair is the same.  The former is harder to fix and sadly more common.
def dependentRepresentativeCheck(roots):
    checkSum = []
    checkProd = []
    index = pairMe(6)
    for i in range(15):
        if roots[index[i][0]]+roots[index[i][1]] == roots[index[i][2]]+roots[index[i][3]] and\
        roots[index[i][0]]+roots[index[i][1]] == roots[index[i][4]]+roots[index[i][5]]:
            checkSum.append(1)
        else:
            checkSum.append(0)
        if roots[index[i][0]]*roots[index[i][1]] == roots[index[i][2]]*roots[index[i][3]] and\
        roots[index[i][0]]*roots[index[i][1]] == roots[index[i][4]]*roots[index[i][5]]:
            checkProd.append(1)
        else:
            checkProd.append(0)
    toggle1 = False
    toggle2 = False
    if 1 in checkSum:
        toggle1 = True
    if 1 in checkProd:
        toggle2 = True
    #We return (False, False) if all is well, (True, False) if we have a sum error, (False, True) for a product error.
    #(True, True) cannot happen unless the roots aren't unique.  This will never happen in our program.
    return(toggle1,toggle2)

#Fixes the above dependent representative problem.
def dependentRepresentativeFix(roots,toggle1,toggle2):
    if toggle1 == False and toggle2 == False:
        return(roots)
    #Fixes a product type problem.
    if toggle2 == True:
        #Add one to every root and the product problem immediately vanishes. i.e. the LFT z+1. ([1,1,0,1]).
        return([n+1 for n in roots])
    
    #Fixes a sum type problem.
    if toggle1 == True:
        toggle3 = True
        n = 0
        m = 1
        #We could get stuck forever if probability isn't in our favor.  Unlikely, but we give the code 100 chances before throwing an error.
        abortCycleKey = 0
        while toggle3 == True:
            #If m = n we don't get an LFT.  If m is a root, then we leave the degree 6 model.
            if n != m and m not in roots:
                #Builds an LFT to try to fix the problem. ([1,-n,1,-m])
                def quickLFT(roots,n,m):
                    return([(roots[i]-n)/(roots[i]-m) for i in range(6)])
                testRoot = quickLFT(roots,n,m)
                #Check if this fixes the problem.
                check = dependentRepresentativeCheck(testRoot)
                if True not in check:
                    #Problem fixed.
                    break
            #Arbitrary but helps vary stuff a bit.
            if m < n:
                m+=1
            else:
                n+=1
                m=0
            abortCycleKey += 1
            if abortCycleKey == 100:
                print("Attempted 100 different alternate represenatives with no luck.  Aborting.")
                toggle3 = False
        #We return the fixed roots, and what fixed them.
        return(testRoot,[1,-n,1,-m])
    
#This procedure determines what LFTs defining symmetries of some of the roots result in full automorphisms.
def automorphismDetector(roots):
    box = []
    #6 Choices for what to pair with 0.
    for i in range(6):
        holdRootsI = roots.copy()
        A = holdRootsI.pop(i)
        #5 Choices for what to pair with 1.
        for j in range(5):
            holdRootsJ = holdRootsI.copy()
            B = holdRootsJ.pop(j)
            #4 Choices for what to pair with -1.
            for k in range(4):
                holdRootsK = holdRootsJ.copy()
                C = holdRootsK.pop(k)
                #Sorts into roots used to generate the LFT and those the LFT acts on.
                index = [A,B,C]
                dualindex = holdRootsK.copy()
                #Uses a specified RRLFT to generate the map.
                try:
                    box.append([specifiedReverseRosenhain6LFT(index[0],index[1],index[2],dualindex[0],dualindex[1],dualindex[2],auto=2),i,j,k,index,dualindex])
                except:
                    box.append(["inf",i,j,k])
    newbox = []
    for i in range(len(box)):
        #If one of the roots goes to infinity, it immediately cannot be an automorphism.
        if box[i][0] != "inf":
            #Check the old set of roots matches the new set.
            if all([a in box[i][0] for a in roots]):
                newbox.append(box[i])
    indexDict = {}
    #Pairs roots with a number.
    for i in range(6):
        indexDict[F(roots[i])] = i
    #We build up the autmorphisms as elements of the symmetric group.
    automorphismsAsGroup = []
    for i in range(len(newbox)):
        string = ""
        for j in range(6):
            A1 = newbox[i][0][6][0]
            B1 = newbox[i][0][6][1]
            C1 = newbox[i][0][6][2]
            D1 = newbox[i][0][6][3]
            string = string + str(j) + str(indexDict[(F(A1)*F(roots[j])+F(B1))/(F(C1)*F(roots[j])+F(D1))])
        #Even pairs of characters in the string tell you where an element is sent.  Always of form 0a1b2c3d4e5f.
        automorphismsAsGroup.append(string)
    
    #Here we check the orders of each element.
    autoOrders = []
    for i in range(len(automorphismsAsGroup)):
        auto = automorphismsAsGroup[i]
        comp = auto
        order = 1
        while comp != "001122334455":
            comp = symmetrySum(comp,auto,6)
            order += 1
        autoOrders.append([auto,order])
    
    orderBox = [autoOrders[i][1] for i in range(len(autoOrders))]
    compileMe = [[newbox[i][0][6],autoOrders[i][0],autoOrders[i][1]] for i in range(len(autoOrders))]
    
    #len(compileMe) returns the size of the automorphism group, compileMe returns a listing of all automorphisms and their orders, groupDetector(OrderBox)
    #returns the name of the group
    return(len(compileMe),compileMe,groupDetector(orderBox))

#Mocks the group operator of the symmetric group on the strings in the abvoe method.
def symmetrySum(string1,string2,length):
    def fString1(x):
        L = len(string1)/2
        for i in range(L):
            if x == int(string1[2*i]):
                return(string1[2*i+1])
        return(str(x))
    
    def fString2(x):
        L = len(string2)/2
        for i in range(L):
            if x == int(string2[2*i]):
                return(string2[2*i+1])
        return(int(x))
    
    outstring = ""
    for i in range(length):
        outstring = outstring + str(i) + str(fString2(int(fString1(i))))
    return(outstring)

#Returns the type and name of the automorphism group based on its orders.
def groupDetector(orders):
    orders.sort()
    if orders == [1]:
        return("C1","0")
    elif orders == [1,2]:
        return("C2,","1")
    elif orders == [1,2,2,2,3,3]:
        return("S3","2")
    elif orders == [1,2,2,2]:
        return("K4","3")
    elif orders == [1,2,2,2,2,2,2,2,3,3,6,6]:
        return("D12","4")
    elif orders == [1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4]:
        return("S4","5")
    elif orders == [1,5,5,5,5]:
        return("C5","6")
    elif len(orders) == 120:
        return("S5","6")
    else:
        #This element shouldn't occur in most cases.
        return("?"+str(len(orders)),"?")

#Finds the unique elements in a set.
def uniqueSetElements(mySet, stringSwitch1 = False):
    outSet = []
    if stringSwitch1 == False:
        for i in range(len(mySet)):
            if mySet[i] not in outSet:
                outSet.append(mySet[i])
    else:
        for i in range(len(mySet)):
            if str(mySet[i]) not in outSet:
                outSet.append(str(mySet[i]))
    return(outSet)

#Counts how many of each unique element in a set there are.
def uniqueSetCounter(mySet,stringSwitch = False):
    uniques = uniqueSetElements(mySet, stringSwitch1 = stringSwitch)
    preppedBox = [[n,0,[]] for n in uniques]
    if stringSwitch == False:
        for i in range(len(uniques)):
            for j in range(len(mySet)):
                if uniques[i] == mySet[j]:
                    preppedBox[i][1] += 1
                    preppedBox[i][2].append(j)
    else:
        for i in range(len(uniques)):
            for j in range(len(mySet)):
                if str(uniques[i]) == str(mySet[j]):
                    preppedBox[i][1] += 1
                    preppedBox[i][2].append(j)
    return(preppedBox)

#Applies Igusa's symmetric polynomials method for calculating the first invariant.
def symmetric2_6root(data,scalar=1):
    sumMe = 0
    parity = pairMe(6)
    for i in range(len(parity)):
        sumMe += scalar*(data[parity[i][0]]-data[parity[i][1]])^2*(data[parity[i][2]]-data[parity[i][3]])^2*(data[parity[i][4]]-data[parity[i][5]])^2
    return(sumMe)

#Applies Igusa's symmetric polynomials method for calculating the second invariant.
def symmetric4_6root(data,scalar=1):
    sumMe = 0
    for i in range(20):
        outindex = [0,1,2,3,4,5]
        index = orderedChoose(6,3,i)
        outindex.remove(index[0])
        outindex.remove(index[1])
        outindex.remove(index[2])
        firstThree = (data[index[0]]-data[index[1]])^2*(data[index[1]]-data[index[2]])^2*(data[index[0]]-data[index[2]])^2
        lastThree = (data[outindex[0]]-data[outindex[1]])^2*(data[outindex[1]]-data[outindex[2]])^2*(data[outindex[0]]-data[outindex[2]])^2
        sumMe += scalar*firstThree*lastThree
    return(sumMe/2)

#Applies Igusa's symmetric polynomials method for calculating the third invariant.
def symmetric6_6root(data,scalar=1):
    sumMe = 0
    for i in range(20):
        outindex = [0,1,2,3,4,5]
        index = orderedChoose(6,3,i)
        outindex.remove(index[0])
        outindex.remove(index[1])
        outindex.remove(index[2])
        firstThree = (data[index[0]]-data[index[1]])^2*(data[index[1]]-data[index[2]])^2*(data[index[0]]-data[index[2]])^2
        midThree = (data[outindex[0]]-data[outindex[1]])^2*(data[outindex[1]]-data[outindex[2]])^2*(data[outindex[0]]-data[outindex[2]])^2
        firstSix = scalar*firstThree*midThree
        sumMe2 = 0
        for j in range(3):
            for k in range(2):
                lowerList = outindex.copy()
                A = lowerList.pop(j)
                B = lowerList.pop(k)
                C = lowerList[0]
                newindex = [A,B,C]
                sumMe2 += (data[index[0]]-data[newindex[0]])^2*(data[index[1]]-data[newindex[1]])^2*(data[index[2]]-data[newindex[2]])^2
        sumMe += firstSix*sumMe2
    return(sumMe/2)

#Applies Igusa's symmetric polynomials method for calculating the fourth invariant, aka the discriminant.
def symmetric10_6root(data,scalar=1):
    prodMe = scalar
    parity = pairMe(6)
    for i in range(15):
        index = orderedChoose(6,2,i)
        prodMe = prodMe*(data[index[0]]-data[index[1]])^2
    return(prodMe)

#Translates Igusa's invariants to Wamelen's.
def IgusaToWamelen(values):
    a = values[0]
    b = values[1]
    c = values[2]
    d = values[3]
    return(a^5/d,a^3*b/d,a^2*c/d)

#Translates Igusa's invariants to Kohel's
def IgusaToKohel(values):
    a = values[0]
    b = values[1]
    c = values[2]
    d = values[3]
    return(b*c/d,a^3*b/d,a^2*c/d)

#Scales Igusa's invariants by a factor k.
def scaleIgusa(values,k):
    a = values[0]
    b = values[1]
    c = values[2]
    d = values[3]
    return(k*a,k^2*b,k^3*c,k^5*d)