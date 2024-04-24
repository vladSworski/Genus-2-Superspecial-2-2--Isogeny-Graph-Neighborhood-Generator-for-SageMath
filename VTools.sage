#VTools file for Local Graph Generator v1.5 (Fully commented.)
#Last updated 08/18/2023.
#Other files in our package: "LGGClass.sage", "LGGUtilNew.sage", "LGGIsog.sage", "LGGDataGen.sage"
#Our files for data production are "LGGScript.sagews", "LGGDataCollector.sagews"
#This file contains methods that the author deems useful outside of this project.  Most are utility methods.
#The author is fond of collapsing loops to be over a single index.

import numpy as np
from pickle import dump, load

#Returns the general isogeny of an elliptic curve for a 2-torsion point (the built in isogeny method doesn't work for general fields)
def veluIsogenyDegree2(equation,xq):
    coeffs = equation.coefficients(sparse=False)
    b2 = coeffs[2]; b4 = coeffs[1]; b6 = coeffs[0]
    gqx = 3*xq^2 + 2*b2*xq + b4
    v = gqx; w = gqx*xq
    return(x^3 + (b4-5*v)*x + (b6-7*w))

#Provides a definition for the order of elements in a finite field of size p^2.
#Def: (a+b*z) < (c+d*z) if b < d or b = d and a < c. (Dictionary order favoring z).
def orderInFF(A,B,z):
    A = F(A)
    B = F(B)
    #This is the best way I could find to separate a and b from (a+b*z).
    Ai = (A^p-A)/(z^p-z)
    Ar = A-z*Ai
    Bi = (B^p-B)/(z^p-z)
    Br = B-z*Bi
    #Conditions as outlined above.
    if Ai < Bi:
        return(True)
    if Ai > Bi:
        return(False)
    if Ai == Bi:
        if Ar < Br:
            return(True)
        if Ar > Br:
            return(False)
        if Ar == Br:
            return(True)

#This method is here for the sake of completeness, but does not run well in an attached/loaded file.
#Please copy/paste this into your main worksheet.
#Saves to pickle file.
def saveToFile(SaveMe,filename):
    file = open(filename, 'wb')
    dump(SaveMe, file)
    file.close()

#This method is here for the sake of completeness, but does not run well in an attached/loaded file.
#Please copy/paste this into your main worksheet.
#Loads from pickle file.
def loadFromFile(filename, Output=0):
    file = open(filename, 'rb')
    data = load(file)
    file.close()
    #See data immediately?  No by default.
    if Output == 1:
        print(data)
    return(data)

#The author believes that the best element to represent the new basis element in a field of size p^2, is i if i is not in the base field.
#Otherwise, the author prefers to use the smallest prime that isn't square in the base field.  This method finds out what the best
#'z' is under those constraints.  Note that this doesn't work well for p=2,3.  But then again, neither does anything else the author has written.
def simplestFieldModulus(char):
    if kronecker(-1,char) == -1:
        return(-1)
    else:
        P = Primes()
        I = 0
        while P[I] < (char-1):
            if kronecker(P[I],char) == -1:
                return(P[I])
            I += 1
        else:
            return("This shouldn't happen.")

#Provides an easy method to iterate over every choice of elements in a set.
#Example use: orderedChoose(5,17,31).  This returns the 31st way to choose 5 elements from a set of 17.  Can run range from 0 to (17 choose 5).
def orderedChoose(objects, choose, iterative):
    triple = [objects, choose, iterative+1]
    finalOut = []
    for i in range(choose):
        chooseList = [binomial(triple[1]+j-1,triple[1]) for j in range(triple[0]-triple[1]+2)]
        adjustList = [chooseList[k] - triple[2] for k in range(len(chooseList))]
        l = 0
        while adjustList[l] < 0:
            l += 1
        index = l+triple[1]-1
        finalOut.append(index)
        newIter = triple[2] - chooseList[l-1]
        triple[1] -= 1
        triple[2] = newIter
    Out = [finalOut[j] - 1 for j in range(len(finalOut))]
    Out.sort()
    return(Out)

#A weird method that returns a different output for each element of equality that x could be equal to.
#The author assumes the length of equality and output are the same and that x is actually one of the elements of equality.
def multiIf(x, equality, output):
    for i in range(len(equality)):
        if x == equality[i]:
            return(output[i])

#This method translates "n" to a factorial basis vector.  "fact" is the first factorial bigger than n.
def factorialBase(n,fact):
    downiter = fact-1
    def downcycle(remainder,downiter):
        innerOut = []
        newRemainder = remainder%factorial(downiter)
        place = (remainder - newRemainder)/factorial(downiter)
        innerOut.append(place)
        downiter -= 1
        if downiter > 1:
            newdata = downcycle(newRemainder,downiter)
            innerOut = innerOut+newdata
        else:
            innerOut.append(newRemainder)
        return(innerOut)
    return(downcycle(n,downiter))

#Provides an easy method to iterate over every choice of symmetries of a set.
#This tool redorders a set of a set size using an ordering on the symmetric group of that size.
#The iterator should be set to less than the factorial of the size of inSet.
def symmetricGroupIndex(inSet,i):
    size = len(inSet)
    Index = factorialBase(i,size)
    usefulCopy = inSet.copy()
    out = []
    for i in range(size-1):
        newValue = usefulCopy.pop(Index[i])
        out.append(newValue)
    out.append(usefulCopy[0])
    return(out)

#Provides a method to produce every way to pair off an even set of objects.  Please enter an even number.
def pairMe(numbers):
    def condense(OutList,InList):
        outData = []
        for i in range(len(InList)):
            setOutData = OutList.copy()
            setInData = InList.copy()
            setOutData.append(setInData.pop(i))
            setOutData.append(setInData.pop(0))
            if len(InList) == 3:
                setOutData.append(setInData.pop(0))
                outData.append(setOutData)
            else:
                processTierDown = condense(setOutData,setInData)
                outData.append(processTierDown)
        return(outData)
    
    def dataCollapse(data,step):
        global outDataX
        if step > 0:
            for i in range(len(data)):
                dataCollapse(data[i],step-1)
        else:
            for i in range(len(data)):
                outDataX.append(data[i])
            return(outDataX)
    
    if numbers%2 == 1:
        return("Error, please input a list of even length")
    else:
        initOutList = [0]
        initInList = [i for i in range(1,numbers)]
        outDo = condense(initOutList,initInList)
        outNumpy = np.array(outDo)
        shape = outNumpy.shape
        dim = len(shape)
        global outDataX
        outDataX = []
        dataCollapse(outNumpy,dim-1)
        outNum = outDataX
        del outDataX
        outVect = [[outNum[numbers*i+j] for j in range(numbers)] for i in range(len(outNum)/numbers)]
        return(outVect)

#Takes a set of unique roots, and returns the simplest polynomial to contain those roots.
def rootsToPoly(Roots):
    factors = [x-n for n in Roots]
    out = 1
    for i in range(len(factors)):
        out = out*factors[i]
    return(out)
    
#Takes a polynomial with unique roots and returns a list of those roots.
def polyToRoots(poly):
    Roots = poly.roots()
    out = [n[0] for n in Roots]
    out.sort()
    return(out)

#Takes an elliptic curve object of the form y^2 = f(x) and return f(x).
def functFromEC(E):
    A = E.a_invariants()
    EN = x^3+A[1]*x^2+A[3]*x+A[4]
    return(A,EN)

#Finds all fifth roots of uniquty present in a finite field.
def find5thRU(p,degree=2):
    global z5
    #p is a prime and congruent to 0 mod 5 iff p IS 5.  In that case, z5 is 2.
    if p%5 == 0:
        #Here p is 5.
        z5 = 2
        return(z5)
    else:
        #Calculate square root of 5.
        fivePoly = x^2-F(5)
        #We only get roots that work for us under these circumstances:
        testA = (p%5 == 1)
        testB = (p%5 == 4 and degree==2)
        testC = (degree>=4)
        if testA or testB or testC:
            #Define square root of 5.
            sqrt5 = fivePoly.roots()[0][0]
            #Define the golden ratio.
            goldenRat = (1+sqrt5)/2
            #z5 is a root fo the following poly.
            sqrtGRM3 = (x^2-goldenRat+3).roots()[0][0]
            z5 = F(-1/2)*(F(sqrtGRM3)+F(goldenRat))
            #Check z5 is a fifth root of unity.
            if z5^5 == 1:
                return(z5)
            else:
                return("An unexpected error has occurred.")
        else:
            raise Exception("Field degree does not match...")
            
#Currently unused in new version of code
def manualIn(a,B):
    #A manual version of the "in" command, because it sometimes doesn't work with unusual types.
    length = len(B)
    truth = False
    for i in range(length):
        if a == B[i]:
            truth = True
    return(truth)

#Currently unused in new version of code
def listRemoval(List,indices):
    Remain = List.copy()
    Removed = []
    for i in range(len(indices)):
        Remain.remove(List[indices[i]])
        Removed.append(List[indices[i]])
    return(Remain,Removed)