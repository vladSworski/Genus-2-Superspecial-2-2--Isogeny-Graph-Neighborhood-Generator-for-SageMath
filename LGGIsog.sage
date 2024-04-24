#LGGIsog file for Local Graph Generator v1.5 (Fully commented.)
#Last updated 08/18/2023.
#Other files in our package: "LGGClass.sage", "LGGUtilNew.sage", "LGGDataGen.sage", "VTools.sage"
#Our files for data production are "LGGScript.sagews", "LGGDataCollector.sagews"
#This file primarily provides the methods for degree 2 isogenies between various kinds of curves.
#In essence, this builds the edges of our graph, but not in a literal sense.

#An isogeny between a Hyperelliptic Class and another Hyperelliptic Class.
def isogenyHH(quadPairing,determinant):
    #These isogenies are constructed from the determinant of the Wronskian of a quadratic splitting.
    H = [diff(quadPairing[(i+1)%3],x)*quadPairing[(i+2)%3] - diff(quadPairing[(i+2)%3],x)\
         *quadPairing[(i+1)%3] for i in range(3)]
    OutCurve = polyToRoots(determinant^(-1)*H[0]*H[1]*H[2])
    return(OutCurve)

#An isogeny between a Hyperelliptic Class and the Class of a product of Elliptic curves.
def isogenyHE(quadPairing):
    #These methods are taken from Benjamin Smith's thesis, where he provides a constructive proof
    #that the constants a_ij, and s_k should exist.
    #Run the algorithm to find the coefficients a_ij. Pull coefficients for processing later.
    quadMiddleCo = [quadPairing[i].coefficients(sparse=False)[1] for i in range(3)]
    quadLastCo = [quadPairing[i].coefficients(sparse=False)[0] for i in range(3)]
    #Relabeling for simplicity.
    g12 = quadMiddleCo[0]; g22 = quadMiddleCo[1]; g32 = quadMiddleCo[2]
    g13 = quadLastCo[0]; g23 = quadLastCo[1]; g33 = quadLastCo[2]
    if (g12 == g22)&(g22 == g32):
        sys.exit('Entered degenerate case. Cannot compute.')
    #Calculate the discriminant of the alpha polynomial.
    AF.<a> = XF[]
    dg = (g12+a*g22)^2-4*(a+1)*(g13+a*g23)
    #Find the roots alpha_1 and alpha_2.
    print(dg)
    roots = dg.roots()
    cA = [roots[i][0] for i in range(2)]
    #Calculate s_1 and s_2.
    cX = [quadPairing[0] + cA[i]*quadPairing[1] for i in range(2)]
    cP = [cX[i]/cX[i].coefficients(sparse=False)[2] for i in range(2)]
    s1 = cP[0].roots()[0][0]; s2 = cP[1].roots()[0][0]
    #Calculate the a_ij using the middle coefficients from earlier.
    a11 = (g12 + 2*s2)/(2*s2-2*s1); a12 = 1 - a11
    a21 = (g22 + 2*s2)/(2*s2-2*s1); a22 = 1 - a21
    a31 = (g32 + 2*s2)/(2*s2-2*s1); a32 = 1 - a31
    #Check that the constructed polynomials are in fact the original quadratic forms: g1,g2,g3.
    check = [expand(a11*(x-s1)^2+a12*(x-s2)^2),expand(a21*(x-s1)^2+a22*(x-s2)^2),\
             expand(a31*(x-s1)^2+a32*(x-s2)^2)]
    if (check[0] == quadPairing[0])&(check[1] == quadPairing[1])&(check[2] == quadPairing[2]):
        #Here co stands for coefficients
        co = [a11,a12,a21,a22,a31,a32]
    else:
        #If the algorithm failed to work, debugging is necessary.  This shouldn't happen.
        sys.exit('Critical failure, contact the author.')
    #Construct f1 and f2. Divide by leading coefficient to make monic.
    f1 = expand((co[0]*x + co[1])*(co[2]*x+co[3])*(co[4]*x+co[5]))/(co[0]*co[2]*co[4])
    f2 = expand((co[1]*x + co[0])*(co[3]*x+co[2])*(co[5]*x+co[4]))/(co[1]*co[3]*co[5])
    OutCurve = [polyToRoots(f1),polyToRoots(f2)]
    return(OutCurve)

#An isogeny between the Class of a product of Elliptic Curves and a Hyperelliptic Class.
def isogenyEH(roots,E,i,polyMode = False):
    #Relabel the roots of the two polynomials.
    R1 = roots[0]; R2Pre = roots[1]
    #Generate a symmetry of the roots of the second polynomial.
    R2 = symmetricGroupIndex(R2Pre,i)
    a1 = F(R1[0]); a2 = F(R1[1]); a3 = F(R1[2])
    b1 = F(R2[0]); b2 = F(R2[1]); b3 = F(R2[2])
    #Calculate constants for isogeny formula.
    A1 = (a3-a2)^2/(b3-b2)+(a2-a1)^2/(b2-b1)+(a1-a3)^2/(b1-b3)
    B1 = (b3-b2)^2/(a3-a2)+(b2-b1)^2/(a2-a1)+(b1-b3)^2/(a1-a3)
    A2 = a1*(b3-b2)+a2*(b1-b3)+a3*(b2-b1)
    B2 = b1*(a3-a2)+b2*(a1-a3)+b3*(a2-a1)
    #If either of the denominators is equal to zero, this isogeny does not exist.
    if A2 == 0 or B2 == 0:
        #The False here, and the True below make it easier to test when this occurs outside the method.
        return(False,0)
    #Pull out coefficients for discriminant calculation.
    C1 = E[0].coefficients(sparse=False)
    C2 = E[1].coefficients(sparse=False)
    #Calculate the discriminants.
    DA = C1[1]^2*C1[2]^2-4*C1[0]*C1[2]^3-4*C1[1]^3*C1[3]-27*C1[0]^2*C1[3]^2+18*C1[0]*C1[1]*C1[2]*C1[3]
    DB = C2[1]^2*C2[2]^2-4*C2[0]*C2[2]^3-4*C2[1]^3*C2[3]-27*C2[0]^2*C2[3]^2+18*C2[0]*C2[1]*C2[2]*C2[3]
    #Calculate the last constants needed for the isogeny formula.
    A = DB*A1/A2; B = DA*B1/B2
    #Finally calculate the polynomial for the resulting image.
    g1 = -1*(A*(a2-a1)*(a1-a3)*x^2+B*(b2-b1)*(b1-b3))
    g2 = (A*(a3-a2)*(a2-a1)*x^2+B*(b3-b2)*(b2-b1))
    g3 = (A*(a1-a3)*(a3-a2)*x^2+B*(b1-b3)*(b3-b2))
    g = g1*g2*g3
    #Calculate and return the roots of the new polynomial.
    if polyMode == True:
        return(g1,g2,g3)
    else:
        roots = polyToRoots(g)
        return(True,roots)

#Generate all the isogenies between the class of an elliptic curve, and the class of another elliptic curve.
#Note that unlike the other methods this one returns multiple isogenies, and only calculates for one half of
#any elliptic product pair.
def isogenyEE(roots,E):
    #For each possible 2-isogeny, take the isogeny root, calculate the isogeny, find the codomain, take its function, then find its roots.
    s0 = E(roots[0],0); I0 = E.isogeny(s0); E0 = I0.codomain(); fcn0 = functFromEC(E0); r0 = fcn0[1].roots()
    s1 = E(roots[1],0); I1 = E.isogeny(s1); E1 = I1.codomain(); fcn1 = functFromEC(E1); r1 = fcn1[1].roots()
    s2 = E(roots[2],0); I2 = E.isogeny(s2); E2 = I2.codomain(); fcn2 = functFromEC(E2); r2 = fcn2[1].roots()
    outRoots = [[r0[i][0] for i in range(len(r0))],[r1[i][0] for i in range(len(r1))],[r2[i][0] for i in range(len(r2))]]
    return(outRoots)

#This returns the determinant of the Wronskian of a quadratic splitting/pairing.  Represented by delta in the literature.
def determinantDelta(quadPairing):
    detMat = []
    for i in range(3):
        Row = quadPairing[i].coefficients(sparse=False)
        for j in range(3):
            detMat.append(Row[j])
    determinant = matrix(3, 3, detMat).det()
    return(determinant)
    
#This method returns all possible quadratic splittings for a set of roots, so they can be iterated over.
def quadraticPairingGenerator(roots):
    Index = pairMe(6)
    monomialBox = [[x-roots[Index[i][j]] for j in range(6)] for i in range(15)]
    polyOutBox = [[expand(monomialBox[i][2*j]*monomialBox[i][2*j+1]) for j in range(3)]\
                  for i in range(15)]
    return(polyOutBox)