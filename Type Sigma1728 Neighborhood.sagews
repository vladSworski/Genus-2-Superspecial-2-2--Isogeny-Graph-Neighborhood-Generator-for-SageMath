︠93d9b2b9-8419-4fe3-acd0-fee6aecd6f5fs︠
#Here we load in the appropriate exterior files for this project.
attach("VTools.sage","LGGUtilNew.sage","LGGClass.sage","LGGIsog.sage")
#We also setup an index for quadratic splittings now.
Index = pairMe(6)
︡8eb2a768-c7d6-49d7-b133-92086a7532c2︡
︠cabba038-af12-4a46-ac6b-8a6fc6be6829s︠
#Our series of calculations begins over Q this time.
polyFieldLevel0.<x> = QQ[]
F = QQ
#We will be investigating the neighborhood graph of the type Sigma 1728 curve, (E11).  This surface consists of the product of two Elliptic Curves with j-invariant 1728.  We verify our definition is correct.
ellip1728 = EllipticCurve(QQ,[0,0,0,-1,0]); ellip1728; ellip1728.j_invariant(); funct1728=x^3-x
︡59384516-0060-442a-abe1-39243a5c65c4︡
︠442627a2-6e90-4e94-86c5-be19fc3d082bs︠
#We run a function call to calculate the first isogeny from our starting curve a Hyperelliptic Curve.
try:
    isogenies1728_0 = isogenyEH([[0,1,-1],[0,1,-1]],[funct1728,funct1728],5,True); quadFactors1728_0 = [factor(n) for n in isogenies1728_0]; quadFactors1728_0
except:
    print("0,1,-1 invalid")
#We observe that we will need two field extensions, one to include sqrt(2), and another to include i.
setupFieldLevel0.<y0> = QQ[]
interFieldLevel0.<a0> = setupFieldLevel0.quotient(y0^2-2)
setupFieldLevel1.<y1> = interFieldLevel0[]
interFieldLevel1.<a1> = setupFieldLevel1.quotient(y1^2+1)
polyFieldLevel1.<x> = interFieldLevel1[]
F = interFieldLevel1
#And now we run our calculations again:
funct1728 = polyFieldLevel1(funct1728)
isogenies1728_0 = isogenyEH([[0,1,-1],[0,1,-1]],[funct1728,funct1728],5);
isogenies1728_0[1]
polyRoots1728_0 = [x-n for n in isogenies1728_0[1]]
prodRoots1728_0 = 1
for i in range(len(polyRoots1728_0)):
    prodRoots1728_0 = prodRoots1728_0*polyRoots1728_0[i]
prodRoots1728_0
Igusa1728_0 = [symmetric2_6root(polyRoots1728_0),symmetric4_6root(polyRoots1728_0),symmetric6_6root(polyRoots1728_0),symmetric10_6root(polyRoots1728_0)]
Kohel1728_0 = IgusaToKohel(Igusa1728_0)
Kohel1728_0
#From here, we will calculate the remaining isogenies from the starting node to hyperelliptic curves.  We check again that we do not need further field extensions.
try:
    isogenies1728_1 = isogenyEH([[0,1,-1],[0,-1,1]],[funct1728,funct1728],5,True); quadFactors1728_1 = [factor(n) for n in isogenies1728_1]; quadFactors1728_1
except:
    print("0,-1,1 invalid")
try:
    isogenies1728_2 = isogenyEH([[0,1,-1],[1,0,-1]],[funct1728,funct1728],5,True); quadFactors1728_2 = [factor(n) for n in isogenies1728_2]; quadFactors1728_2
except:
    print("1,0,-1, invalid")
try:
    isogenies1728_3 = isogenyEH([[0,1,-1],[1,-1,0]],[funct1728,funct1728],5,True); quadFactors1728_3 = [factor(n) for n in isogenies1728_3]; quadFactors1728_3
except:
    print("1,-1,0 invalid")
try:
    isogenies1728_4 = isogenyEH([[0,1,-1],[-1,0,1]],[funct1728,funct1728],5,True); quadFactors1728_4 = [factor(n) for n in isogenies1728_4]; quadFactors1728_4
except:
    print("-1,0,1 invalid")
try:
    isogenies1728_5 = isogenyEH([[0,1,-1],[-1,1,0]],[funct1728,funct1728],5,True); quadFactors1728_5 = [factor(n) for n in isogenies1728_5]; quadFactors1728_5
except:
    print("-1,1,0 invalid")
#Observation:  Everything is factorable.  Further,
#of six possible isogenies, only 4 are valid.  This is to be expected for type sigma curves.  The ones that are valid are as follows:
#
# 1) [0,1,-1] returns (x - a0) * (x + a0) * (x - a1) * (x + a1) * (x - 1/2*a0) * (x + 1/2*a0)
# 2) [0,-1,1] returns (x - a0*a1) * (x + a0*a1) * (x - 1) * (x + 1) * (x - 1/2*a0*a1) * (x + 1/2*a0*a1)
# 3) [1,0,-1] returns (x - a0*a1) * (x + a0*a1) * (x - 1/2*a0*a1) * (x + 1/2*a0*a1) * (x - 1) * (x + 1)
# 4) [-1,0,1] returns (x - a0) * (x + a0) * (x - 1/2*a0) * (x + 1/2*a0) * (x - a1) * (x + a1)
#
# We now run the remainder of our calculations to identify unique curves, etc.  We note the invalidity occurs for index 3 and 5.
# P.S. Florit and Smith denote the two failed isogenies as a weight two self-loop separate from the weight one self loop discussed in the next section.
print("")
isogenies1728_1 = isogenyEH([[0,1,-1],[0,-1,1]],[funct1728,funct1728],5);
isogenies1728_1[1]
polyRoots1728_1 = [x-n for n in isogenies1728_1[1]]
prodRoots1728_1 = 1
for i in range(len(polyRoots1728_1)):
    prodRoots1728_1 = prodRoots1728_1*polyRoots1728_1[i]
prodRoots1728_1
Igusa1728_1 = [symmetric2_6root(polyRoots1728_1),symmetric4_6root(polyRoots1728_1),symmetric6_6root(polyRoots1728_1),symmetric10_6root(polyRoots1728_1)]
Kohel1728_1 = IgusaToKohel(Igusa1728_1)
Kohel1728_1
print("")
isogenies1728_2 = isogenyEH([[0,1,-1],[1,0,-1]],[funct1728,funct1728],5);
isogenies1728_2[1]
polyRoots1728_2 = [x-n for n in isogenies1728_2[1]]
prodRoots1728_2 = 1
for i in range(len(polyRoots1728_2)):
    prodRoots1728_2 = prodRoots1728_2*polyRoots1728_2[i]
prodRoots1728_2
Igusa1728_2 = [symmetric2_6root(polyRoots1728_2),symmetric4_6root(polyRoots1728_2),symmetric6_6root(polyRoots1728_2),symmetric10_6root(polyRoots1728_2)]
Kohel1728_2 = IgusaToKohel(Igusa1728_2)
Kohel1728_2
print("")
isogenies1728_3 = isogenyEH([[0,1,-1],[-1,0,1]],[funct1728,funct1728],5);
isogenies1728_3[1]
polyRoots1728_3 = [x-n for n in isogenies1728_3[1]]
prodRoots1728_3 = 1
for i in range(len(polyRoots1728_3)):
    prodRoots1728_3 = prodRoots1728_3*polyRoots1728_3[i]
prodRoots1728_3
Igusa1728_3 = [symmetric2_6root(polyRoots1728_3),symmetric4_6root(polyRoots1728_3),symmetric6_6root(polyRoots1728_3),symmetric10_6root(polyRoots1728_3)]
Kohel1728_3 = IgusaToKohel(Igusa1728_3)
Kohel1728_3
#I am aware that the above code could be massively reduced in size.
#
#Observations:
#All four curves have the same invariants and hence are the same curve.  There is thus a degree 4 edge from (1728,1728) to (751/8, 778688/27, 3178232/81).  Moving forward we will use x^6 - 3/2*x^4 - 3/2*x^2 + 1 to represent this curve.
#We lastly observe the automorphism type of this curve:
automorphismDetector([n for n in specifiedReverseRosenhain6LFT(a0, -a0, a1, -a1, 1/2*a0, -1/2*a0)])
#This curve is a type 3 curve as it has an automorphism group of "K4".  In addition, the Orbit-Stabilizer theorem tells us that the reverse edge must have a weight of 1.
︡8d8de454-1904-48d5-87bf-25f97cc7e6db︡
︠cd8d6f05-316d-445c-a92a-c8ae14431d5bs︠
#Next, we need to calculate the elliptic neighbors of (1728,1728)
ellip1728_0 = ellip1728
im1728_1 = ellip1728_0.isogeny(ellip1728_0(0,0))
ellip1728_1 = im1728_1.codomain(); ellip1728_1; ellip1728_1.j_invariant(); funct1728_1 = functFromEC(ellip1728_1); roots1728_1 = funct1728_1[1].roots(); roots1728_1
im1728_2 = ellip1728_0.isogeny(ellip1728_0(1,0)); ellip1728_2 = im1728_2.codomain(); ellip1728_2; ellip1728_2.j_invariant(); funct1728_2 = functFromEC(ellip1728_2); roots1728_2 = funct1728_2[1].roots(); roots1728_2
im1728_3 = ellip1728_0.isogeny(ellip1728_0(-1,0)); ellip1728_3 = im1728_3.codomain(); ellip1728_3; ellip1728_3.j_invariant(); funct1728_3 = functFromEC(ellip1728_3); roots1728_3 = funct1728_3[1].roots(); roots1728_3
#Observations:
#
#Two of these curves are the same curve: j = 287496.  We will represent this curve with the equation y^2 = x^3 - 11*x - 14  The other is a self loop: j = 1728.
#Combinatorially, we have 3 outgoing edges to elliptic products:
#We have a weight 4 edge to a curve of type Sigma (E22) with invariants: (287496, 287496).  Since a type sigma curve has autormorphism group of order 4, we can conclude the dual edge has weight 1.
#We have a weight 4 edge to a curve of type Pi 1728 (E12) with invariants (1728, 287496).  Since a type Pi 1728 curve has automorphism group of order 4, we can conclude the dual edge has weight 1.
#We have a weight 1 self loop.
#
#At this point we have drawn the entirety of the radius 1 neighborhood.  We have 3 neighboring curves to consider.
︡b0a9d986-e1ca-41fa-9498-2e200b3a9a86︡
︠dac180dc-7738-4768-85f5-e4b0181d9a49s︠
#Investigating the type 3 curve.
#
#We begin by calculating the Hyperelliptic Neighbors of the type 3 curve.  By Florit and Smith's Atlas, we expect 13 of these: one of which should be a self loop.

#Calculating Richelot Isogenies relies on a Lie Bracket, we define that here so we can call it quickly whenever.  Note that the "d" function takes derivatives and the "Lie" function calculates the Lie Bracket for a quadratic pairing.
def d(poly):
    return(diff(poly, x))

def Lie(A,B):
    return([B[1]*A[2]-B[2]*A[1],B[2]*A[0]-B[0]*A[2],B[0]*A[1]-B[1]*A[0]])
︡d7b02893-b6f2-42a5-9317-4859285c95f9︡
︠010c247d-3de2-4108-ba71-68920d0ed4ads︠
#We establish a box to store all of the quadratic splittings for building kernels.
rootsT3 = [(x+a0),(x+1/2*a0),(x-1/2*a0),(x-a0),(x+a1),(x-a1)]
kernelBoxT3 = [[rootsT3[n[0]]*rootsT3[n[1]],rootsT3[n[2]]*rootsT3[n[3]],rootsT3[n[4]]*rootsT3[n[5]]] for n in Index]
determinantsT3 = [determinantDelta(n) for n in kernelBoxT3]
#We take the derivative
d_kernelBoxT3 = [[d(n) for n in m] for m in kernelBoxT3]
#And then we take the image of the Richelot Isogeny, where possible.
im_kernelBoxT3 = []
factoredBoxT3 = []
for i in range(len(kernelBoxT3)):
    if determinantsT3[i] == 0:
        im_kernelBoxT3.append(["Elliptic"])
        factoredBoxT3.append(["Elliptic"])
    else:
        im_kernelBoxT3.append(Lie(kernelBoxT3[i],d_kernelBoxT3[i]))
        factoredBoxT3.append([factor(n) for n in im_kernelBoxT3[i]])
factoredBoxT3
#Here are our observations thus far:
#
#0) An elliptic curve.  This pairing we will return to later.  Pairing code: [0,1,2,3,4,5]
#1) Does not fully factor.  Factors present are:
#     (x - a1 - a0) * (x + 1/3*a1 - 1/3*a0), (x^2 + (-8/9*a1 + 2/9*a0)*x - 2/3*a0*a1 - 2/3), (x^2 + (2/3*a1 + 1/3*a0)*x + 1/2*a0*a1 - 1/2)
#2) Does not fully factor.  Factors present are:
#     (x - 1/3*a1 - 1/3*a0) * (x + a1 - a0), (x^2 + (8/9*a1 + 2/9*a0)*x + 2/3*a0*a1 - 2/3), (x^2 + (-2/3*a1 + 1/3*a0)*x - 1/2*a0*a1 - 1/2)
#3) Fully factors.  Factors present are:
#     (x + 2*a0 - 3) * (x + 2*a0 + 3), (x - 2*a0 - 3) * (x - 2*a0 + 3), (x - a1) * (x + a1)
#4) Does not fully factor.  Factors present are:
#     (x + (-9/17*a0 - 3/17)*a1 - 2/17*a0 - 12/17) * (x + (9/17*a0 - 3/17)*a1 - 2/17*a0 + 12/17), (x^2 + (-8/11*a1 - 10/11*a0)*x - 2/11*a0*a1 + 6/11), (x^2 + (-2*a1 + a0)*x - 1/2*a0*a1 + 3/2)
#5) Does not fully factor.  Factors present are:
#     (x + (-9/17*a0 + 3/17)*a1 - 2/17*a0 + 12/17) * (x + (9/17*a0 + 3/17)*a1 - 2/17*a0 - 12/17), (x^2 + (8/11*a1 - 10/11*a0)*x + 2/11*a0*a1 + 6/11), (x^2 + (2*a1 + a0)*x + 1/2*a0*a1 + 3/2)
#6) An elliptic curve.  This pairing we will return to later.  Pairing code: [0,3,1,2,4,5]
#7) Does not fully factor.  Factors present are:
#     (x^2 - 1/2*a0*a1), (x + (-a0 + 1)*a1 - a0 + 1) * (x + (a0 + 1)*a1 - a0 - 1), (x + (-a0 - 1)*a1 + a0 + 1) * (x + (a0 - 1)*a1 + a0 - 1)
#8) Does not fully factor.  Factors present are:
#     (x^2 + 1/2*a0*a1), (x + (-a0 - 1)*a1 - a0 - 1) * (x + (a0 - 1)*a1 - a0 + 1), (x + (-a0 + 1)*a1 + a0 - 1) * (x + (a0 + 1)*a1 + a0 + 1)
#9) Does not fully factor.  Factors present are:
#     (x + (-1/2*a0 - 1/2)*a1 - 1/2*a0 - 1/2) * (x + (1/2*a0 - 1/2)*a1 - 1/2*a0 + 1/2), (x^2 - a0*a1), (x + (-1/2*a0 + 1/2)*a1 + 1/2*a0 - 1/2) * (x + (1/2*a0 + 1/2)*a1 + 1/2*a0 + 1/2)
#10) Does not fully factor.  Factors present are:
#     (x^2 + (2*a1 - a0)*x - 1/2*a0*a1 + 3/2), (x + (-9/17*a0 + 3/17)*a1 + 2/17*a0 - 12/17) * (x + (9/17*a0 + 3/17)*a1 + 2/17*a0 + 12/17), (x^2 + (8/11*a1 + 10/11*a0)*x - 2/11*a0*a1 + 6/11)
#11) Does not fully factor.  Factors present are:
#     (x^2 + (-2/3*a1 - 1/3*a0)*x + 1/2*a0*a1 - 1/2), (a1 + 5/2*a0) * (x^2 + (8/9*a1 - 2/9*a0)*x - 2/3*a0*a1 - 2/3), (-2*a1 - 1/2*a0) * (x - 1/3*a1 + 1/3*a0) * (x + a1 + a0)
#12) Does not fully factor.  Factors present are:
#     (x + (-1/2*a0 + 1/2)*a1 - 1/2*a0 + 1/2) * (x + (1/2*a0 + 1/2)*a1 - 1/2*a0 - 1/2), (x^2 + a0*a1), (x + (-1/2*a0 - 1/2)*a1 + 1/2*a0 + 1/2) * (x + (1/2*a0 - 1/2)*a1 + 1/2*a0 - 1/2)
#13) Does not fully factor.  Factors present are:
#     (x^2 + (-2*a1 - a0)*x + 1/2*a0*a1 + 3/2), (x + (-9/17*a0 - 3/17)*a1 + 2/17*a0 + 12/17) * (x + (9/17*a0 - 3/17)*a1 + 2/17*a0 - 12/17), (x^2 + (-8/11*a1 + 10/11*a0)*x + 2/11*a0*a1 + 6/11)
#14) Does not fully factor.  Factors present are:
#     (x^2 + (2/3*a1 - 1/3*a0)*x - 1/2*a0*a1 - 1/2), (x^2 + (-8/9*a1 - 2/9*a0)*x + 2/3*a0*a1 - 2/3), (x - a1 + a0) * (x + 1/3*a1 + 1/3*a0)
︡e17004de-998c-4af8-9729-96f3d8188719︡
︠c0ca017e-a2a4-4f58-b9a5-6238741ef2c1s︠
#An easy extensions is the fourth root of 2.  We extend our field to include this.
setupFieldLevel2.<y2> = interFieldLevel1[]
interFieldLevel2.<a2> = setupFieldLevel2.quotient(y2^2-a0)
polyFieldLevel2.<x> = interFieldLevel2[]
F = interFieldLevel2
#And we retest our factorizations.
rootsT3new = [polyFieldLevel2(x+a0),polyFieldLevel2(x+1/2*a0),polyFieldLevel2(x-1/2*a0),polyFieldLevel2(x-a0),polyFieldLevel2(x+a1),polyFieldLevel2(x-a1)]
kernelBoxT3new = [[rootsT3new[n[0]]*rootsT3new[n[1]],rootsT3new[n[2]]*rootsT3new[n[3]],rootsT3new[n[4]]*rootsT3new[n[5]]] for n in Index]
#We take the derivative
d_kernelBoxT3new = [[d(n) for n in m] for m in kernelBoxT3new]
im_kernelBoxT3new = []
factoredBoxT3new = []
for i in range(len(kernelBoxT3new)):
    if determinantsT3[i] == 0:
        im_kernelBoxT3new.append(["Elliptic"])
        factoredBoxT3new.append(["Elliptic"])
    else:
        im_kernelBoxT3new.append(Lie(kernelBoxT3new[i],d_kernelBoxT3new[i]))
        factoredBoxT3new.append([factor(n) for n in im_kernelBoxT3new[i]])
factoredBoxT3new
#Everything fully factors within the field of sqrt(2), (2)^(1/4) and i.
︡399246ad-a2f1-41c3-b3df-39f21641dbda︡
︠137b5ef6-1444-498a-b312-dd4eaa92d94bs︠
#Now we find the Igusa invariants of each curve.
invarsT3Box = []
for i in range(len(kernelBoxT3new)):
    if determinantsT3[i] == 0:
        invarsT3Box.append(["Elliptic"])
    else:
        invarsT3Box.append(HyperellipticCurve(expand(factoredBoxT3new[i][0]*factoredBoxT3new[i][1]*factoredBoxT3new[i][2])).absolute_igusa_invariants_kohel())
uniqueIgusasT3 = uniqueSetCounter(invarsT3Box)
uniqueIgusasT3
#Here are our observations:
#
#There are two Elliptic outgoing edges as stated before.
#There is an outgoing edge of weight 2 with invariants:
#(277864546323/234256*a0*a1 - 1437440713439/468512, 74992314819190916/1929229929*a0*a1 - 987617438131435438/5787689787, 49062146505649193/5787689787*a0*a1 - 1432983554086604327/34726138722)
#There is an outgoing edge of weight 2 with invariants:
#(-277864546323/234256*a0*a1 - 1437440713439/468512, -74992314819190916/1929229929*a0*a1 - 987617438131435438/5787689787, -49062146505649193/5787689787*a0*a1 - 1432983554086604327/34726138722)
#There is an outgoing edge of weight 1 with invariants:
#(751/8, 778688/27, 3178232/81)
#This is a self loop.
#There is an outgoing edge of weight 4 with invariants:
#(14774210551/76832, 12374781614502014/155649627, 27205955244928399/933897762)
#There is an outgoing edge of weight 4 with invariants:
#(-17197, -14848000/27, -7590400/81)
#
#One of these two outgoing weight 4 edges is unexpected and is the merging of 2 weight 2 edges in Florit and Smith's graph for the generic type 3 node.  Everything else is as expected thus far.
︡ba265b72-914e-49c4-ab99-0a73ec2ba194︡
︠29bd6ac4-9bc6-4327-9d87-f378186f3ff6s︠
#Next we seek to calculate the two elliptic nodes from the type 3 curve.  These curves appear to be dependent representatives, but as it would turn out - the determinant map is 0 regardless of choice and hence these curves truly are elliptic nodes.  We begin by finding representatives that don't cause errors however.
XF = polyFieldLevel2
#
#The Following section is directly copied from LGGIsog.  Sadly, due to the generalized nature of this proof, some lines needed to be tweaked and hence it was easier to do this from here.  See the LGGIsog file for a full description of what is happening here.
#
quadPairing1 = kernelBoxT3new[0]
quadMiddleCo = [quadPairing1[i].coefficients(sparse=False)[1] for i in range(3)]
quadLastCo = [quadPairing1[i].coefficients(sparse=False)[0] for i in range(3)]
g12 = quadMiddleCo[0]; g22 = quadMiddleCo[1]; g32 = quadMiddleCo[2]
g13 = quadLastCo[0]; g23 = quadLastCo[1]; g33 = quadLastCo[2]
if (g12 == g22)&(g22 == g32):
    sys.exit('Entered degenerate case. Cannot compute.')
AF.<a> = XF[]
dg = (g12+a*g22)^2-4*(a+1)*(g13+a*g23)
print(dg)
#Normally a call of roots = dg.roots() is made here, but due to errors, we consruct dg manually from the output.  The polynomial returned is essentially a^2 -34a +1.  This polynomial does not need any further extensions to be factored.
roots = [[17+12*a0,1],[17-12*a0,1]]
cA = [roots[i][0] for i in range(2)]
cX = [quadPairing1[0] + cA[i]*quadPairing1[1] for i in range(2)]
cP = [cX[i]/cX[i].coefficients(sparse=False)[2] for i in range(2)]
s1 = cP[0].roots()[0][0]; s2 = cP[1].roots()[0][0]
a11 = (g12 + 2*s2)/(2*s2-2*s1); a12 = 1 - a11
a21 = (g22 + 2*s2)/(2*s2-2*s1); a22 = 1 - a21
a31 = (g32 + 2*s2)/(2*s2-2*s1); a32 = 1 - a31
check = [expand(a11*(x-s1)^2+a12*(x-s2)^2),expand(a21*(x-s1)^2+a22*(x-s2)^2),\
         expand(a31*(x-s1)^2+a32*(x-s2)^2)]
if (check[0] == quadPairing1[0])&(check[1] == quadPairing1[1])&(check[2] == quadPairing1[2]):
    co = [a11,a12,a21,a22,a31,a32]
else:
    sys.exit('Critical failure, contact the author.')
f1 = expand((co[0]*x + co[1])*(co[2]*x+co[3])*(co[4]*x+co[5]))/(co[0]*co[2]*co[4])
f2 = expand((co[1]*x + co[0])*(co[3]*x+co[2])*(co[5]*x+co[4]))/(co[1]*co[3]*co[5])
f1; f2
ellipticOutputT3_1 = [polyToRoots(f1),polyToRoots(f2)]
ellipticOutputT3_1
ellipticOutputT3_1_invariants = [EllipticCurve([0,-33,0,-33,1]).j_invariant(), EllipticCurve([0,-33,0,-33,1]).j_invariant()]
ellipticOutputT3_1_invariants
#Observations
#This edge goes to the surface [287496,287496].  Note that we only have two edges going to elliptic curves.  Since one of them must go back to [1728,1728] that should be the other one.  Hence we don't need to calculate it.  Since the node [287496,287496] is a type sigma node, the dual edge should have weight 1 as well.  Florit and Smith *do* recognize this edge in their graph.
︡c567f14f-f6e2-46af-a707-1040766b5761︡
︠010e50e4-1949-49ba-8fad-020cc2731675s︠
#Next we return to the Pi 1728 and Sigma nodes at a distance of one from the origin.  As these nodes are dependent on Elliptic curves 1728, and 287496, we begin by finding the isogenies of curve 287496.
roots287496_0 = [n[0] for n in roots1728_2]; roots287496_0
ellip287496_0 = EllipticCurve(interFieldLevel2,[0,0,0,-11,-14])
im287496_1 = ellip287496_0.isogeny(x-roots287496_0[0]); ellip287496_1 = im287496_1.codomain(); ellip287496_1; ellip287496_1.j_invariant(); funct287496_1 = functFromEC(ellip287496_1); roots287496_1 = funct287496_1[1].roots(); roots287496_1
im287496_2 = ellip287496_0.isogeny(x-roots287496_0[1]);
ellip287496_2 = im287496_2.codomain(); ellip287496_2; ellip287496_2.j_invariant(); funct287496_2 = functFromEC(ellip287496_2); roots287496_2 = funct287496_2[1].roots(); roots287496_2
im287496_3 = ellip287496_0.isogeny(x-roots287496_0[2]); ellip287496_3 = im287496_3.codomain(); ellip287496_3; ellip287496_3.j_invariant(); funct287496_3 = functFromEC(ellip287496_3); roots287496_3 = funct287496_3[1].roots(); roots287496_3
#Observations:
#
#As we would expect, one of these isogenies goes back to the curve j = 1728.  The other two go to j = 29071392966*a0 + 41113158120 and j = -29071392966*a0 + 41113158120.
#We will represent these two curves with the equations y^2 = x^3 + (-60*a0-91)*x + (-308*a0-462) and y^2 = x^3 + (60*a0-91)*x + (308*a0-462) respectively.
#
#Combinatorially from (1728,287496), we have 6 outgoing edges to elliptic products:
#We have a weight 1 edge to a curve of type Sigma 1728 (E11) with invariants: (1728, 1728).  This is our dual edge back to the origin.
#We have a weight 2 self-loop.
#We have a weight 1 edge to a curve of type Pi 1728 (E12) with invariants: (1728,29071392966*a0 + 41113158120).  Since a type Pi 1728 curve has automorphism group of order 4, we can conclude the dual edge has weight 1.
#We have a weight 1 edge to a curve of type Pi 1728 (E12) with invariants: (1728,-29071392966*a0 + 41113158120).  Since a type Pi 1728 curve has automorphism group of order 4, we can conclude the dual edge has weight 1.
#We have a weight 2 edge to a curve of type Pi (E23) with invariants: (287496,29071392966*a0 + 41113158120).  Since a type Pi curve has automorphism group of order 2, we can conclude the dual edge has weight 1.
#We have a weight 2 edge to a curve of type Pi (E23) with invariants: (287496,-29071392966*a0 + 41113158120).  Since a type Pi curve has automorphism group of order 2, we can conclude the dual edge has weight 1.
#
#
#Combinatorially from (287496,287496), we have 6 outgoing edges to elliptic products:
#We have a weight 1 edge to a curve of type Sigma 1728 (E11) with invariants: (1728, 1728).  This is our dual edge back to the origin.
#We have a weight 2 edge to a curve of type Pi 1728 (E12) with invariants: (1728,29071392966*a0 + 41113158120).  Since a type Pi 1728 curve has automorphism group of order 4, we can conclude the dual edge has weight 2.
#We have a weight 2 edge to a curve of type Pi 1728 (E12) with invariants: (1728,-29071392966*a0 + 41113158120).  Since a type Pi 1728 curve has automorphism group of order 4, we can conclude the dual edge has weight 2.
#We have a weight 1 edge to a curve of type Sigma (E22) with invariants: (29071392966*a0 + 41113158120,29071392966*a0 + 41113158120). Since a type Sigma curve has automorphism group of order 4, we can conclude the dual edge has weight 1.
#We have a weight 1 edge to a curve of type Sigma (E22) with invariants: (-29071392966*a0 + 41113158120,-29071392966*a0 + 41113158120). Since a type Sigma curve has automorphism group of order 4, we can conclude the dual edge has weight 1.
#We have a weight 2 edge to a curve of type Pi (E23) with invariants: (29071392966*a0 + 41113158120,-29071392966*a0 + 41113158120). Since a type Pi curve has automorphism group of order 2, we can conclude the dual edge has weight 1.
#
#Another observation is that the nodes: (1728,29071392966*a0 + 41113158120), and (1728,-29071392966*a0 + 41113158120) are common neighbors to both (1728, 287496) and (287496,287496).
#
#All of this matches our expectations from Florit and Smith of the neighborhoods of type Sigma and type Pi 1728 curves, with a moderate amount of specification.
︡25a3adc9-9af9-4446-bcf3-5d5354194aa7︡
︠7c85d45d-f1d1-43ec-a823-4669530ecd96s︠
#Next we calculate the Hyperelliptic Neighbors for (287496,287496)
try:
    isogenies287496_0 = isogenyEH([[-2, 2*a0 + 1, -2*a0 + 1],[-2, 2*a0 + 1, -2*a0 + 1]],[x^3 - 11*x - 14,x^3 - 11*x - 14],5,True); quadFactors287496_0 = [factor(n) for n in isogenies287496_0]; quadFactors287496_0
except:
    print("-2, 2*a0 + 1, -2*a0 + 1 invalid")
try:
    isogenies287496_1 = isogenyEH([[-2, 2*a0 + 1, -2*a0 + 1],[-2, -2*a0 + 1, 2*a0 + 1]],[x^3 - 11*x - 14,x^3 - 11*x - 14],5,True); quadFactors287496_1 = [factor(n) for n in isogenies287496_1]; quadFactors287496_1
except:
    print("-2, -2*a0 + 1, 2*a0 + 1 invalid")
try:
    isogenies287496_2 = isogenyEH([[-2, 2*a0 + 1, -2*a0 + 1],[2*a0 + 1, -2, -2*a0 + 1]],[x^3 - 11*x - 14,x^3 - 11*x - 14],5,True); quadFactors287496_2 = [factor(n) for n in isogenies287496_2]; quadFactors287496_2
except:
    print("2*a0 + 1, -2, -2*a0 + 1 invalid")
try:
    isogenies287496_3 = isogenyEH([[-2, 2*a0 + 1, -2*a0 + 1],[2*a0 + 1, -2*a0 + 1, -2]],[x^3 - 11*x - 14,x^3 - 11*x - 14],5,True); quadFactors287496_3 = [factor(n) for n in isogenies287496_3]; quadFactors287496_3
except:
    print("2*a0 + 1, -2*a0 + 1, -2 invalid")
try:
    isogenies287496_4 = isogenyEH([[-2, 2*a0 + 1, -2*a0 + 1],[-2*a0 + 1, -2, 2*a0 + 1]],[x^3 - 11*x - 14,x^3 - 11*x - 14],5,True); quadFactors287496_4 = [factor(n) for n in isogenies287496_4]; quadFactors287496_4
except:
    print("-2*a0 + 1, -2, 2*a0 + 1 invalid")
try:
    isogenies287496_5 = isogenyEH([[-2, 2*a0 + 1, -2*a0 + 1],[-2*a0 + 1, 2*a0 + 1, -2]],[x^3 - 11*x - 14,x^3 - 11*x - 14],5,True); quadFactors287496_5 = [factor(n) for n in isogenies287496_5]; quadFactors287496_5
except:
    print("-2*a0 + 1, 2*a0 + 1, -2 invalid")
#Observation:  Everything is factorable.  Further,
#of six possible isogenies, only 5 are valid.  The ones that are valid are as follows:
#
# 1) [-2, 2*a0 + 1, -2*a0 + 1] returns (x + (-2*a0 + 2)*a2) * (x + (2*a0 - 2)*a2) * (x - a1) * (x + a1) * (x + (-1/4*a0 - 1/2)*a2) * (x + (1/4*a0 + 1/2)*a2)
# 2) [-2, -2*a0 + 1, 2*a0 + 1] returns (x + (-2*a0 - 2)*a2) * (x + (2*a0 + 2)*a2) * (x + (-2*a0 + 3)*a1) * (x + (2*a0 - 3)*a1) * (x + ((-1/4*a0 - 1/2)*a1)*a2) * (x + ((1/4*a0 + 1/2)*a1)*a2)
# 3) [2*a0 + 1, -2, -2*a0 + 1] returns (x + ((-2*a0 + 2)*a1)*a2) * (x + ((2*a0 - 2)*a1)*a2) * (x + (-1/4*a0 + 1/2)*a2) * (x + (1/4*a0 - 1/2)*a2) * (x + (-2*a0 - 3)*a1) * (x + (2*a0 + 3)*a1)
# 4) [2*a0 + 1, -2*a0 + 1, -2] returns (x - a1) * (x + a1) * (x - 2*a0 + 3) * (x + 2*a0 - 3) * (x - 2*a0 - 3) * (x + 2*a0 + 3)
# 5) [-2*a0 + 1, -2, 2*a0 + 1] returns (x + ((-2*a0 - 2)*a1)*a2) * (x + ((2*a0 + 2)*a1)*a2) * (x + ((-1/4*a0 + 1/2)*a1)*a2) * (x + ((1/4*a0 - 1/2)*a1)*a2) * (x - a1) * (x + a1)
#
# We now run the remainder of our calculations to identify unique curves, etc.  We note the invalidity occurs for index 5
print("")
isogenies287496_0 = isogenyEH([[-2, 2*a0 + 1, -2*a0 + 1],[-2, 2*a0 + 1, -2*a0 + 1]],[x^3 - 11*x - 14,x^3 - 11*x - 14],5);
isogenies287496_0[1]
polyRoots287496_0 = [x-n for n in isogenies287496_0[1]]
prodRoots287496_0 = 1
for i in range(len(polyRoots287496_0)):
    prodRoots287496_0 = prodRoots287496_0*polyRoots287496_0[i]
prodRoots287496_0
Igusa287496_0 = [symmetric2_6root(polyRoots287496_0),symmetric4_6root(polyRoots287496_0),symmetric6_6root(polyRoots287496_0),symmetric10_6root(polyRoots287496_0)]
Kohel287496_0 = IgusaToKohel(Igusa287496_0)
Kohel287496_0
print("")
isogenies287496_1 = isogenyEH([[-2, 2*a0 + 1, -2*a0 + 1],[-2, -2*a0 + 1, 2*a0 + 1]],[x^3 - 11*x - 14,x^3 - 11*x - 14],5);
isogenies287496_1[1]
polyRoots287496_1 = [x-n for n in isogenies287496_1[1]]
prodRoots287496_1 = 1
for i in range(len(polyRoots287496_1)):
    prodRoots287496_1 = prodRoots287496_1*polyRoots287496_1[i]
prodRoots287496_1
Igusa287496_1 = [symmetric2_6root(polyRoots287496_1),symmetric4_6root(polyRoots287496_1),symmetric6_6root(polyRoots287496_1),symmetric10_6root(polyRoots287496_1)]
Kohel287496_1 = IgusaToKohel(Igusa287496_1)
Kohel287496_1
print("")
isogenies287496_2 = isogenyEH([[-2, 2*a0 + 1, -2*a0 + 1],[2*a0 + 1, -2, -2*a0 + 1]],[x^3 - 11*x - 14,x^3 - 11*x - 14],5);
isogenies287496_2[1]
polyRoots287496_2 = [x-n for n in isogenies287496_2[1]]
prodRoots287496_2 = 1
for i in range(len(polyRoots287496_2)):
    prodRoots287496_2 = prodRoots287496_2*polyRoots287496_2[i]
prodRoots287496_2
Igusa287496_2 = [symmetric2_6root(polyRoots287496_2),symmetric4_6root(polyRoots287496_2),symmetric6_6root(polyRoots287496_2),symmetric10_6root(polyRoots287496_2)]
Kohel287496_2 = IgusaToKohel(Igusa287496_2)
Kohel287496_2
print("")
isogenies287496_3 = isogenyEH([[-2, 2*a0 + 1, -2*a0 + 1],[2*a0 + 1, -2*a0 + 1, -2]],[x^3 - 11*x - 14,x^3 - 11*x - 14],5);
isogenies287496_3[1]
polyRoots287496_3 = [x-n for n in isogenies287496_3[1]]
prodRoots287496_3 = 1
for i in range(len(polyRoots287496_3)):
    prodRoots287496_3 = prodRoots287496_3*polyRoots287496_3[i]
prodRoots287496_3
Igusa287496_3 = [symmetric2_6root(polyRoots287496_3),symmetric4_6root(polyRoots287496_3),symmetric6_6root(polyRoots287496_3),symmetric10_6root(polyRoots287496_3)]
Kohel287496_3 = IgusaToKohel(Igusa287496_3)
Kohel287496_3
print("")
isogenies287496_4 = isogenyEH([[-2, 2*a0 + 1, -2*a0 + 1],[-2*a0 + 1, -2, 2*a0 + 1]],[x^3 - 11*x - 14,x^3 - 11*x - 14],5);
isogenies287496_4[1]
polyRoots287496_4 = [x-n for n in isogenies287496_4[1]]
prodRoots287496_4 = 1
for i in range(len(polyRoots287496_4)):
    prodRoots287496_4 = prodRoots287496_4*polyRoots287496_4[i]
prodRoots287496_4
Igusa287496_4 = [symmetric2_6root(polyRoots287496_4),symmetric4_6root(polyRoots287496_4),symmetric6_6root(polyRoots287496_4),symmetric10_6root(polyRoots287496_4)]
Kohel287496_4 = IgusaToKohel(Igusa287496_4)
Kohel287496_4

#Observations:
#
#There is an outgoing edge of weight 1 to the curve with invariants: (6523884707703/38416*a0 + 36905627140531/153664, 1504854025601701280/51883209*a0 + 6384680530861497176/155649627, 1693972751480432852/155649627*a0 + 7187037942295351235/466948881)
#There is an outgoing edge of weight 1 to the curve with invariants: (-6523884707703/38416*a0 + 36905627140531/153664, -1504854025601701280/51883209*a0 + 6384680530861497176/155649627, -1693972751480432852/155649627*a0 + 7187037942295351235/466948881)
#There is an outgoing edge of weight 2 to the curve with invariants: (-101155/64, -12888615920536/47832147, -13901603081995/143496441)
#There is an outgoing edge of weight 1 to the curve with invariants: (751/8, 778688/27, 3178232/81).  We note that this last edge is the type 3 curve we investigated earlier and this is the dual edge we expected.
︡e055af68-497c-4819-870f-b48944d34834︡
︠6524b0af-57a2-4364-be8c-affedd777fa6s︠
#Next we calculate the Hyperelliptic Neighbors for (1728,287496)
try:
    isogenies1728_287496_0 = isogenyEH([[0, 1, -1],[-2, 2*a0 + 1, -2*a0 + 1]],[x^3 - x,x^3 - 11*x - 14],5,True); quadFactors1728_287496_0 = [factor(n) for n in isogenies1728_287496_0]; quadFactors1728_287496_0
except:
    print("-2, 2*a0 + 1, -2*a0 + 1 invalid")
try:
    isogenies1728_287496_1 = isogenyEH([[0, 1, -1],[-2, -2*a0 + 1, 2*a0 + 1]],[x^3 - x,x^3 - 11*x - 14],5,True); quadFactors1728_287496_1 = [factor(n) for n in isogenies1728_287496_1]; quadFactors1728_287496_1
except:
    print("-2, -2*a0 + 1, 2*a0 + 1 invalid")
try:
    isogenies1728_287496_2 = isogenyEH([[0, 1, -1],[2*a0 + 1, -2, -2*a0 + 1]],[x^3 - x,x^3 - 11*x - 14],5,True); quadFactors1728_287496_2 = [factor(n) for n in isogenies1728_287496_2]; quadFactors1728_287496_2
except:
    print("2*a0 + 1, -2, -2*a0 + 1 invalid")
try:
    isogenies1728_287496_3 = isogenyEH([[0, 1, -1],[2*a0 + 1, -2*a0 + 1, -2]],[x^3 - x,x^3 - 11*x - 14],5,True); quadFactors1728_287496_3 = [factor(n) for n in isogenies1728_287496_3]; quadFactors1728_287496_3
except:
    print("2*a0 + 1, -2*a0 + 1, -2 invalid")
try:
    isogenies1728_287496_4 = isogenyEH([[0, 1, -1],[-2*a0 + 1, -2, 2*a0 + 1]],[x^3 - x,x^3 - 11*x - 14],5,True); quadFactors1728_287496_4 = [factor(n) for n in isogenies1728_287496_4]; quadFactors1728_287496_4
except:
    print("-2*a0 + 1, -2, 2*a0 + 1 invalid")
try:
    isogenies1728_287496_5 = isogenyEH([[0, 1, -1],[-2*a0 + 1, 2*a0 + 1, -2]],[x^3 - x,x^3 - 11*x - 14],5,True); quadFactors1728_287496_5 = [factor(n) for n in isogenies1728_287496_5]; quadFactors1728_287496_5
except:
    print("-2*a0 + 1, 2*a0 + 1, -2 invalid")
#Observation:  Everything is factorable.  Further,
#of six possible isogenies, all are valid.
#
# 1) [-2, 2*a0 + 1, -2*a0 + 1] returns (x - a0 + 2) * (x + a0 - 2) * (x - a0 - 1) * (x + a0 + 1) * (x - 1/4*a0*a2) * (x + 1/4*a0*a2)
# 2) [-2, -2*a0 + 1, 2*a0 + 1] returns (x - a0 - 2) * (x + a0 + 2) * (x - a0 + 1) * (x + a0 - 1) * (x - 1/4*a0*a1*a2) * (x + 1/4*a0*a1*a2)
# 3) [2*a0 + 1, -2, -2*a0 + 1] returns (x + (-a0 + 2)*a1) * (x + (a0 - 2)*a1) * (x - 1/4*a0*a1*a2) * (x + 1/4*a0*a1*a2) * (x + (-a0 - 1)*a1) * (x + (a0 + 1)*a1)
# 4) [2*a0 + 1, -2*a0 + 1, -2] returns (x - 1/2*a1*a2) * (x + 1/2*a1*a2) * (x + (-a0 + 1)*a1) * (x + (a0 - 1)*a1) * (x - a0 - 1) * (x + a0 + 1)
# 5) [-2*a0 + 1, -2, 2*a0 + 1] returns (x + (-a0 - 2)*a1) * (x + (a0 + 2)*a1) * (x - 1/4*a0*a2) * (x + 1/4*a0*a2) * (x + (-a0 + 1)*a1) * (x + (a0 - 1)*a1)
# 6) [-2*a0 + 1, 2*a0 + 1, -2] returns (x - 1/2*a2) * (x + 1/2*a2) * (x + (-a0 - 1)*a1) * (x + (a0 + 1)*a1) * (x - a0 + 1) * (x + a0 - 1)
#
# We now run the remainder of our calculations to identify unique curves, etc.
print("")
isogenies1728_287496_0 = isogenyEH([[0,1,-1],[-2, 2*a0 + 1, -2*a0 + 1]],[x^3 - x,x^3 - 11*x - 14],5);
isogenies1728_287496_0[1]
polyRoots1728_287496_0 = [x-n for n in isogenies1728_287496_0[1]]
prodRoots1728_287496_0 = 1
for i in range(len(polyRoots1728_287496_0)):
    prodRoots1728_287496_0 = prodRoots1728_287496_0*polyRoots1728_287496_0[i]
prodRoots1728_287496_0
Igusa1728_287496_0 = [symmetric2_6root(polyRoots1728_287496_0),symmetric4_6root(polyRoots1728_287496_0),symmetric6_6root(polyRoots1728_287496_0),symmetric10_6root(polyRoots1728_287496_0)]
Kohel1728_287496_0 = IgusaToKohel(Igusa1728_287496_0)
Kohel1728_287496_0
print("")
isogenies1728_287496_1 = isogenyEH([[0,1,-1],[-2, -2*a0 + 1, 2*a0 + 1]],[x^3 - x,x^3 - 11*x - 14],5);
isogenies1728_287496_1[1]
polyRoots1728_287496_1 = [x-n for n in isogenies1728_287496_1[1]]
prodRoots1728_287496_1 = 1
for i in range(len(polyRoots1728_287496_1)):
    prodRoots1728_287496_1 = prodRoots1728_287496_1*polyRoots1728_287496_1[i]
prodRoots1728_287496_1
Igusa1728_287496_1 = [symmetric2_6root(polyRoots1728_287496_1),symmetric4_6root(polyRoots1728_287496_1),symmetric6_6root(polyRoots1728_287496_1),symmetric10_6root(polyRoots1728_287496_1)]
Kohel1728_287496_1 = IgusaToKohel(Igusa1728_287496_1)
Kohel1728_287496_1
print("")
isogenies1728_287496_2 = isogenyEH([[0,1,-1],[2*a0 + 1, -2, -2*a0 + 1]],[x^3 - x,x^3 - 11*x - 14],5);
isogenies1728_287496_2[1]
polyRoots1728_287496_2 = [x-n for n in isogenies1728_287496_2[1]]
prodRoots1728_287496_2 = 1
for i in range(len(polyRoots1728_287496_2)):
    prodRoots1728_287496_2 = prodRoots1728_287496_2*polyRoots1728_287496_2[i]
prodRoots1728_287496_2
Igusa1728_287496_2 = [symmetric2_6root(polyRoots1728_287496_2),symmetric4_6root(polyRoots1728_287496_2),symmetric6_6root(polyRoots1728_287496_2),symmetric10_6root(polyRoots1728_287496_2)]
Kohel1728_287496_2 = IgusaToKohel(Igusa1728_287496_2)
Kohel1728_287496_2
print("")
isogenies1728_287496_3 = isogenyEH([[0,1,-1],[2*a0 + 1, -2*a0 + 1, -2]],[x^3 - x,x^3 - 11*x - 14],5);
isogenies1728_287496_3[1]
polyRoots1728_287496_3 = [x-n for n in isogenies1728_287496_3[1]]
prodRoots1728_287496_3 = 1
for i in range(len(polyRoots1728_287496_3)):
    prodRoots1728_287496_3 = prodRoots1728_287496_3*polyRoots1728_287496_3[i]
prodRoots1728_287496_3
Igusa1728_287496_3 = [symmetric2_6root(polyRoots1728_287496_3),symmetric4_6root(polyRoots1728_287496_3),symmetric6_6root(polyRoots1728_287496_3),symmetric10_6root(polyRoots1728_287496_3)]
Kohel1728_287496_3 = IgusaToKohel(Igusa1728_287496_3)
Kohel1728_287496_3
print("")
isogenies1728_287496_4 = isogenyEH([[0,1,-1],[-2*a0 + 1, -2, 2*a0 + 1]],[x^3 - x,x^3 - 11*x - 14],5);
isogenies1728_287496_4[1]
polyRoots1728_287496_4 = [x-n for n in isogenies1728_287496_4[1]]
prodRoots1728_287496_4 = 1
for i in range(len(polyRoots1728_287496_4)):
    prodRoots1728_287496_4 = prodRoots1728_287496_4*polyRoots1728_287496_4[i]
prodRoots1728_287496_4
Igusa1728_287496_4 = [symmetric2_6root(polyRoots1728_287496_4),symmetric4_6root(polyRoots1728_287496_4),symmetric6_6root(polyRoots1728_287496_4),symmetric10_6root(polyRoots1728_287496_4)]
Kohel1728_287496_4 = IgusaToKohel(Igusa1728_287496_4)
Kohel1728_287496_4
print("")
isogenies1728_287496_5 = isogenyEH([[0,1,-1],[-2*a0 + 1, 2*a0 + 1, -2]],[x^3 - x,x^3 - 11*x - 14],5);
isogenies1728_287496_5[1]
polyRoots1728_287496_5 = [x-n for n in isogenies1728_287496_5[1]]
prodRoots1728_287496_5 = 1
for i in range(len(polyRoots1728_287496_5)):
    prodRoots1728_287496_5 = prodRoots1728_287496_5*polyRoots1728_287496_5[i]
prodRoots1728_287496_5
Igusa1728_287496_5 = [symmetric2_6root(polyRoots1728_287496_5),symmetric4_6root(polyRoots1728_287496_5),symmetric6_6root(polyRoots1728_287496_5),symmetric10_6root(polyRoots1728_287496_5)]
Kohel1728_287496_5 = IgusaToKohel(Igusa1728_287496_5)
Kohel1728_287496_5
#Observations:
#
#There is an outgoing edge of weight 2 to the curve with invariants: (3118654011/76832*a0 + 1103115731/19208, 291101769555344/51883209*a0 + 1235001447887552/155649627, 250932640456322/155649627*a0 + 1064612285857304/466948881)
#There is an outgoing edge of weight 2 to the curve with invariants: (-3118654011/76832*a0 + 1103115731/19208, -291101769555344/51883209*a0 + 1235001447887552/155649627, -250932640456322/155649627*a0 + 1064612285857304/466948881)
#There is an outgoing edge of weight 2 to the curve with invariants: (-17197, -14848000/27, -7590400/81)  (We note that this is also a neighbor of the type 3 node.)
︡f330a2bb-cf54-4fca-9c28-e451cab2f8bc︡
︠ac055b4a-c98f-4c7b-a30b-f9236849bf83s︠
#Finally, we check if any of Hyperelliptic edges we have found are common to (287496,287496), (1728, 287496) and/or (751/8, 778688/27, 3178232/81).
Kohels287496 = [n[0] for n in uniqueSetCounter([Kohel287496_0,Kohel287496_1,Kohel287496_2,Kohel287496_3,Kohel287496_4])]
Kohels1728_287496 = [n[0] for n in uniqueSetCounter([Kohel1728_287496_0,Kohel1728_287496_1,Kohel1728_287496_2,Kohel1728_287496_3,Kohel1728_287496_4,Kohel1728_287496_5])]
KohelsT3 = [n[0] for n in uniqueSetCounter([n[0] for n in uniqueIgusasT3][1:])] #Index 1 onwards removes the Elliptic Label
#We now ensure for checking purposes that all of the invariants are in the same field.
Kohels287496Re = []
for i in range(len(Kohels287496)):
    Kohels287496Re.append([])
    for j in range(3):
        Kohels287496Re[i].append(str(Kohels287496[i][j]))
Kohels1728_287496Re = []
for i in range(len(Kohels1728_287496)):
    Kohels1728_287496Re.append([])
    for j in range(3):
        Kohels1728_287496Re[i].append(str(Kohels1728_287496[i][j]))
KohelsT3Re = []
for i in range(len(KohelsT3)):
    KohelsT3Re.append([])
    for j in range(3):
        KohelsT3Re[i].append(str(KohelsT3[i][j]))
#We can now check for common nodes.
KohelsAB = Kohels287496Re + Kohels1728_287496Re
KohelsAC = Kohels287496Re + KohelsT3Re
KohelsBC = Kohels1728_287496Re + KohelsT3Re
AllKohels = Kohels287496Re + Kohels1728_287496Re + KohelsT3Re
uniqueAB = uniqueSetCounter(KohelsAB)
uniqueAC = uniqueSetCounter(KohelsAC)
uniqueBC = uniqueSetCounter(KohelsBC)
uniqueABC = uniqueSetCounter(AllKohels)
print("Common to A and B:")
for i in range(len(uniqueAB)):
    if uniqueAB[i][1] > 1:
        uniqueAB[i]
print("Common to A and C:")
for i in range(len(uniqueAC)):
    if uniqueAC[i][1] > 1:
        uniqueAC[i]
print("Common to B and C:")
for i in range(len(uniqueBC)):
    if uniqueBC[i][1] > 1:
        uniqueBC[i]
print("Common to All:")
for i in range(len(uniqueABC)):
    if uniqueABC[i][1] > 2:
        uniqueABC[i]
#We see the nodes that appear more than once are the following:
#(751/8, 778688/27, 3178232/81).
#But this isn't new information.  This is just saying that the tier 3 node borders the tier 3 node (self-loop) and (287496,287496) which was previously established.
#(-17197, -14848000/27, -7590400/81)
︡a3cdadc2-302b-4c41-ac6c-7470e39d5037︡
︠02a899f1-fb30-4888-8e7e-81e0e476e309s︠
#We now calculate the type of the curve with the invariants (-17197, -14848000/27, -7590400/81).
aRoots = isogenies1728_287496_5[1]
automorphismDetector([n for n in specifiedReverseRosenhain6LFT(aRoots[0],aRoots[1],aRoots[2],aRoots[3],aRoots[4],aRoots[5])])
#This is a type 1 curve, and hence the dual edge to (1728,287496) has weight 1 and the dual edge to the type 3 curve has weight 2.
︡f29a1d58-da5d-4c12-80c1-4a1c8c706143︡









