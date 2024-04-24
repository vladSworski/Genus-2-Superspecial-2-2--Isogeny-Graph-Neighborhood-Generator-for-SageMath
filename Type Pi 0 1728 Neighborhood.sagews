︠394cbad7-f184-4c02-b42f-21f2b4aa11cds︠
#Here we load in the appropriate exterior files for this project.
attach("VTools.sage","LGGUtilNew.sage","LGGClass.sage","LGGIsog.sage")
#We also setup an index for quadratic splittings now.
Index = pairMe(6)
︡91e25a97-c934-4822-a9ab-a071f6f70bee︡
︠9008b5ba-170b-4d88-b949-876630c2315fs︠
#Because we are working with the curve j=0, we also need third roots of unity.  We can easily build these from the field extensions for sqrt(3) and i.  Let's establish our field to include those now.
setupFieldLevel0.<y0> = QQ[]
interFieldLevel0.<a0> = setupFieldLevel0.quotient(y0^2-3)
setupFieldLevel1.<y1> = interFieldLevel0[]
interFieldLevel1.<a1> = setupFieldLevel1.quotient(y1^2+1)
polyFieldLevel0.<x> = interFieldLevel1[]
F = interFieldLevel1
#We will also make a shortcut for ourselves to call the third root of unity:
r3 = -1/2+1/2*a0*a1
#We will be investigating the neighborhood graph of the type Pi0_1728 curve, (E01).  This surface consists of the product of two Elliptic Curves with j-invariants 0 and 1728.  We verify our definition is correct.
ellip1728 = EllipticCurve(QQ,[0,0,0,-1,0]); ellip1728; ellip1728.j_invariant(); funct1728=x^3-x
ellip0 = EllipticCurve(QQ,[0,0,0,0,-1]); ellip0; ellip0.j_invariant(); funct0 = x^3 - 1
︡81aeb11c-f41b-4acb-ae78-7441e807ef90︡
︠a5df0973-078b-44f0-b7c6-7cb3b4bd112as︠
#We run a function call to calculate the first isogeny from our starting curve to a Hyperelliptic Curve.
try:
    isogenies0_1728_0 = isogenyEH([[1,r3,r3^2],[0,1,-1]],[funct0,funct1728],5,True); quadFactors0_1728_0 = [factor(n) for n in isogenies0_1728_0]; quadFactors0_1728_0
except:
    print("0,1,-1 invalid")
#We observe that we will need another two field extensions, to include the fourth root of -3 and sqrt(a0*a1-3).
setupFieldLevel2.<y2> = interFieldLevel1[]
interFieldLevel2.<a2> = setupFieldLevel2.quotient(y2^2-a0*a1)
setupFieldLevel3.<y3> = interFieldLevel2[]
interFieldLevel3.<a3> = setupFieldLevel3.quotient(y3^2-a0*a1+3)
polyFieldLevel1.<x> = interFieldLevel3[]
F = interFieldLevel3
#And now we run our calculations again:
funct0 = polyFieldLevel1(funct0)
funct1728 = polyFieldLevel1(funct1728)
isogenies0_1728_0 = isogenyEH([[1,r3,r3^2],[0,1,-1]],[funct0,funct1728],5)
isogenies0_1728_0[1]
polyRoots0_1728_0 = [x-n for n in isogenies0_1728_0[1]]
prodRoots0_1728_0 = 1
for i in range(len(polyRoots0_1728_0)):
    prodRoots0_1728_0 = prodRoots0_1728_0*polyRoots0_1728_0[i]
prodRoots0_1728_0
Igusa0_1728_0 = [symmetric2_6root(polyRoots0_1728_0),symmetric4_6root(polyRoots0_1728_0),symmetric6_6root(polyRoots0_1728_0),symmetric10_6root(polyRoots0_1728_0)]
Kohel0_1728_0 = IgusaToKohel(Igusa0_1728_0)
Kohel0_1728_0
polyRoots0_1728_0
#From here, we will calculate the remaining isogenies from the starting node to hyperelliptic curves.  We check again that we do not need further field extensions.
try:
    isogenies0_1728_1 = isogenyEH([[1,r3,r3^2],[0,-1,1]],[funct0,funct1728],5,True); quadFactors0_1728_1 = [factor(n) for n in isogenies0_1728_1]; quadFactors0_1728_1
except:
    print("0,-1,1 invalid")
try:
    isogenies0_1728_2 = isogenyEH([[1,r3,r3^2],[1,0,-1]],[funct0,funct1728],5,True); quadFactors0_1728_2 = [factor(n) for n in isogenies0_1728_2]; quadFactors0_1728_2
except:
    print("1,0,-1, invalid")
try:
    isogenies0_1728_3 = isogenyEH([[1,r3,r3^2],[1,-1,0]],[funct0,funct1728],5,True); quadFactors0_1728_3 = [factor(n) for n in isogenies0_1728_3]; quadFactors0_1728_3
except:
    print("1,-1,0 invalid")
try:
    isogenies0_1728_4 = isogenyEH([[1,r3,r3^2],[-1,0,1]],[funct0,funct1728],5,True); quadFactors0_1728_4 = [factor(n) for n in isogenies0_1728_4]; quadFactors0_1728_4
except:
    print("-1,0,1 invalid")
try:
    isogenies0_1728_5 = isogenyEH([[1,r3,r3^2],[-1,1,0]],[funct0,funct1728],5,True); quadFactors0_1728_5 = [factor(n) for n in isogenies0_1728_5]; quadFactors0_1728_5
except:
    print("-1,1,0 invalid")
#Observation:  Everything is factorable.  Further,
#of six possible isogenies, all 6 are valid.
#
# 1) [0,1,-1] returns (x + a2) * (x - a2) * (x + (1/2*a0*a1 + 1/2)*a2) * (x + (-1/2*a0*a1 - 1/2)*a2) * (x + 1/2*a3) * (x - 1/2*a3)
# 2) [0,-1,1] returns (x - a1*a2) * (x + a1*a2) * (x + (-1/2*a1 + 1/2*a0)*a2) * (x + (1/2*a1 - 1/2*a0)*a2) * (x - 1/2*a1*a3) * (x + 1/2*a1*a3)
# 3) [1,0,-1] returns (x - a1*a2) * (x + a1*a2) * (x + (-1/4*a1 - 1/4*a0)*a3) * (x + (1/4*a1 + 1/4*a0)*a3) * (x + (-1/2*a1 - 1/2*a0)*a2) * (x + (1/2*a1 + 1/2*a0)*a2)
# 4) [1,-1,0] returns (x + (-1/4*a0*a1 - 1/4)*a3) * (x + (1/4*a0*a1 + 1/4)*a3) * (x + (-1/2*a0*a1 - 1/2)*a2) * (x + (1/2*a0*a1 + 1/2)*a2) * (x + (-1/2*a0*a1 + 1/2)*a2) * (x + (1/2*a0*a1 - 1/2)*a2)
# 5) [-1,0,1] returns (x - a2) * (x + a2) * (x + (-1/4*a0*a1 + 1/4)*a3) * (x + (1/4*a0*a1 - 1/4)*a3) * (x + (-1/2*a0*a1 + 1/2)*a2) * (x + (1/2*a0*a1 - 1/2)*a2)
# 6) [-1,1,0] returns (x + (-1/4*a1 + 1/4*a0)*a3) * (x + (1/4*a1 - 1/4*a0)*a3) * (x + (-1/2*a1 + 1/2*a0)*a2) * (x + (1/2*a1 - 1/2*a0)*a2) * (x + (-1/2*a1 - 1/2*a0)*a2) * (x + (1/2*a1 + 1/2*a0)*a2)
#
# We now run the remainder of our calculations to identify unique curves, etc.
print("")
isogenies0_1728_1 = isogenyEH([[1,r3,r3^2],[0,-1,1]],[funct0,funct1728],5);
isogenies0_1728_1[1]
polyRoots0_1728_1 = [x-n for n in isogenies0_1728_1[1]]
prodRoots0_1728_1 = 1
for i in range(len(polyRoots0_1728_1)):
    prodRoots0_1728_1 = prodRoots0_1728_1*polyRoots0_1728_1[i]
prodRoots0_1728_1
Igusa0_1728_1 = [symmetric2_6root(polyRoots0_1728_1),symmetric4_6root(polyRoots0_1728_1),symmetric6_6root(polyRoots0_1728_1),symmetric10_6root(polyRoots0_1728_1)]
Kohel0_1728_1 = IgusaToKohel(Igusa0_1728_1)
Kohel0_1728_1
print("")
isogenies0_1728_2 = isogenyEH([[1,r3,r3^2],[1,0,-1]],[funct0,funct1728],5);
isogenies0_1728_2[1]
polyRoots0_1728_2 = [x-n for n in isogenies0_1728_2[1]]
prodRoots0_1728_2 = 1
for i in range(len(polyRoots0_1728_2)):
    prodRoots0_1728_2 = prodRoots0_1728_2*polyRoots0_1728_2[i]
prodRoots0_1728_2
Igusa0_1728_2 = [symmetric2_6root(polyRoots0_1728_2),symmetric4_6root(polyRoots0_1728_2),symmetric6_6root(polyRoots0_1728_2),symmetric10_6root(polyRoots0_1728_2)]
Kohel0_1728_2 = IgusaToKohel(Igusa0_1728_2)
Kohel0_1728_2
print("")
isogenies0_1728_3 = isogenyEH([[1,r3,r3^2],[1,-1,0]],[funct0,funct1728],5);
isogenies0_1728_3[1]
polyRoots0_1728_3 = [x-n for n in isogenies0_1728_3[1]]
prodRoots0_1728_3 = 1
for i in range(len(polyRoots0_1728_3)):
    prodRoots0_1728_3 = prodRoots0_1728_3*polyRoots0_1728_3[i]
prodRoots0_1728_3
Igusa0_1728_3 = [symmetric2_6root(polyRoots0_1728_3),symmetric4_6root(polyRoots0_1728_3),symmetric6_6root(polyRoots0_1728_3),symmetric10_6root(polyRoots0_1728_3)]
Kohel0_1728_3 = IgusaToKohel(Igusa0_1728_3)
Kohel0_1728_3
print("")
isogenies0_1728_4 = isogenyEH([[1,r3,r3^2],[-1,0,1]],[funct0,funct1728],5);
isogenies0_1728_4[1]
polyRoots0_1728_4 = [x-n for n in isogenies0_1728_4[1]]
prodRoots0_1728_4 = 1
for i in range(len(polyRoots0_1728_4)):
    prodRoots0_1728_4 = prodRoots0_1728_4*polyRoots0_1728_4[i]
prodRoots0_1728_4
Igusa0_1728_4 = [symmetric2_6root(polyRoots0_1728_4),symmetric4_6root(polyRoots0_1728_4),symmetric6_6root(polyRoots0_1728_4),symmetric10_6root(polyRoots0_1728_4)]
Kohel0_1728_4 = IgusaToKohel(Igusa0_1728_4)
Kohel0_1728_4
print("")
isogenies0_1728_5 = isogenyEH([[1,r3,r3^2],[-1,1,0]],[funct0,funct1728],5);
isogenies0_1728_5[1]
polyRoots0_1728_5 = [x-n for n in isogenies0_1728_5[1]]
prodRoots0_1728_5 = 1
for i in range(len(polyRoots0_1728_5)):
    prodRoots0_1728_5 = prodRoots0_1728_5*polyRoots0_1728_5[i]
prodRoots0_1728_5
Igusa0_1728_5 = [symmetric2_6root(polyRoots0_1728_5),symmetric4_6root(polyRoots0_1728_5),symmetric6_6root(polyRoots0_1728_5),symmetric10_6root(polyRoots0_1728_5)]
Kohel0_1728_5 = IgusaToKohel(Igusa0_1728_5)
Kohel0_1728_5
#Observations:
#All six curves have the same invariants and hence are the same curve.  There is thus a degree 6 edge from (0,1728) to (20565/4, 4218240, 1235728).  Moving forward we will use x^6 - 3/2*a0*a1*x^4 - 9/2*x^2 + 3/2*a0*a1 to represent this curve.
#We lastly observe the automorphism type of this curve:
automorphismDetector([n for n in specifiedReverseRosenhain6LFT((-1/2*a1 - 1/2*a0)*a2, (-1/2*a1 + 1/2*a0)*a2, (1/2*a1 - 1/2*a0)*a2, (1/2*a1 + 1/2*a0)*a2, (-1/4*a1 + 1/4*a0)*a3, (1/4*a1 - 1/4*a0)*a3)])
#This curve is a type 1 curve as it has an automorphism group of "C2".  In addition, the Orbit-Stabilizer theorem tells us that the reverse edge must have a weight of 1.
︡a62065bb-2dc9-4845-a2c2-54f37f8b8d62︡
︠bd82694e-69b4-4392-b185-3914b6c66efcs︠
#Next, we need to calculate the elliptic neighbors of (0,1728)
ellip1728_0 = ellip1728
im1728_1 = ellip1728_0.isogeny(ellip1728_0(0,0))
ellip1728_1 = im1728_1.codomain(); ellip1728_1; ellip1728_1.j_invariant(); funct1728_1 = functFromEC(ellip1728_1); roots1728_1 = funct1728_1[1].roots(); roots1728_1
im1728_2 = ellip1728_0.isogeny(ellip1728_0(1,0)); ellip1728_2 = im1728_2.codomain(); ellip1728_2; ellip1728_2.j_invariant(); funct1728_2 = functFromEC(ellip1728_2); roots1728_2 = funct1728_2[1].roots(); roots1728_2
im1728_3 = ellip1728_0.isogeny(ellip1728_0(-1,0)); ellip1728_3 = im1728_3.codomain(); ellip1728_3; ellip1728_3.j_invariant(); funct1728_3 = functFromEC(ellip1728_3); roots1728_3 = funct1728_3[1].roots(); roots1728_3
print("")
ellip0_0 = ellip0
funct0_0 = functFromEC(ellip0_0)[1]
funct0_1 = veluIsogenyDegree2(funct0_0,funct0_0.roots()[0][0]); coeffs0_1 = funct0_1.coefficients(sparse=False); ellip0_1 = EllipticCurve([0,coeffs0_1[2],0,coeffs0_1[1],coeffs0_1[0]]); roots0_1 = funct0_1.roots(); roots0_1; ellip0_1; ellip0_1.j_invariant()
funct0_2 = veluIsogenyDegree2(funct0_0,funct0_0.roots()[1][0]); coeffs0_2 = funct0_2.coefficients(sparse=False); ellip0_2 = EllipticCurve([0,coeffs0_2[2],0,coeffs0_2[1],coeffs0_2[0]]); roots0_2 = funct0_2.roots(); roots0_2; ellip0_2; ellip0_2.j_invariant()
funct0_3 = veluIsogenyDegree2(funct0_0,funct0_0.roots()[2][0]); coeffs0_3 = funct0_3.coefficients(sparse=False); ellip0_3 = EllipticCurve([0,coeffs0_3[2],0,coeffs0_3[1],coeffs0_3[0]]); roots0_3 = funct0_3.roots(); roots0_3; ellip0_3; ellip0_3.j_invariant()
#Observations:
#
#For 1728, two of these curves are the same curve: j = 287496.  We will represent this curve with the equation y^2 = x^3 - 11*x - 14.  The other is a self loop: j = 1728.
#For 0, all of the curves are the same curve: j = 54000.  We will represent this curve with the equation y^2 = x^3 - 15*x - 22.
#Combinatorially, we have 2 outgoing edges to elliptic products:
#We have a weight 3 edge to a curve of type Pi 1728 (E12) with invariants: (1728, 54000).  Since a type Pi 1728 curve has autormorphism group of order 4, we can conclude the dual edge has weight 1.
#We have a weight 6 edge to a curve of type Pi (E23) with invariants (54000, 287496).  Since a type Pi curve has automorphism group of order 2, we can conclude the dual edge has weight 1.
#
#At this point we have drawn the entirety of the radius 1 neighborhood.  We have 3 neighboring curves to consider.
︡8eb4c4c4-0553-4917-a374-d14606faec16︡
︠cd72ba4f-bcf6-4527-be1e-377d605f55a1s︠
#Investigating the type 1 curve.
#
#We begin by calculating the Hyperelliptic Neighbors of the type 1 curve.  By Florit and Smith's Atlas, we expect 14 of these.  The only exception will be the dual edge back to (0, 1728) and hence we don't need to bother calculating it.

#Calculating Richelot Isogenies relies on a Lie Bracket, we define that here so we can call it quickly whenever.  Note that the "d" function takes derivatives and the "Lie" function calculates the Lie Bracket for a quadratic pairing.
def d(poly):
    return(diff(poly, x))

def Lie(A,B):
    return([B[1]*A[2]-B[2]*A[1],B[2]*A[0]-B[0]*A[2],B[0]*A[1]-B[1]*A[0]])
︡4e104df2-58e9-47b0-9240-cd253c24b8ea︡
︠b64059c5-22c8-4056-b9f3-4bc788760eb9s︠
#Here are our observations thus far:
#0) Fully factor, factors present are:
#     (x + ((-1/2*a0 + 1/2)*a1)*a2), (x + ((1/2*a0 + 1/2)*a1)*a2), (x + ((-1/2*a0 - 1/2)*a1)*a2), (x + ((1/2*a0 - 1/2)*a1)*a2), (x - a1*a2), (x + a1*a2)
#1) Does not fully factor.  Factors present are:
#     (x^2 + ((-3/2*a1 - 1/2*a0)*a3 - 2*a1*a2)*x + ((-1/4*a0*a1 + 3/4)*a2)*a3 + 1/2*a0*a1), (x^2 + ((16/43*a1 - 12/43*a0)*a3 + (-29/43*a1 + 11/43*a0)*a2)*x + ((6/43*a0*a1 + 8/43)*a2)*a3 + 57/86*a0*a1 + 33/86), (x^2 + ((-10/43*a1 + 14/43*a0)*a3 + (-29/43*a1 - 11/43*a0)*a2)*x + ((-7/43*a0*a1 - 5/43)*a2)*a3 + 57/86*a0*a1 - 33/86)
#2) Does not fully factor.  Factors present are:
#     (x^2 + ((3/2*a1 + 1/2*a0)*a3 - 2*a1*a2)*x + ((1/4*a0*a1 - 3/4)*a2)*a3 + 1/2*a0*a1), (x^2 + ((-16/43*a1 + 12/43*a0)*a3 + (-29/43*a1 + 11/43*a0)*a2)*x + ((-6/43*a0*a1 - 8/43)*a2)*a3 + 57/86*a0*a1 + 33/86), (x^2 + ((10/43*a1 - 14/43*a0)*a3 + (-29/43*a1 - 11/43*a0)*a2)*x + ((7/43*a0*a1 + 5/43)*a2)*a3 + 57/86*a0*a1 - 33/86)
#3) Fully factors.  Factors present are:
#     (x + (1/2*a0 - 1/2)*a2), (x + (1/2*a0 + 1/2)*a2), (x + (-1/2*a0 - 1/2)*a2), (x + (-1/2*a0 + 1/2)*a2), (x - a2), (x + a2)
#4) Does not fully factor.  Factors present are:
#     (x^2 + ((1/2*a1 + 1/6*a0)*a3 + 2/3*a0*a2)*x + ((1/4*a0*a1 + 1/4)*a2)*a3 + 1/2*a0*a1), (x^2 + ((4/19*a1 - 16/57*a0)*a3 + (3/19*a1 - 31/57*a0)*a2)*x + ((-2/19*a0*a1 + 8/19)*a2)*a3 - 7/38*a0*a1 + 9/38), (x^2 + ((-6/19*a1 + 14/57*a0)*a3 + (-3/19*a1 - 31/57*a0)*a2)*x + ((3/19*a0*a1 - 7/19)*a2)*a3 - 7/38*a0*a1 - 9/38)
#5) Does not fully factor.  Factors present are:
#     (x^2 + ((-1/2*a1 - 1/6*a0)*a3 + 2/3*a0*a2)*x + ((-1/4*a0*a1 - 1/4)*a2)*a3 + 1/2*a0*a1), (x^2 + ((-4/19*a1 + 16/57*a0)*a3 + (3/19*a1 - 31/57*a0)*a2)*x + ((2/19*a0*a1 - 8/19)*a2)*a3 - 7/38*a0*a1 + 9/38), (x^2 + ((6/19*a1 - 14/57*a0)*a3 + (-3/19*a1 - 31/57*a0)*a2)*x + ((-3/19*a0*a1 + 7/19)*a2)*a3 - 7/38*a0*a1 - 9/38)
#6) Is elliptic and hence the dual edge back to (0, 1728)
#7) Does not fully factor.  Factors present are:
#     (x^2 + ((1/4*a0*a1 - 1/4)*a2)*a3), (x + ((-1/2*a0 + 1/2)*a1 - 1/2*a0 + 1/2)*a3 + (-1/2*a1 + 1/2*a0 - 1)*a2), (x + ((1/2*a0 + 1/2)*a1 - 1/2*a0 - 1/2)*a3 + (-1/2*a1 + 1/2*a0 + 1)*a2), (x + ((-1/2*a0 - 1/2)*a1 + 1/2*a0 + 1/2)*a3 + (1/2*a1 - 1/2*a0 - 1)*a2), (x + ((1/2*a0 - 1/2)*a1 + 1/2*a0 - 1/2)*a3 + (1/2*a1 - 1/2*a0 + 1)*a2)
#8) Does not fully factor.  Factors present are:
#     (x^2 + ((-1/4*a0*a1 + 1/4)*a2)*a3), (x + ((-1/2*a0 - 1/2)*a1 + 1/2*a0 + 1/2)*a3 + (-1/2*a1 + 1/2*a0 + 1)*a2), (x + ((1/2*a0 - 1/2)*a1 + 1/2*a0 - 1/2)*a3 + (-1/2*a1 + 1/2*a0 - 1)*a2), (x + ((-1/2*a0 + 1/2)*a1 - 1/2*a0 + 1/2)*a3 + (1/2*a1 - 1/2*a0 + 1)*a2), (x + ((1/2*a0 + 1/2)*a1 - 1/2*a0 - 1/2)*a3 + (1/2*a1 - 1/2*a0 - 1)*a2)
#9) Does not fully factor.  Factors present are:
#     (x + (1/2*a1 - 1/2*a0 - 1)*a3 + (-1/2*a1 - 1/2*a0 - 1)*a2), (x + (1/2*a1 - 1/2*a0 + 1)*a3 + (-1/2*a1 - 1/2*a0 + 1)*a2), (x^2 + 1/2*a2*a3), (x + (-1/2*a1 + 1/2*a0 - 1)*a3 + (1/2*a1 + 1/2*a0 - 1)*a2), (x + (-1/2*a1 + 1/2*a0 + 1)*a3 + (1/2*a1 + 1/2*a0 + 1)*a2)
#10) Does not fully factor.  Factors present are:
#     (x^2 + ((6/19*a1 - 14/57*a0)*a3 + (3/19*a1 + 31/57*a0)*a2)*x + ((3/19*a0*a1 - 7/19)*a2)*a3 - 7/38*a0*a1 - 9/38), (x^2 + ((-1/2*a1 - 1/6*a0)*a3 - 2/3*a0*a2)*x + ((1/4*a0*a1 + 1/4)*a2)*a3 + 1/2*a0*a1), (x^2 + ((-4/19*a1 + 16/57*a0)*a3 + (-3/19*a1 + 31/57*a0)*a2)*x + ((-2/19*a0*a1 + 8/19)*a2)*a3 - 7/38*a0*a1 + 9/38)
#11) Does not fully factor.  Factors present are:
#     (x^2 + ((10/43*a1 - 14/43*a0)*a3 + (29/43*a1 + 11/43*a0)*a2)*x + ((-7/43*a0*a1 - 5/43)*a2)*a3 + 57/86*a0*a1 - 33/86), (x^2 + ((-16/43*a1 + 12/43*a0)*a3 + (29/43*a1 - 11/43*a0)*a2)*x + ((6/43*a0*a1 + 8/43)*a2)*a3 + 57/86*a0*a1 + 33/86), (x^2 + ((3/2*a1 + 1/2*a0)*a3 + 2*a1*a2)*x + ((-1/4*a0*a1 + 3/4)*a2)*a3 + 1/2*a0*a1)
#12) Does not fully factor.  Factors present are:
#     (x + (-1/2*a1 + 1/2*a0 - 1)*a3 + (-1/2*a1 - 1/2*a0 + 1)*a2), (x + (-1/2*a1 + 1/2*a0 + 1)*a3 + (-1/2*a1 - 1/2*a0 - 1)*a2), (x^2 - 1/2*a2*a3), (x + (1/2*a1 - 1/2*a0 - 1)*a3 + (1/2*a1 + 1/2*a0 + 1)*a2), (x + (1/2*a1 - 1/2*a0 + 1)*a3 + (1/2*a1 + 1/2*a0 - 1)*a2)
#13) Does not fully factor.  Factors present are:
#     (x^2 + ((-6/19*a1 + 14/57*a0)*a3 + (3/19*a1 + 31/57*a0)*a2)*x + ((-3/19*a0*a1 + 7/19)*a2)*a3 - 7/38*a0*a1 - 9/38), (x^2 + ((1/2*a1 + 1/6*a0)*a3 - 2/3*a0*a2)*x + ((-1/4*a0*a1 - 1/4)*a2)*a3 + 1/2*a0*a1), (x^2 + ((4/19*a1 - 16/57*a0)*a3 + (-3/19*a1 + 31/57*a0)*a2)*x + ((2/19*a0*a1 - 8/19)*a2)*a3 - 7/38*a0*a1 + 9/38)
#14) Does not fully factor.  Factors present are:
#     (x^2 + ((-10/43*a1 + 14/43*a0)*a3 + (29/43*a1 + 11/43*a0)*a2)*x + ((7/43*a0*a1 + 5/43)*a2)*a3 + 57/86*a0*a1 - 33/86), (x^2 + ((16/43*a1 - 12/43*a0)*a3 + (29/43*a1 - 11/43*a0)*a2)*x + ((-6/43*a0*a1 - 8/43)*a2)*a3 + 57/86*a0*a1 + 33/86), (x^2 + ((-3/2*a1 - 1/2*a0)*a3 + 2*a1*a2)*x + ((1/4*a0*a1 - 3/4)*a2)*a3 + 1/2*a0*a1)
︡5ff75b7b-82aa-4c9d-b132-12108f90eb50︡
︠34c02f58-3d52-4f11-8207-4093740be525s︠
#Next, we need to calculate the elliptic neighbors of 54000
roots54000_0 = [n[0] for n in roots0_1]; roots54000_0
ellip54000_0 = EllipticCurve(interFieldLevel3,[0,0,0,-15,-22])
im54000_1 = ellip54000_0.isogeny(x-roots54000_0[0]); ellip54000_1 = im54000_1.codomain(); ellip54000_1; ellip54000_1.j_invariant(); funct54000_1 = functFromEC(ellip54000_1); roots54000_1 = funct54000_1[1].roots(); roots54000_1
im54000_2 = ellip54000_0.isogeny(x-roots54000_0[1]);
ellip54000_2 = im54000_2.codomain(); ellip54000_2; ellip54000_2.j_invariant(); funct54000_2 = functFromEC(ellip54000_2); roots54000_2 = funct54000_2[1].roots(); roots54000_2
im54000_3 = ellip54000_0.isogeny(x-roots54000_0[2]); ellip54000_3 = im54000_3.codomain(); ellip54000_3; ellip54000_3.j_invariant(); funct54000_3 = functFromEC(ellip54000_3); roots54000_3 = funct54000_3[1].roots(); roots54000_3
︡30f32777-6c81-4f52-ba43-ec5b1afb6271︡
︠b7a906de-f87f-4bd9-8402-5bbb9ee213cfs︠
funct54000 = x^3 - 15*x -22
try:
    isogenies1728_54000_0 = isogenyEH([[0,1,-1],[-2,-2*a0+1,2*a0+1]],[funct1728,funct54000],5,True); quadFactors1728_54000_0 = [factor(n) for n in isogenies1728_54000_0]; quadFactors1728_54000_0
except:
    print("-2,-2*a0+1,2*a0+1 invalid")
try:
    isogenies1728_54000_1 = isogenyEH([[0,1,-1],[-2,2*a0+1,-2*a0+1]],[funct1728,funct54000],5,True); quadFactors1728_54000_1 = [factor(n) for n in isogenies1728_54000_1]; quadFactors1728_54000_1
except:
    print("-2,2*a0+1,-2*a0+1 invalid")
try:
    isogenies1728_54000_2 = isogenyEH([[0,1,-1],[2*a0+1,-2,-2*a0+1]],[funct1728,funct54000],5,True); quadFactors1728_54000_2 = [factor(n) for n in isogenies1728_54000_2]; quadFactors1728_54000_2
except:
    print("2*a0+1,-2,-2*a0+1, invalid")
try:
    isogenies1728_54000_3 = isogenyEH([[0,1,-1],[2*a0+1,-2*a0+1,-2]],[funct1728,funct54000],5,True); quadFactors1728_54000_3 = [factor(n) for n in isogenies1728_54000_3]; quadFactors1728_54000_3
except:
    print("2*a0+1,-2*a0+1,-2 invalid")
try:
    isogenies1728_54000_4 = isogenyEH([[0,1,-1],[-2*a0+1,-2,2*a0+1]],[funct1728,funct54000],5,True); quadFactors1728_54000_4 = [factor(n) for n in isogenies1728_54000_4]; quadFactors1728_54000_4
except:
    print("-2*a0+1,-2,2*a0+1 invalid")
try:
    isogenies1728_54000_5 = isogenyEH([[0,1,-1],[-2*a0+1,2*a0+1,-2]],[funct1728,funct54000],5,True); quadFactors1728_54000_5 = [factor(n) for n in isogenies1728_54000_5]; quadFactors1728_54000_5
except:
    print("-2*a0+1,2*a0+1,-2 invalid")
︡89bdd197-da45-46fe-a40d-2c1878b66a52︡
︠74593b6e-17e0-425b-8935-008efb9f5a5ds︠
# We now run the remainder of our calculations to identify unique curves, etc.
isogenies1728_54000_0 = isogenyEH([[0,1,-1],[-2,-2*a0+1,2*a0+1]],[funct1728,funct54000],5);
isogenies1728_54000_0[1]
polyRoots1728_54000_0 = [x-n for n in isogenies1728_54000_0[1]]
prodRoots1728_54000_0 = 1
for i in range(len(polyRoots1728_54000_0)):
    prodRoots1728_54000_0 = prodRoots1728_54000_0*polyRoots1728_54000_0[i]
prodRoots1728_54000_0
Igusa1728_54000_0 = [symmetric2_6root(polyRoots1728_54000_0),symmetric4_6root(polyRoots1728_54000_0),symmetric6_6root(polyRoots1728_54000_0),symmetric10_6root(polyRoots1728_54000_0)]
Kohel1728_54000_0 = IgusaToKohel(Igusa1728_54000_0)
Kohel1728_54000_0
print("")
isogenies1728_54000_1 = isogenyEH([[0,1,-1],[-2,2*a0+1,-2*a0+1]],[funct1728,funct54000],5);
isogenies1728_54000_1[1]
polyRoots1728_54000_1 = [x-n for n in isogenies1728_54000_1[1]]
prodRoots1728_54000_1 = 1
for i in range(len(polyRoots1728_54000_1)):
    prodRoots1728_54000_1 = prodRoots1728_54000_1*polyRoots1728_54000_1[i]
prodRoots1728_54000_1
Igusa1728_54000_1 = [symmetric2_6root(polyRoots1728_54000_1),symmetric4_6root(polyRoots1728_54000_1),symmetric6_6root(polyRoots1728_54000_1),symmetric10_6root(polyRoots1728_54000_1)]
Kohel1728_54000_1 = IgusaToKohel(Igusa1728_54000_1)
Kohel1728_54000_1
print("")
isogenies1728_54000_2 = isogenyEH([[0,1,-1],[-2*a0+1,-2,2*a0+1]],[funct1728,funct54000],5);
isogenies1728_54000_2[1]
polyRoots1728_54000_2 = [x-n for n in isogenies1728_54000_2[1]]
prodRoots1728_54000_2 = 1
for i in range(len(polyRoots1728_54000_2)):
    prodRoots1728_54000_2 = prodRoots1728_54000_2*polyRoots1728_54000_2[i]
prodRoots1728_54000_2
Igusa1728_54000_2 = [symmetric2_6root(polyRoots1728_54000_2),symmetric4_6root(polyRoots1728_54000_2),symmetric6_6root(polyRoots1728_54000_2),symmetric10_6root(polyRoots1728_54000_2)]
Kohel1728_54000_2 = IgusaToKohel(Igusa1728_54000_2)
Kohel1728_54000_2
print("")
isogenies1728_54000_3 = isogenyEH([[0,1,-1],[-2*a0+1,2*a0+1,-2]],[funct1728,funct54000],5);
isogenies1728_54000_3[1]
polyRoots1728_54000_3 = [x-n for n in isogenies1728_54000_3[1]]
prodRoots1728_54000_3 = 1
for i in range(len(polyRoots1728_54000_3)):
    prodRoots1728_54000_3 = prodRoots1728_54000_3*polyRoots1728_54000_3[i]
prodRoots1728_54000_3
Igusa1728_54000_3 = [symmetric2_6root(polyRoots1728_54000_3),symmetric4_6root(polyRoots1728_54000_3),symmetric6_6root(polyRoots1728_54000_3),symmetric10_6root(polyRoots1728_54000_3)]
Kohel1728_54000_3 = IgusaToKohel(Igusa1728_54000_3)
Kohel1728_54000_3
print("")
isogenies1728_54000_4 = isogenyEH([[0,1,-1],[2*a0+1,-2,-2*a0+1]],[funct1728,funct54000],5);
isogenies1728_54000_4[1]
polyRoots1728_54000_4 = [x-n for n in isogenies1728_54000_4[1]]
prodRoots1728_54000_4 = 1
for i in range(len(polyRoots1728_54000_4)):
    prodRoots1728_54000_4 = prodRoots1728_54000_4*polyRoots1728_54000_4[i]
prodRoots1728_54000_4
Igusa1728_54000_4 = [symmetric2_6root(polyRoots1728_54000_4),symmetric4_6root(polyRoots1728_54000_4),symmetric6_6root(polyRoots1728_54000_4),symmetric10_6root(polyRoots1728_54000_4)]
Kohel1728_54000_4 = IgusaToKohel(Igusa1728_54000_4)
Kohel1728_54000_4
print("")
isogenies1728_54000_5 = isogenyEH([[0,1,-1],[2*a0+1,-2*a0+1,-2]],[funct1728,funct54000],5);
isogenies1728_54000_5[1]
polyRoots1728_54000_5 = [x-n for n in isogenies1728_54000_5[1]]
prodRoots1728_54000_5 = 1
for i in range(len(polyRoots1728_54000_5)):
    prodRoots1728_54000_5 = prodRoots1728_54000_5*polyRoots1728_54000_5[i]
prodRoots1728_54000_5
Igusa1728_54000_5 = [symmetric2_6root(polyRoots1728_54000_5),symmetric4_6root(polyRoots1728_54000_5),symmetric6_6root(polyRoots1728_54000_5),symmetric10_6root(polyRoots1728_54000_5)]
Kohel1728_54000_5 = IgusaToKohel(Igusa1728_54000_5)
Kohel1728_54000_5
︡4a672800-fdde-4404-b951-c417dfeac452︡
︠fa1b616b-eef4-4884-b5cd-839b356a5fbas︠
uniqueSetCounter([Kohel1728_54000_0,Kohel1728_54000_1,Kohel1728_54000_2,Kohel1728_54000_3,Kohel1728_54000_4,Kohel1728_54000_5])
︡df1b66d8-1550-4043-aae9-430fd0d98d40︡
︠d8cd1003-43db-4505-8c06-3bf8df982079s︠
funct287496 = x^3 - 11*x - 14
funct287496.roots()
try:
    isogenies54000_287496_0 = isogenyEH([[-2,-2*a0+1,2*a0+1],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5,True); quadFactors54000_287496_0 = [factor(n) for n in isogenies54000_287496_0]; quadFactors54000_287496_0
except:
    print("-2,-2*a0+1,2*a0+1 invalid")
try:
    isogenies54000_287496_1 = isogenyEH([[-2,2*a0+1,-2*a0+1],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5,True); quadFactors54000_287496_1 = [factor(n) for n in isogenies54000_287496_1]; quadFactors54000_287496_1
except:
    print("-2,2*a0+1,-2*a0+1 invalid")
try:
    isogenies54000_287496_2 = isogenyEH([[2*a0+1,-2,-2*a0+1],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5,True); quadFactors54000_287496_2 = [factor(n) for n in isogenies54000_287496_2]; quadFactors54000_287496_2
except:
    print("2*a0+1,-2,-2*a0+1, invalid")
try:
    isogenies54000_287496_3 = isogenyEH([[2*a0+1,-2*a0+1,-2],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5,True); quadFactors54000_287496_3 = [factor(n) for n in isogenies54000_287496_3]; quadFactors54000_287496_3
except:
    print("2*a0+1,-2*a0+1,-2 invalid")
try:
    isogenies54000_287496_4 = isogenyEH([[-2*a0+1,-2,2*a0+1],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5,True); quadFactors54000_287496_4 = [factor(n) for n in isogenies54000_287496_4]; quadFactors54000_287496_4
except:
    print("-2*a0+1,-2,2*a0+1 invalid")
try:
    isogenies54000_287496_5 = isogenyEH([[-2*a0+1,2*a0+1,-2],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5,True); quadFactors54000_287496_5 = [factor(n) for n in isogenies54000_287496_5]; quadFactors54000_287496_5
except:
    print("-2*a0+1,2*a0+1,-2 invalid")
︡526cb0ca-c1aa-47e7-880e-f44ad747472a︡
︠55993d52-589c-4a82-a1d3-38c276728375s︠
setupFieldLevel4.<y4> = interFieldLevel3[]
interFieldLevel4.<a4> = setupFieldLevel4.quotient(y4^2 + ((1/4*a0*a1 + 1/4)*a2)*a3)
polyFieldLevel3.<x> = interFieldLevel4[]
F = interFieldLevel4
︡5e3d33df-99a6-4386-9812-05e6c2e0b13c︡
︠c7de7450-3b25-4ea5-acb7-0727d08dc267s︠
funct287496 = x^3 - 11*x - 14
funct287496.roots()
try:
    isogenies54000_287496_0 = isogenyEH([[-2,-2*a0+1,2*a0+1],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5,True); quadFactors54000_287496_0 = [factor(n) for n in isogenies54000_287496_0]; quadFactors54000_287496_0
except:
    print("-2,-2*a0+1,2*a0+1 invalid")
try:
    isogenies54000_287496_1 = isogenyEH([[-2,2*a0+1,-2*a0+1],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5,True); quadFactors54000_287496_1 = [factor(n) for n in isogenies54000_287496_1]; quadFactors54000_287496_1
except:
    print("-2,2*a0+1,-2*a0+1 invalid")
try:
    isogenies54000_287496_2 = isogenyEH([[2*a0+1,-2,-2*a0+1],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5,True); quadFactors54000_287496_2 = [factor(n) for n in isogenies54000_287496_2]; quadFactors54000_287496_2
except:
    print("2*a0+1,-2,-2*a0+1, invalid")
try:
    isogenies54000_287496_3 = isogenyEH([[2*a0+1,-2*a0+1,-2],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5,True); quadFactors54000_287496_3 = [factor(n) for n in isogenies54000_287496_3]; quadFactors54000_287496_3
except:
    print("2*a0+1,-2*a0+1,-2 invalid")
try:
    isogenies54000_287496_4 = isogenyEH([[-2*a0+1,-2,2*a0+1],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5,True); quadFactors54000_287496_4 = [factor(n) for n in isogenies54000_287496_4]; quadFactors54000_287496_4
except:
    print("-2*a0+1,-2,2*a0+1 invalid")
try:
    isogenies54000_287496_5 = isogenyEH([[-2*a0+1,2*a0+1,-2],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5,True); quadFactors54000_287496_5 = [factor(n) for n in isogenies54000_287496_5]; quadFactors54000_287496_5
except:
    print("-2*a0+1,2*a0+1,-2 invalid")
︡68437375-92ff-4735-9c54-2a8807aaae6b︡
︠a42342e9-baa6-47b5-b1bc-2013a4603390s︠
# We now run the remainder of our calculations to identify unique curves, etc.
isogenies54000_287496_0 = isogenyEH([[-2,-2*a0+1,2*a0+1],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5);
isogenies54000_287496_0[1]
polyRoots54000_287496_0 = [x-n for n in isogenies54000_287496_0[1]]
prodRoots54000_287496_0 = 1
for i in range(len(polyRoots54000_287496_0)):
    prodRoots54000_287496_0 = prodRoots54000_287496_0*polyRoots54000_287496_0[i]
prodRoots54000_287496_0
Igusa54000_287496_0 = [symmetric2_6root(polyRoots54000_287496_0),symmetric4_6root(polyRoots54000_287496_0),symmetric6_6root(polyRoots54000_287496_0),symmetric10_6root(polyRoots54000_287496_0)]
Kohel54000_287496_0 = IgusaToKohel(Igusa54000_287496_0)
Kohel54000_287496_0
print("")
isogenies54000_287496_1 = isogenyEH([[-2,2*a0+1,-2*a0+1],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5);
isogenies54000_287496_1[1]
polyRoots54000_287496_1 = [x-n for n in isogenies54000_287496_1[1]]
prodRoots54000_287496_1 = 1
for i in range(len(polyRoots54000_287496_1)):
    prodRoots54000_287496_1 = prodRoots54000_287496_1*polyRoots54000_287496_1[i]
prodRoots54000_287496_1
Igusa54000_287496_1 = [symmetric2_6root(polyRoots54000_287496_1),symmetric4_6root(polyRoots54000_287496_1),symmetric6_6root(polyRoots54000_287496_1),symmetric10_6root(polyRoots54000_287496_1)]
Kohel54000_287496_1 = IgusaToKohel(Igusa54000_287496_1)
Kohel54000_287496_1
print("")
isogenies54000_287496_2 = isogenyEH([[-2*a0+1,-2,2*a0+1],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5);
isogenies54000_287496_2[1]
polyRoots54000_287496_2 = [x-n for n in isogenies54000_287496_2[1]]
prodRoots54000_287496_2 = 1
for i in range(len(polyRoots54000_287496_2)):
    prodRoots54000_287496_2 = prodRoots54000_287496_2*polyRoots54000_287496_2[i]
prodRoots54000_287496_2
Igusa54000_287496_2 = [symmetric2_6root(polyRoots54000_287496_2),symmetric4_6root(polyRoots54000_287496_2),symmetric6_6root(polyRoots54000_287496_2),symmetric10_6root(polyRoots54000_287496_2)]
Kohel54000_287496_2 = IgusaToKohel(Igusa54000_287496_2)
Kohel54000_287496_2
print("")
isogenies54000_287496_3 = isogenyEH([[-2*a0+1,2*a0+1,-2],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5);
isogenies54000_287496_3[1]
polyRoots54000_287496_3 = [x-n for n in isogenies54000_287496_3[1]]
prodRoots54000_287496_3 = 1
for i in range(len(polyRoots54000_287496_3)):
    prodRoots54000_287496_3 = prodRoots54000_287496_3*polyRoots54000_287496_3[i]
prodRoots54000_287496_3
Igusa54000_287496_3 = [symmetric2_6root(polyRoots54000_287496_3),symmetric4_6root(polyRoots54000_287496_3),symmetric6_6root(polyRoots54000_287496_3),symmetric10_6root(polyRoots54000_287496_3)]
Kohel54000_287496_3 = IgusaToKohel(Igusa54000_287496_3)
Kohel54000_287496_3
print("")
isogenies54000_287496_4 = isogenyEH([[2*a0+1,-2,-2*a0+1],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5);
isogenies54000_287496_4[1]
polyRoots54000_287496_4 = [x-n for n in isogenies54000_287496_4[1]]
prodRoots54000_287496_4 = 1
for i in range(len(polyRoots54000_287496_4)):
    prodRoots54000_287496_4 = prodRoots54000_287496_4*polyRoots54000_287496_4[i]
prodRoots54000_287496_4
Igusa54000_287496_4 = [symmetric2_6root(polyRoots54000_287496_4),symmetric4_6root(polyRoots54000_287496_4),symmetric6_6root(polyRoots54000_287496_4),symmetric10_6root(polyRoots54000_287496_4)]
Kohel54000_287496_4 = IgusaToKohel(Igusa54000_287496_4)
Kohel54000_287496_4
print("")
isogenies54000_287496_5 = isogenyEH([[2*a0+1,-2*a0+1,-2],[-2,((a1 + 1/3*a0)*a2)*a3 + 1,((-a1 - 1/3*a0)*a2)*a3 + 1]],[funct54000,funct287496],5);
isogenies54000_287496_5[1]
polyRoots54000_287496_5 = [x-n for n in isogenies54000_287496_5[1]]
prodRoots54000_287496_5 = 1
for i in range(len(polyRoots54000_287496_5)):
    prodRoots54000_287496_5 = prodRoots54000_287496_5*polyRoots54000_287496_5[i]
prodRoots54000_287496_5
Igusa54000_287496_5 = [symmetric2_6root(polyRoots54000_287496_5),symmetric4_6root(polyRoots54000_287496_5),symmetric6_6root(polyRoots54000_287496_5),symmetric10_6root(polyRoots54000_287496_5)]
Kohel54000_287496_5 = IgusaToKohel(Igusa54000_287496_5)
Kohel54000_287496_5
︡7d0b8e19-ede3-4ac0-9af0-8cd3bf14e9a1︡
︠6a4f9912-eb6c-4cb2-8cff-7af1dbc25906s︠
#What neighbors appear in our graph?
uniqueSetCounter([Kohel54000_287496_0,Kohel54000_287496_1,Kohel54000_287496_2,Kohel54000_287496_3,Kohel54000_287496_4,Kohel54000_287496_5])
︡37ec0b48-cede-4c5d-a7a0-c8a21aa61283︡
︠ce750ecf-9cc6-4271-885e-66bd142141a2s︠
#We now prepare to calculate the isogenies of the Type 1 curve.
rootSet = (x^6 - 3/2*a0*a1*x^4 - 9/2*x^2 + 3/2*a0*a1).roots()
rootsT1 = [x-n[0] for n in rootSet]
kernelBoxT1 = [[rootsT1[n[0]]*rootsT1[n[1]],rootsT1[n[2]]*rootsT1[n[3]],rootsT1[n[4]]*rootsT1[n[5]]] for n in Index]
determinantsT1 = [determinantDelta(n) for n in kernelBoxT1]
#We take the derivative
d_kernelBoxT1 = [[d(n) for n in m] for m in kernelBoxT1]
#And then we take the image of the Richelot Isogeny, where possible.
im_kernelBoxT1 = []
factoredBoxT1 = []
for i in range(len(kernelBoxT1)):
    if determinantsT1[i] == 0:
        im_kernelBoxT1.append(["Elliptic"])
        factoredBoxT1.append(["Elliptic"])
    else:
        im_kernelBoxT1.append(Lie(kernelBoxT1[i],d_kernelBoxT1[i]))
        factoredBoxT1.append([factor(n) for n in im_kernelBoxT1[i]])
︡0238dab0-8a25-44a9-98ce-6f2c750eaf9c︡
︠7c034325-2b89-4dbf-823c-acd97ae00369s︠
#Print the above results.
for i in range(len(factoredBoxT1)):
    print(i, factoredBoxT1[i])
    print("")
︡51fff2f4-4cc3-485a-8f63-7ba8403a896d︡
︠80948aa2-c21a-462f-82ff-f349dd294fdas︠
#We add once of the missing polynomials as a field extension.  It will turn out we only need one.
setupFieldLevel5.<y5> = interFieldLevel4[]
interFieldLevel5.<a5> = setupFieldLevel5.quotient(y5^2 + ((-3/2*a1 - 1/2*a0)*a3 - 2*a1*a2)*y5 + ((-1/4*a0*a1 + 3/4)*a2)*a3 + 1/2*a0*a1)
polyFieldLevel4.<x> = interFieldLevel5[]
F = interFieldLevel5
︡8d4250f2-c79f-4d30-9cf3-00871c4ec55e︡
︠1faa8b37-dc9c-42de-b2a5-7f67bbc36d2fs︠
#And we retest our factorizations.
rootsT1new = [x-n[0] for n in rootSet]
kernelBoxT1new = [[rootsT1new[n[0]]*rootsT1new[n[1]],rootsT1new[n[2]]*rootsT1new[n[3]],rootsT1new[n[4]]*rootsT1new[n[5]]] for n in Index]
#We take the derivative
d_kernelBoxT1new = [[d(n) for n in m] for m in kernelBoxT1new]
im_kernelBoxT1new = []
factoredBoxT1new = []
for i in range(len(kernelBoxT1new)):
    if determinantsT1[i] == 0:
        im_kernelBoxT1new.append(["Elliptic"])
        factoredBoxT1new.append(["Elliptic"])
    else:
        im_kernelBoxT1new.append(Lie(kernelBoxT1new[i],d_kernelBoxT1new[i]))
        factoredBoxT1new.append([factor(n) for n in im_kernelBoxT1new[i]])
factoredBoxT1new
︡cbb36513-6427-4bc7-aab4-3e878a10f014︡
︠a94e06bb-727d-4aae-8b42-04fedd3d136ds︠
#Now we find the Igusa invariants of each curve.
invarsT1Box = []
for i in range(len(kernelBoxT1new)):
    if determinantsT1[i] == 0:
        invarsT1Box.append(["Elliptic"])
    else:
        invarsT1Box.append(HyperellipticCurve(expand(factoredBoxT1new[i][0]*factoredBoxT1new[i][1]*factoredBoxT1new[i][2])).absolute_igusa_invariants_kohel())
uniqueIgusasT1 = uniqueSetCounter(invarsT1Box)
uniqueIgusasT1
︡6ce5589c-3c74-4b29-a3e2-21aa9bcacfc0︡
︠f59c1d67-6093-436b-833b-54c935cc8a20s︠
#Here we take all the neighbors of the radius 1 nodes and prep them to see if there are any common nodes by removing duplicates.
boxOfIgusas = [n[0] for n in uniqueIgusasT1]
input1728_54000 = uniqueSetCounter([Kohel1728_54000_0,Kohel1728_54000_1,Kohel1728_54000_2,Kohel1728_54000_3,Kohel1728_54000_4,Kohel1728_54000_5])
input54000_287496 = uniqueSetCounter([Kohel54000_287496_0,Kohel54000_287496_1,Kohel54000_287496_2,Kohel54000_287496_3,Kohel54000_287496_4,Kohel54000_287496_5])
boxOfIgusas_1728_54000 = [n[0] for n in input1728_54000]
boxOfIgusas_54000_287496 = [n[0] for n in input54000_287496]
boxOfIgusasRe = []
#The following three for-loops ensure everything is in the same field so that equality can be properly tested without types getting in the way.
for i in range(len(boxOfIgusas)):
    boxOfIgusasRe.append([])
    for j in range(3):
        try:
            boxOfIgusasRe[i].append(F(boxOfIgusas[i][j]))
        except:
            pass
boxOfIgusas_1728_54000Re = []
for i in range(len(boxOfIgusas_1728_54000)):
    boxOfIgusas_1728_54000Re.append([])
    for j in range(3):
        try:
            boxOfIgusas_1728_54000Re[i].append(F(boxOfIgusas_1728_54000[i][j]))
        except:
            pass
boxOfIgusas_54000_287496Re = []
for i in range(len(boxOfIgusas_54000_287496)):
    boxOfIgusas_54000_287496Re.append([])
    for j in range(3):
        try:
            boxOfIgusas_54000_287496Re[i].append(F(boxOfIgusas_54000_287496[i][j]))
        except:
            pass
︡94a9b3b4-c847-45b0-89ee-8093cc3dff0c︡
︠a8a051b2-2e43-4c84-9112-1606306a50bes︠
#We count occurences of each distance two hyperelliptic node to see if there are any in common.
compare123 = uniqueSetCounter(boxOfIgusasRe + boxOfIgusas_1728_54000Re + boxOfIgusas_54000_287496Re)
compare12 = uniqueSetCounter(boxOfIgusasRe + boxOfIgusas_1728_54000Re)
compare13 = uniqueSetCounter(boxOfIgusasRe + boxOfIgusas_54000_287496Re)
compare23 = uniqueSetCounter(boxOfIgusas_1728_54000Re + boxOfIgusas_54000_287496Re)
︡47c965aa-5e9a-46bc-8d7d-2b8ad35b0844︡
︠7e396622-c252-4733-b45f-314a73389b1fs︠
#We see that there is a common hyperelliptic node (104895,31610880,9746688), between our distance 1 nodes.
for i in range(len(compare123)):
    if int(compare123[i][1]) > 1 and len(compare123[i][0]) > 0:
        print(compare123[i])
︡b2ac5f89-b164-4608-9557-992f912b5ea6︡
︠185e2082-cad7-4812-ac6d-c3e4c2c90281s︠
#We see that the above common node is shared between the (1728,54000) and type 1 nodes.
for i in range(len(compare12)):
    if int(compare12[i][1]) > 1 and len(compare12[i][0]) > 0:
        print(compare12[i])
︡fdee4f94-8560-4abf-a2c4-579f7395972c︡
︠b8382481-cc80-4506-89ab-0432be1c648fs︠
#No common hyperelliptic nodes between (54000,287496) and the type 1 node.
for i in range(len(compare13)):
    if int(compare13[i][1]) > 1 and len(compare13[i][0]) > 0:
        print(compare13[i])
︡9ebfea84-c358-4689-b19b-13481623d253︡
︠9704d54b-d307-49ad-8b94-52a5027adfd1s︠
#No common hyperelliptic nodes between (54000,287496) and (1728,54000)
for i in range(len(compare23)):
    if int(compare23[i][1]) > 1 and len(compare23[i][0]) > 0:
        print(compare23[i])
︡3ff3d4ab-af6c-4807-8b73-f50d7e5249df︡
︠9ad0d100-0ca4-4cc7-a458-ed7ef3551bfcs︠
#Since [104895, 31610880, 9746688] appears in a four-cycle, we need its RA group so we can calculate dual edge weights.  We observe it has RA of C2, and hence is a type 1 curve.  As such, the dual edge to the original type one curve has weight 2, and the dual edge to (1728,54000) has weight 1.
FR = [n[0] for n in prodRoots1728_54000_5.roots()]
automorphismDetector([n for n in specifiedReverseRosenhain6LFT(FR[0],FR[1],FR[2],FR[3],FR[4],FR[5])])
︡c6233b3a-d877-4a13-a3a2-30a416cdde42︡
︠ef857c9a-c560-4e3f-a526-2c269a566f65s︠
#We now do the same for the two curves that are neighbors to the Type 1 curve and (54000,287496).
FR = [n[0] for n in prodRoots54000_287496_3.roots()]
automorphismDetector([n for n in specifiedReverseRosenhain6LFT(FR[0],FR[1],FR[2],FR[3],FR[4],FR[5])])
FR = [n[0] for n in prodRoots54000_287496_5.roots()]
automorphismDetector([n for n in specifiedReverseRosenhain6LFT(FR[0],FR[1],FR[2],FR[3],FR[4],FR[5])])
︡48e2f423-c313-4ee3-940e-8bf8478767c6︡
︠837a4dfc-0f79-4c4f-90ec-51eb6541bc53s︠
FR = [n[0] for n in prodRoots54000_287496_0.roots()]
automorphismDetector([n for n in specifiedReverseRosenhain6LFT(FR[0],FR[1],FR[2],FR[3],FR[4],FR[5])])
FR = [n[0] for n in prodRoots54000_287496_1.roots()]
automorphismDetector([n for n in specifiedReverseRosenhain6LFT(FR[0],FR[1],FR[2],FR[3],FR[4],FR[5])])
FR = [n[0] for n in prodRoots54000_287496_2.roots()]
automorphismDetector([n for n in specifiedReverseRosenhain6LFT(FR[0],FR[1],FR[2],FR[3],FR[4],FR[5])])
FR = [n[0] for n in prodRoots54000_287496_4.roots()]
automorphismDetector([n for n in specifiedReverseRosenhain6LFT(FR[0],FR[1],FR[2],FR[3],FR[4],FR[5])])
︡ce32fab9-cc41-4e62-80be-1577010e3c3a︡
︠1dba4d3a-aad6-493c-8819-68a5bf67f3e6︠









