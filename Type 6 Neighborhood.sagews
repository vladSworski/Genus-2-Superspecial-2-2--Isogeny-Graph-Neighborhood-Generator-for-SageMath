︠6da29cc0-dc06-4b2a-8cbc-e33245db1914s︠
#Here we load in the appropriate exterior files for this project.
attach("VTools.sage","LGGUtilNew.sage","LGGClass.sage","LGGIsog.sage")
#We also setup an index for quadratic splittings now.
Index = pairMe(6)
︡fb96dc1c-2a6e-4df9-b8fa-576b1d6e4276︡
︠c08fc458-f0c5-45b8-81c9-e56d0079e83fs︠
#The top level field is best described with the notion of the square root of 5.
setupFieldLevel0.<y0> = QQ[]
#Define an intermediate field with a0 as the square root of 5.
interFieldLevel0.<a0> = setupFieldLevel0.quotient(y0^2-5)
#The golden ratio and its conjugate can be defined from the square root of 5.
phi = (1+a0)/2; cphi = (1-a0)/2
#We now define the top necessary field, which includes the fifth roots of unity.
setupFieldLevel1.<y1> = interFieldLevel0[]
interFieldLevel1.<a1> = setupFieldLevel1.quotient(y1^2+phi*y1+1)

#Calculating Richelot Isogenies relies on a Lie Bracket, we define that here so we can call it quickly whenever.  Note that the "d" function takes derivatives and the "Lie" function calculates the Lie Bracket for a quadratic pairing.
def d(poly):
    return(diff(poly, x))

def Lie(A,B):
    return([B[1]*A[2]-B[2]*A[1],B[2]*A[0]-B[0]*A[2],B[0]*A[1]-B[1]*A[0]])
︡5fd59a0a-8c4c-432a-8462-e6ecef51ce9a︡
︠a9100f54-b693-4790-8353-b21f71bd2285s︠
#We setup a polynomial field in x over the previous intermediate field.
polyFieldLevel1.<x> = interFieldLevel1[]
#We define the first kernel for consideration.  This one pairs a1 with a1^2 and a1^3 with a1^4.  R0 implies this is for the node at radius 0 from the origin node and N0 implies we are considering this to be the 0th such kernel.
kerR0N0 = [(x-1),(x-a1)*(x-a1^2),(x-a1^3)*(x-a1^4)]
#We take the derivative
d_kerR0N0 = [d(n) for n in kerR0N0]
#And then we take the image of the Richelot Isogeny.
im_kerR0N0 = Lie(kerR0N0,d_kerR0N0)
#When this cell is run, the below call will print the three quadratic factors of the image, factoring out their leading coefficients.
[factor(n) for n in im_kerR0N0]
#We record our three factors as:
#    x^2 + (a0+3)*x + 1
#    x^2 - 2*x + a0*a1 + 1
#    x^2 - 2*x - a0*a1 - 1/2*a0 - 3/2
︡5ecb3251-4380-4d0f-8c8f-c9bc57f170bd︡
︠48e8ddcf-3bd7-44c8-8b50-0068fef1b86es︠
#We now need to factor these into linear terms.  We establish a couple field extensions to do this.
setupFieldLevel2.<y2> = interFieldLevel1[]
interFieldLevel2.<a2> = setupFieldLevel2.quotient(y2^2+(a0+3)*y2+1)
setupFieldLevel3.<y3> = interFieldLevel2[]
interFieldLevel3.<a3> = setupFieldLevel3.quotient(y3^2 - 2*y3 + a0*a1 + 1)

#We then expand our polynomial field to encompass this new info:
polyFieldLevel2.<x> = interFieldLevel3[]

#And attempt factoring into linear terms again:
[factor(polyFieldLevel2(n)) for n in im_kerR0N0]
︡5375112a-5387-4ba2-b049-fc977f1529a2︡
︠03e776ef-9639-44d1-bc73-c383dbce01fas︠
#Now that we know it is possible, we prepare the linear roots for being used in quadratic splittings.
prod_kerR0N0 = polyFieldLevel2(im_kerR0N0[0])*polyFieldLevel2(im_kerR0N0[1])*polyFieldLevel2(im_kerR0N0[2])
rawRoots_kerR0N0 = prod_kerR0N0.roots()
roots_kerR0N0 = [n[0] for n in rawRoots_kerR0N0]
roots_kerR0N0

#We also briefly check what the igusa invariant of our curve is.
hyper_kerR0N0 = HyperellipticCurve(prod_kerR0N0)
igusa_kerR0N0 = hyper_kerR0N0.absolute_igusa_invariants_kohel()
igusa_kerR0N0
#This curve has igusa invariants: (-2035611/4*a0 + 4598775/4, -199559376/5*a0 + 439744464/5, -257829804/25*a0 + 113452812/5)
#  We will refer to this curve by the designation: A_bar_0
︡1cb29d0b-769d-4dc8-9ae4-79c77f53cb9f︡
︠bd6bfead-00fb-4bdc-aaa3-96fc47b7b836s︠
#And we now prepare to construct all the isogenies from this curve.
kersR1N0 = [[(x-roots_kerR0N0[n[0]])*(x-roots_kerR0N0[n[1]]),(x-roots_kerR0N0[n[2]])*(x-roots_kerR0N0[n[3]]),(x-roots_kerR0N0[n[4]])*(x-roots_kerR0N0[n[5]])] for n in Index]
d_kersR1N0 = [[d(m) for m in n] for n in kersR1N0]
im_kersR1N0 = [[Lie(kersR1N0[i],d_kersR1N0[i])] for i in range(len(kersR1N0))]
prod_kersR1N0 = [[n[0][0]*n[0][1]*n[0][2]] for n in im_kersR1N0]

#Now we calculate the Igusa invariants of these images so we can differentiate them.
hypers_kersR1N0 = [HyperellipticCurve(n[0]) for n in prod_kersR1N0]
igusas_kersR1N0 = [n.absolute_igusa_invariants_kohel() for n in hypers_kersR1N0]
#Pulls out unique instances of the igusa invariants.  Each entry is a three-tuple containing the invariants, the number of instances of that invariant, and the indices for those instances.
uniqueIgusas_kersR1N0 = uniqueSetCounter(igusas_kersR1N0)
#This list only includes the unique invariants.
justIgusas_uniqueIgusas_kersR1N0 = [n[0] for n in uniqueIgusas_kersR1N0]
#This list includes just the number and instances.
simplified_uniqueIgusas_kersR1N0 = [n[1:] for n in uniqueIgusas_kersR1N0]
#Here we both list out the above data telling us how many instances of unique curves there are and their outgoing weights, but also we check the weight of the backwards edge.
simplified_uniqueIgusas_kersR1N0; uniqueIgusas_kersR1N0[0]
#Test whether there are any self-loops.  True means there is a self-loop, false means there isn't.
(-2035611/4*a0 + 4598775/4, -199559376/5*a0 + 439744464/5, -257829804/25*a0 + 113452812/5) in justIgusas_uniqueIgusas_kersR1N0

#Takeaways: Exactly one curve returns to the starting type 6 node.  There are no self loops.
#  Otherwise there are 10 image curves with a weight one edge
#        and there are  2 image curves with a weight two edge
︡585de778-3991-4a12-831b-8005238607e8︡
︠dbd2392e-a518-4081-b60a-12e127f4f131s︠
#We define the second kernel for consideration.  This one pairs a1 with a1^3 and a1^2 with a1^4.  R0 implies this is for the node at radius 0 from the origin node and N1 implies we are considering this to be the 1st such kernel.
kerR0N1 = [(x-1),(x-a1)*(x-a1^3),(x-a1^2)*(x-a1^4)]
d_kerR0N1 = [d(n) for n in kerR0N1]
im_kerR0N1 = Lie(kerR0N1,d_kerR0N1)
#When this cell is run, the below call will print the three quadratic factors of the image, factoring out their leading coefficients.
[factor(n) for n in im_kerR0N1]
#We record our three factors as:
#    (x + (-a1 - 1/2*a0 + 1/2)*a3 + a1 + 1) * (x + (a1 + 1/2*a0 - 1/2)*a3 - a1 - a0 + 2)
#    (x + ((-1/2*a0 + 1/2)*a1 - 1)*a2 + (-1/2*a0 - 1/2)*a1 - 1/2*a0 - 5/2) * (x + ((1/2*a0 - 1/2)*a1 + 1)*a2 + (1/2*a0 + 1/2)*a1 + 1/2*a0 + 1/2)
#    (x + ((-1/2*a0 + 1/2)*a1)*a2 + (-1/2*a0 - 1/2)*a1 - 1) * (x + ((1/2*a0 - 1/2)*a1)*a2 + (1/2*a0 + 1/2)*a1 - 1)
#Rather noticably, these are all linear factors, and hence we don't need any new field extensions.
︡b27a59b3-84a3-4f65-aba6-627f90e51001︡
︠cd62556c-769f-4294-85b4-9c3ccdabb762s︠
#We prepare the linear roots for being used in quadratic splittings.
prod_kerR0N1 = polyFieldLevel2(im_kerR0N1[0])*polyFieldLevel2(im_kerR0N1[1])*polyFieldLevel2(im_kerR0N1[2])
rawRoots_kerR0N1 = prod_kerR0N1.roots()
roots_kerR0N1 = [n[0] for n in rawRoots_kerR0N1]
roots_kerR0N1

#We also briefly check what the igusa invariant of our curve is.
hyper_kerR0N1 = HyperellipticCurve(prod_kerR0N1)
igusa_kerR0N1 = hyper_kerR0N1.absolute_igusa_invariants_kohel()
igusa_kerR0N1
#This curve has igusa invariants: (2035611/4*a0 + 4598775/4, 199559376/5*a0 + 439744464/5, 257829804/25*a0 + 113452812/5)
#  We will refer to this curve by the designation: A_bar_1.  Note that it is not the same curve as A_bar_0
︡5d74ab65-3f81-4b2a-908c-6650c694f084︡
︠49840908-fdbc-4fa1-b10d-7287d54099bds︠
#And we now prepare to construct all the isogenies from this curve.
kersR1N1 = [[(x-roots_kerR0N1[n[0]])*(x-roots_kerR0N1[n[1]]),(x-roots_kerR0N1[n[2]])*(x-roots_kerR0N1[n[3]]),(x-roots_kerR0N1[n[4]])*(x-roots_kerR0N1[n[5]])] for n in Index]
d_kersR1N1 = [[d(m) for m in n] for n in kersR1N1]
im_kersR1N1 = [[Lie(kersR1N1[i],d_kersR1N1[i])] for i in range(len(kersR1N1))]
prod_kersR1N1 = [[n[0][0]*n[0][1]*n[0][2]] for n in im_kersR1N1]

#Now we calculate the Igusa invariants of these images so we can differentiate them.
hypers_kersR1N1 = [HyperellipticCurve(n[0]) for n in prod_kersR1N1]
igusas_kersR1N1 = [n.absolute_igusa_invariants_kohel() for n in hypers_kersR1N1]
#Pulls out unique instances of the igusa invariants.  Each entry is a three-tuple containing the invariants, the number of instances of that invariant, and the indices for those instances.
uniqueIgusas_kersR1N1 = uniqueSetCounter(igusas_kersR1N1)
#This list only includes the unique invariants.
justIgusas_uniqueIgusas_kersR1N1 = [n[0] for n in uniqueIgusas_kersR1N1]
#This list includes just the number and instances.
simplified_uniqueIgusas_kersR1N1 = [n[1:] for n in uniqueIgusas_kersR1N1]
#Here we both list out the above data telling us how many instances of unique curves there are and their outgoing weights, but also we check the weight of the backwards edge.
simplified_uniqueIgusas_kersR1N1; uniqueIgusas_kersR1N1[6]
#Test whether there are any self-loops.  True means there is a self-loop, false means there isn't.
(2035611/4*a0 + 4598775/4, 199559376/5*a0 + 439744464/5, 257829804/25*a0 + 113452812/5) in justIgusas_uniqueIgusas_kersR1N1
#Also test whether the two type A_bar nodes connect.
(-2035611/4*a0 + 4598775/4, -199559376/5*a0 + 439744464/5, -257829804/25*a0 + 113452812/5) in justIgusas_uniqueIgusas_kersR1N1

#Takeaways: Exactly one curve returns to the starting type 6 node.  There are no self loops.  The two type A_bar nodes don't connect.
#  Otherwise there are 10 image curves with a weight one edge
#        and there are  2 image curves with a weight two edge
#In fact, the above setup is how we define a type A_bar node in this graph.
︡f8a1d2e2-85ee-4854-8224-91de036583bc︡
︠2eee38ac-c829-44d3-a238-68a8a8069411s︠
#We define the third kernel for consideration.  This one pairs a1 with a1^4 and a1^2 with a1^3.  R0 implies this is for the node at radius 0 from the origin node and N2 implies we are considering this to be the 2nd such kernel.
kerR0N2 = [(x-1),(x-a1)*(x-a1^4),(x-a1^2)*(x-a1^3)]
d_kerR0N2 = [d(n) for n in kerR0N2]
im_kerR0N2 = Lie(kerR0N2,d_kerR0N2)
#When this cell is run, the below call will print the three quadratic factors of the image, factoring out their leading coefficients.
[factor(n) for n in im_kerR0N2]
#We record our three factors as:
#    (x - 1) * (x + 1)
#    (x + (((-3/10*a0 + 1/2)*a1 - 1/5*a0)*a2 - 1/5*a0*a1 - 3/10*a0 - 1/2)*a3 + ((3/10*a0 - 1/2)*a1 + 1/5*a0)*a2 + 1/5*a0*a1 + 3/10*a0 - 1/2) * (x + (((3/10*a0 - 1/2)*a1 + 1/5*a0)*a2 + 1/5*a0*a1 + 3/10*a0 + 1/2)*a3 + ((-3/10*a0 + 1/2)*a1 - 1/5*a0)*a2 - 1/5*a0*a1 - 3/10*a0 - 3/2)
#    (x + (((-1/10*a0 + 1/2)*a1 + 1/10*a0 + 1/2)*a2 + (1/10*a0 + 1/2)*a1 + 2/5*a0 + 1)*a3 + ((1/10*a0 - 1/2)*a1 - 1/10*a0 - 1/2)*a2 + (-1/10*a0 - 1/2)*a1 - 2/5*a0 - 2) * (x + (((1/10*a0 - 1/2)*a1 - 1/10*a0 - 1/2)*a2 + (-1/10*a0 - 1/2)*a1 - 2/5*a0 - 1)*a3 + ((-1/10*a0 + 1/2)*a1 + 1/10*a0 + 1/2)*a2 + (1/10*a0 + 1/2)*a1 + 2/5*a0)
#Rather noticably, these are all linear factors, and hence we don't need any new field extensions.
︡f3642e12-80d6-47c7-9e2d-965790cf6d72︡
︠3465267b-4e68-4608-9ba2-19e2e54b4f64s︠
#We prepare the linear roots for being used in quadratic splittings.
prod_kerR0N2 = polyFieldLevel2(im_kerR0N2[0])*polyFieldLevel2(im_kerR0N2[1])*polyFieldLevel2(im_kerR0N2[2])
rawRoots_kerR0N2 = prod_kerR0N2.roots()
roots_kerR0N2 = [n[0] for n in rawRoots_kerR0N2]
roots_kerR0N2

#We also briefly check what the igusa invariant of our curve is.
hyper_kerR0N2 = HyperellipticCurve(prod_kerR0N2)
igusa_kerR0N2 = hyper_kerR0N2.absolute_igusa_invariants_kohel()
igusa_kerR0N2
#This curve has igusa invariants: (2824875/2, 656100000, 203391000)
#  We will refer to this curve by the designation: A_hat.
︡0f57d50e-18b1-44ce-8c4f-4d495ecc8c76︡
︠f28b4bfa-03f8-48ef-afb7-63323182e575s︠
#And we now prepare to construct all the isogenies from this curve.
kersR1N2 = [[(x-roots_kerR0N2[n[0]])*(x-roots_kerR0N2[n[1]]),(x-roots_kerR0N2[n[2]])*(x-roots_kerR0N2[n[3]]),(x-roots_kerR0N2[n[4]])*(x-roots_kerR0N2[n[5]])] for n in Index]
d_kersR1N2 = [[d(m) for m in n] for n in kersR1N2]
im_kersR1N2 = [[Lie(kersR1N2[i],d_kersR1N2[i])] for i in range(len(kersR1N2))]
prod_kersR1N2 = [[n[0][0]*n[0][1]*n[0][2]] for n in im_kersR1N2]

#Now we calculate the Igusa invariants of these images so we can differentiate them.
hypers_kersR1N2 = [HyperellipticCurve(n[0]) for n in prod_kersR1N2]
igusas_kersR1N2 = [n.absolute_igusa_invariants_kohel() for n in hypers_kersR1N2]
#Pulls out unique instances of the igusa invariants.  Each entry is a three-tuple containing the invariants, the number of instances of that invariant, and the indices for those instances.
uniqueIgusas_kersR1N2 = uniqueSetCounter(igusas_kersR1N2)
#This list only includes the unique invariants.
justIgusas_uniqueIgusas_kersR1N2 = [n[0] for n in uniqueIgusas_kersR1N2]
#This list includes just the number and instances.
simplified_uniqueIgusas_kersR1N2 = [n[1:] for n in uniqueIgusas_kersR1N2]
#Here we both list out the above data telling us how many instances of unique curves there are and their outgoing weights, but also we check the weight of the backwards edge.
simplified_uniqueIgusas_kersR1N2; uniqueIgusas_kersR1N2[2]
#Test whether there are any self-loops.  True means there is a self-loop, false means there isn't.
(2824875/2, 656100000, 203391000) in justIgusas_uniqueIgusas_kersR1N2
#Also test whether the type A_hat node connects to either type A_bar node.
(2035611/4*a0 + 4598775/4, 199559376/5*a0 + 439744464/5, 257829804/25*a0 + 113452812/5) in justIgusas_uniqueIgusas_kersR1N2
(-2035611/4*a0 + 4598775/4, -199559376/5*a0 + 439744464/5, -257829804/25*a0 + 113452812/5) in justIgusas_uniqueIgusas_kersR1N2

#Takeaways: Exactly one curve returns to the starting type 6 node.  There are no self loops.  The type A_bar nodes do not connect to the type A_hat node.
#  Otherwise there are 14 image curves with a weight one edge
#        and there are  0 image curves with a weight two edge
#This is different than the type A_bar nodes, hence the different name.
︡b28cbff8-c189-4575-9e96-fa9302fcb3d9︡
︠a6a67a33-1c2c-4f1c-a532-e93f6aa673c5s︠
#We now check for common distance-two nodes in the graph.
uniqueInvarMergeList_01 = justIgusas_uniqueIgusas_kersR1N0 + justIgusas_uniqueIgusas_kersR1N1
uniqueInvarMergeList_02 = justIgusas_uniqueIgusas_kersR1N0 + justIgusas_uniqueIgusas_kersR1N2
uniqueInvarMergeList_12 = justIgusas_uniqueIgusas_kersR1N1 + justIgusas_uniqueIgusas_kersR1N2
overlappedCounts_01 = uniqueSetCounter(uniqueInvarMergeList_01)
overlappedCounts_02 = uniqueSetCounter(uniqueInvarMergeList_02)
overlappedCounts_12 = uniqueSetCounter(uniqueInvarMergeList_12)
print("Here is the data on the overlap of the two A_bar nodes.")
for i in range(len(overlappedCounts_01)):
    if overlappedCounts_01[i][1] > 1 and overlappedCounts_01[i][0] != (0,0,0):
        firstCount01 = simplified_uniqueIgusas_kersR1N0[overlappedCounts_01[i][2][0]]
        secondCount01 = simplified_uniqueIgusas_kersR1N1[overlappedCounts_01[i][2][1]-len(justIgusas_uniqueIgusas_kersR1N0)]
        if igusas_kersR1N0[firstCount01[1][0]] == igusas_kersR1N1[secondCount01[1][0]] and igusas_kersR1N0[firstCount01[1][0]] == overlappedCounts_01[i][0] and igusas_kersR1N1[secondCount01[1][0]] == overlappedCounts_01[i][0]:
            
            print(firstCount01, secondCount01, overlappedCounts_01[i][0])
        else:
            print("Error, double check code.")
            
print("Here is the data on the overlap of the A_bar_0 and A_hat nodes.")
for i in range(len(overlappedCounts_02)):
    if overlappedCounts_02[i][1] > 1 and overlappedCounts_02[i][0] != (0,0,0):
        firstCount02 = simplified_uniqueIgusas_kersR1N0[overlappedCounts_02[i][2][0]]
        secondCount02 = simplified_uniqueIgusas_kersR1N2[overlappedCounts_02[i][2][1]-len(justIgusas_uniqueIgusas_kersR1N0)]
        if igusas_kersR1N0[firstCount02[1][0]] == igusas_kersR1N2[secondCount02[1][0]] and igusas_kersR1N0[firstCount02[1][0]] == overlappedCounts_02[i][0] and igusas_kersR1N2[secondCount02[1][0]] == overlappedCounts_02[i][0]:
            
            print(firstCount02, secondCount02, overlappedCounts_02[i][0])
        else:
            print("Error, double check code.")
            
print("Here is the data on the overlap of the A_bar_1 and A_hat nodes.")
for i in range(len(overlappedCounts_12)):
    if overlappedCounts_12[i][1] > 1 and overlappedCounts_12[i][0] != (0,0,0):
        firstCount12 = simplified_uniqueIgusas_kersR1N1[overlappedCounts_12[i][2][0]]
        secondCount12 = simplified_uniqueIgusas_kersR1N2[overlappedCounts_12[i][2][1]-len(justIgusas_uniqueIgusas_kersR1N1)]
        if igusas_kersR1N1[firstCount12[1][0]] == igusas_kersR1N2[secondCount12[1][0]] and igusas_kersR1N1[firstCount12[1][0]] == overlappedCounts_12[i][0] and igusas_kersR1N2[secondCount12[1][0]] == overlappedCounts_12[i][0]:
            
            print(firstCount12, secondCount12, overlappedCounts_12[i][0])
        else:
            print("Error, double check code.")

#Here's our observations:
#  There are two nodes which are adjacent to both A_bars and the A_hat node.  Every outgoing edge to them is weight 1:
#     [1, [1]] [1, [3]], [1,[1]]
#     [1, [2]] [1, [0]], [1,[0]]
#  We call these type "sigma" nodes.
#
#  There are two nodes which are adjacent to just A_bar_0 and A_hat.  The outgoing edges from A_bar_0 are weight two, and the outgoing
#  edges from A_hat are weight one.
#     [2, [3, 11]] [1, [5]]
#     [2, [6, 14]] [1, [12]]
#  There are two nodes which are adjacent to just A_bar_1 and A_hat.  The outgoing edges from A_bar_1 are weight two, and the outgoing
#  edges from A_hat are weight one.
#     [2, [7, 9]] [1, [10]]
#     [2, [8, 12]] [1, [7]]
#  We call these type "kappa" nodes.
#
#  In the same order they appear above, these nodes invariants are:
#  1) (56531034*a0 + 126411030, 13692229632*a0 + 153086920704/5, 120920174784/25*a0 + 10815607296)
#  2) (-56531034*a0 + 126411030, -13692229632*a0 + 153086920704/5, -120920174784/25*a0 + 10815607296)
#  3) (((((-227571687/2*a0 + 508770801/2)*a1 - 43486065*a0 + 97113492)*a2 + (-43486065*a0 + 97113492)*a1 - 33344703/2*a0 + 73910151/2)*a3 + ((227571687/2*a0 - 508770801/2)*a1 + 43486065*a0 - 97113492)*a2 + (43486065*a0 - 97113492)*a1 + 403167375/2*a0 - 900212211/2, (((321021646848/25*a0 - 143574415488/5)*a1 + 122596431552/25*a0 - 54850799808/5)*a2 + (122596431552/25*a0 - 54850799808/5)*a1 + 46767647808/25*a0 - 20977983936/5)*a3 + ((-321021646848/25*a0 + 143574415488/5)*a1 - 122596431552/25*a0 + 54850799808/5)*a2 + (-122596431552/25*a0 + 54850799808/5)*a1 - 568120454208/25*a0 + 50861711808, (((138527711784/25*a0 - 309785943096/25)*a1 + 52898596128/25*a0 - 118359635184/25)*a2 + (52898596128/25*a0 - 118359635184/25)*a1 + 806723064*a0 - 45292962456/25)*a3 + ((-138527711784/25*a0 + 309785943096/25)*a1 - 52898596128/25*a0 + 118359635184/25)*a2 + (-52898596128/25*a0 + 118359635184/25)*a1 - 245161807704/25*a0 + 548700821496/25))
#  4) ((((227571687/2*a0 - 508770801/2)*a1 + 43486065*a0 - 97113492)*a2 + (43486065*a0 - 97113492)*a1 + 33344703/2*a0 - 73910151/2)*a3 + ((-227571687/2*a0 + 508770801/2)*a1 - 43486065*a0 + 97113492)*a2 + (-43486065*a0 + 97113492)*a1 + 336477969/2*a0 - 752391909/2, (((-321021646848/25*a0 + 143574415488/5)*a1 - 122596431552/25*a0 + 54850799808/5)*a2 + (-122596431552/25*a0 + 54850799808/5)*a1 - 46767647808/25*a0 + 20977983936/5)*a3 + ((321021646848/25*a0 - 143574415488/5)*a1 + 122596431552/25*a0 - 54850799808/5)*a2 + (122596431552/25*a0 - 54850799808/5)*a1 - 474585158592/25*a0 + 212352591168/5, (((-138527711784/25*a0 + 309785943096/25)*a1 - 52898596128/25*a0 + 118359635184/25)*a2 + (-52898596128/25*a0 + 118359635184/25)*a1 - 806723064*a0 + 45292962456/25)*a3 + ((138527711784/25*a0 - 309785943096/25)*a1 + 52898596128/25*a0 - 118359635184/25)*a2 + (52898596128/25*a0 - 118359635184/25)*a1 - 204825654504/25*a0 + 458114896584/25)
#  5) ((((-53627427/2*a0 - 120316833/2)*a1 - 140599557/2*a0 - 314543817/2)*a2 + (-140599557/2*a0 - 314543817/2)*a1 - 184085622*a0 - 411657309)*a3 + ((53627427/2*a0 + 120316833/2)*a1 + 140599557/2*a0 + 314543817/2)*a2 + (140599557/2*a0 + 314543817/2)*a1 - 825714*a0 - 1493721, (((75828783744/25*a0 + 33872815872/5)*a1 + 198425215296/25*a0 + 17744723136)*a2 + (198425215296/25*a0 + 17744723136)*a1 + 519446862144/25*a0 + 232298031168/5)*a3 + ((-75828783744/25*a0 - 33872815872/5)*a1 - 198425215296/25*a0 - 17744723136)*a2 + (-198425215296/25*a0 - 17744723136)*a1 + 1905944256/25*a0 + 1032543936/5, (((32730519528/25*a0 + 73066672728/25)*a1 + 85629115656/25*a0 + 191426307912/25)*a2 + (85629115656/25*a0 + 191426307912/25)*a1 + 44831365488/5*a0 + 501212251008/25)*a3 + ((-32730519528/25*a0 - 73066672728/25)*a1 - 85629115656/25*a0 - 191426307912/25)*a2 + (-85629115656/25*a0 - 191426307912/25)*a1 + 836903664/25*a0 + 2195608032/25)
#  6) ((((53627427/2*a0 + 120316833/2)*a1 + 140599557/2*a0 + 314543817/2)*a2 + (140599557/2*a0 + 314543817/2)*a1 + 184085622*a0 + 411657309)*a3 + ((-53627427/2*a0 - 120316833/2)*a1 - 140599557/2*a0 - 314543817/2)*a2 + (-140599557/2*a0 - 314543817/2)*a1 - 368996958*a0 - 824808339, (((-75828783744/25*a0 - 33872815872/5)*a1 - 198425215296/25*a0 - 17744723136)*a2 + (-198425215296/25*a0 - 17744723136)*a1 - 519446862144/25*a0 - 232298031168/5)*a3 + ((75828783744/25*a0 + 33872815872/5)*a1 + 198425215296/25*a0 + 17744723136)*a2 + (198425215296/25*a0 + 17744723136)*a1 + 1040799668544/25*a0 + 465628606272/5, (((-32730519528/25*a0 - 73066672728/25)*a1 - 85629115656/25*a0 - 191426307912/25)*a2 + (-85629115656/25*a0 - 191426307912/25)*a1 - 44831365488/5*a0 - 501212251008/25)*a3 + ((32730519528/25*a0 + 73066672728/25)*a1 + 85629115656/25*a0 + 191426307912/25)*a2 + (85629115656/25*a0 + 191426307912/25)*a1 + 449150558544/25*a0 + 1004620110048/25)
#
#With this, we now have the weights of every important edge in the graph except for the dual edges from sigma and kappa nodes.
︡2e9bae85-3927-417c-a7d7-684dc25339a7︡
︠797e4ab3-7455-4fca-b1cc-75f42f47887cs︠
#We now establish the roots for the two type sigma nodes.  This may take a few minutes.
sig1Roots = [n[0] for n in prod_kersR1N0[1][0].roots()]
sig2Roots = [n[0] for n in prod_kersR1N0[2][0].roots()]
#And we calculate the automorphism groups for these curves.  Note that the automorphismDetector function requires a defined field called "F."  We want polyFieldLevel2 to fulfill this role.
F = polyFieldLevel2
automorphismDetector([n for n in specifiedReverseRosenhain6LFT(sig1Roots[0],sig1Roots[1],sig1Roots[2],sig1Roots[3],sig1Roots[4],sig1Roots[5])])
automorphismDetector([n for n in specifiedReverseRosenhain6LFT(sig2Roots[0],sig2Roots[1],sig2Roots[2],sig2Roots[3],sig2Roots[4],sig2Roots[5])])

#Here's our observation:
# The first automorphism group is C1
# The second automorphism group is C1
# Hence both curves are type "0" or type "A" curves.
#
# By the orbit stabilizer theorem, this forces the dual edges back to A_bar nodes to be weight 1 and the dual edges back to A_hat to also be weight 1.
︡be40e0f5-4012-417c-9233-c072a9499ea6︡
︠46eaa551-3b06-421a-a068-046628f7a2a9s︠
#We now attempt to do the same for the four type kappa nodes.
kap1Roots = [n[0] for n in prod_kersR1N0[3][0].roots()]; len(kap1Roots)
kap2Roots = [n[0] for n in prod_kersR1N0[6][0].roots()]; len(kap2Roots)
kap3Roots = [n[0] for n in prod_kersR1N1[7][0].roots()]; len(kap3Roots)
kap4Roots = [n[0] for n in prod_kersR1N1[8][0].roots()]; len(kap4Roots)
prod_kersR1N0[3][0].factor()
prod_kersR1N0[6][0].factor()
prod_kersR1N1[7][0].factor()
prod_kersR1N1[8][0].factor()
#Here we notice that we have quadratic terms that won't factor once more.  All 4 of these cases have 4 linear factors and a quadratic factor.
# 1) (x^2 + (((-2/5*a0*a1 + 2/5*a0)*a2 + (-3/5*a0 - 1)*a1 + 3/5*a0 + 1)*a3 + (2/5*a0*a1 - 2/5*a0)*a2 + (8/5*a0 + 4)*a1 - 3/5*a0 - 1)*x + (((-1/10*a0 - 1/2)*a1 + 1/10*a0 - 1/2)*a2 + (-2/5*a0 - 1)*a1 - 1/10*a0 - 1/2)*a3 + ((1/10*a0 + 1/2)*a1 - 1/10*a0 + 1/2)*a2 + (7/5*a0 + 3)*a1 + 3/5*a0 + 2)
# 2) (x^2 + (((2/5*a0*a1 - 2/5*a0)*a2 + (3/5*a0 + 1)*a1 - 3/5*a0 - 1)*a3 + (-2/5*a0*a1 + 2/5*a0)*a2 + (2/5*a0 + 2)*a1 + 3/5*a0 + 1)*x + (((1/10*a0 + 1/2)*a1 - 1/10*a0 + 1/2)*a2 + (2/5*a0 + 1)*a1 + 1/10*a0 + 1/2)*a3 + ((-1/10*a0 - 1/2)*a1 + 1/10*a0 - 1/2)*a2 + (3/5*a0 + 1)*a1 + 2/5*a0 + 1)
# 3) (x^2 + ((((-1/5*a0 + 1)*a1 - 4/5*a0 + 2)*a2 + (1/5*a0 + 1)*a1 - 1/5*a0 + 1)*a3 + ((1/5*a0 - 1)*a1 + 4/5*a0 - 2)*a2 - 6/5*a0*a1 + 6/5*a0 - 4)*x + ((1/5*a0*a1 - 1/5*a0 + 1)*a2 + (3/10*a0 + 1/2)*a1 + 1/5*a0 + 1)*a3 + (-1/5*a0*a1 + 1/5*a0 - 1)*a2 + (-4/5*a0 + 1)*a1 + 3/10*a0 - 3/2)
# 4) (x^2 + ((((1/5*a0 - 1)*a1 + 4/5*a0 - 2)*a2 + (-1/5*a0 - 1)*a1 + 1/5*a0 - 1)*a3 + ((-1/5*a0 + 1)*a1 - 4/5*a0 + 2)*a2 + (-4/5*a0 + 2)*a1 + 4/5*a0 - 2)*x + ((-1/5*a0*a1 + 1/5*a0 - 1)*a2 + (-3/10*a0 - 1/2)*a1 - 1/5*a0 - 1)*a3 + (1/5*a0*a1 - 1/5*a0 + 1)*a2 + (-1/5*a0 + 2)*a1 + 7/10*a0 + 1/2)
︡0d63542e-5bb6-4705-904e-b5cef835a7ca︡
︠6ce397a2-9e40-48b2-9061-7dde24a96672s︠
#We now set up for our final field extension.
setupFieldLevel4.<y4> = interFieldLevel3[]
interFieldLevel4.<a4> = setupFieldLevel4.quotient(y4^2 + (((-2/5*a0*a1 + 2/5*a0)*a2 + (-3/5*a0 - 1)*a1 + 3/5*a0 + 1)*a3 + (2/5*a0*a1 - 2/5*a0)*a2 + (8/5*a0 + 4)*a1 - 3/5*a0 - 1)*y4 + (((-1/10*a0 - 1/2)*a1 + 1/10*a0 - 1/2)*a2 + (-2/5*a0 - 1)*a1 - 1/10*a0 - 1/2)*a3 + ((1/10*a0 + 1/2)*a1 - 1/10*a0 + 1/2)*a2 + (7/5*a0 + 3)*a1 + 3/5*a0 + 2)
polyFieldLevel3.<x> = interFieldLevel4[]
F = polyFieldLevel3
︡24990696-2b9a-400b-86ed-f71d082ee4b1︡
︠3447b85b-5696-4293-b7d5-f34ad94d56afr︠
#And we try to factor the type kappa nodes again.  This step will take a while.
kap1Roots = [n[0] for n in F(prod_kersR1N0[3][0]).roots()]
kap2Roots = [n[0] for n in F(prod_kersR1N0[6][0]).roots()]
kap3Roots = [n[0] for n in F(prod_kersR1N1[7][0]).roots()]
kap4Roots = [n[0] for n in F(prod_kersR1N1[8][0]).roots()]
automorphismDetector([n for n in specifiedReverseRosenhain6LFT(kap1Roots[0],kap1Roots[1],kap1Roots[2],kap1Roots[3],kap1Roots[4],kap1Roots[5])])
automorphismDetector([n for n in specifiedReverseRosenhain6LFT(kap2Roots[0],kap2Roots[1],kap2Roots[2],kap2Roots[3],kap2Roots[4],kap2Roots[5])])
automorphismDetector([n for n in specifiedReverseRosenhain6LFT(kap3Roots[0],kap3Roots[1],kap3Roots[2],kap3Roots[3],kap3Roots[4],kap3Roots[5])])
automorphismDetector([n for n in specifiedReverseRosenhain6LFT(kap4Roots[0],kap4Roots[1],kap4Roots[2],kap4Roots[3],kap4Roots[4],kap4Roots[5])])

#Here's our observation:
# The first automorphism group is C1
# The second automorphism group is C1
# The third automorphism group is C1
# The fourth automorphism group is C1
# Hence all curves are type "0" or type "A" curves.
#
# By the orbit stabilizer theorem, this forces the dual edges back to A_bar nodes to be weight 2 and the dual edges back to A_hat to be weight 1.
#This is enough information to draw the entire crab graph from the paper.
︡a661262f-fde3-4651-940d-1d2a24796279︡
︠6f2518c2-da06-4002-b228-fecca16b7233︠









