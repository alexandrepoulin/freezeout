import numpy as np
import random

##################################################
##################################################
###These are helper functions for randomization###
##################################################
##################################################

##returns a sigmavijkl filled with zeros. nf is the number of fields
def zeroSigv(nf):
    return np.zeros((Nf,Nf,Nf,Nf))

##returns a gammaijk filled with zeros. nf is the number of fields
def zeroGam(nf):
    return np.zeros((Nf,Nf,Nf))

##randomly fills the sigmavijkl with random numbers between
##minX and maxX.  nf is the number of fields
def randomlyFilledSigv(nf, minX=0, maxX=1):
    sigv = zeroSigv(nf)
    for i in range(nf):
        for j in range(i,nf):
            for k in range(nf):
                for l in range(k,nf):
                    init = set([i,j])
                    fin = set([k,l])
                    if init == fin:
                        continue
                    sigv[i][j][k][l] = random.uniform(minX,maxX)
                    if i!=j:
                        sigv[j][i][k][l] = 0
                        sigv[j][i][l][k] = 0
                    if l!=k:
                        sigv[i][j][l][k] = 0
                    
    return sigv

##randomly fills the sigmavijkl with random numbers between
##minX and maxX keeping the cross section 0 if the initial
##state has a SM field.  nf is the number of fields
def randomlyFillSigvNoInitialSM(nf, minX=0, maxX=1):
    sigv = zeroSigv(nf)
    for i in range(nf-1):
        for j in range(i,nf-1):
            for k in range(nf):
                for l in range(k,nf):
                    init = set([i,j])
                    fin = set([k,l])
                    if init == fin:
                        continue
                    sigv[i][j][k][l] = random.uniform(minX,maxX)
                    if i!=j:
                        sigv[j][i][k][l] = 0
                        sigv[j][i][l][k] = 0
                    if l!=k:
                        sigv[i][j][l][k] = 0
    return sigv

def randomlyFillMasses(nf, minM=0, maxM=1):
    m=[]
    for i in range(nf-1):
        m.append(random.uniform(minM,maxM))
    return m

##################################################
##################################################
##################### INPUTS #####################
##################################################
##################################################


##For Nx, enter the number of Dark Matter Fields.
##
##For sigvijkl, use one of the following functions (Nf should be passed
##as the nf argument of the functions):
##    -zeroSigv(nf) for everything to be zero (you can then enter
##        your own values after
##    -randomlyFilledSigv(nf, minX=0, maxX=1) if you want every
##        entry to be random between minX and maxX
##    -randomlyFillSigvNoInitialSM(nf, minX=0, maxX=1) if you
##        want every entry to be random, except for initial state
##        SM fields which will give 0.
##                                             
##For the masses, you can manually enter a list of lenght Nx or
##you can use randomlyFillMasses(nf, minM=0, maxM=1) which will
##return a list of lenght nf of random numbers between minM and maxM.
##The g factor is the number of degrees of freedom for each DM species
##For the step size, you want a smaller step size if you have large
##sigmav or many fields. Usually between 1e-4 to 1e-5 (maybe 1e-6, but
##that will take a long time).
##                                             
##For usePreCalc, there should be a file called YeqVals.txt. This file
##contains precalculated values of Yeq. The list of values are as followed:
##[ Yeq(0.05),
##  Yeq( 250*0.5( 1/(2500-(0)*1e-4)+1/(2500-(1)*1e-4) ) ),
##  ...,
##  Yeq( 250*0.5( 1/(2500-(n)*1e-4)+1/(2500-(n+1)*1e-4) ) ),
##  ...,
##  Yeq( 250*0.5( 1/(2500-(24996676)*1e-4)+1/(2500-(24996677)*1e-4) ) ),
##  0.0]
##The last bin is zero because at that point, a double cannot hold anything
##with a smaller absolute value. The first bin should also never be used.
                                             

outputFile=r"..\inputFiles\inputDataSingleCases9.txt"
wave =1

pathToSigmav = r"C:/Users/AlexandrePoulin/Documents/relic/freezeout_v2.5/inputFiles/mathematicaOutMod.txt"
                    
Nx          = 4##number dark matter fields
Nf          = Nx+1   ##total number of fields
sigvijkl    = zeroSigv(Nf) ##thermally average cross section (s-wave)
sigvijklp   = zeroSigv(Nf) ##thermally average cross section (p-wave)
gammaijk    = zeroGam(Nf) ##decay widths to final states
masses      = randomlyFillMasses(Nf, minM=50, maxM=200) ##masses for the dark matter fields
masses      = [180,150,120,100]
nDiff       = [0,0,0,0]
##masses      = [4.3]
gFactors    = [1,1,1,1]
stepSize1   = 0.000002##The step size for the numerical procedure (x<secStepCut)
stepSize2   = 0.0001##The step size for the numerical procedure (x>secStepCut)
smallestX   = 15

usePreCalc  = 0 ## 1 for True, 0 for false. Determines whether to use
                ## the precalculated values for Yeq. Importing the file
                ## takes about 19 or so seconds, so it's not worth it if you
                ## don't have a small step size or many DM fields. It
                ## also is an approximation (a good one though).
useRK4      = 0 ## If 1, it will use the Runga-Kutta method,
                ## if 0, it will use the Euler method.

useW        = 0


##enter values for sigma manually if you want in pb. For example:
##sigvijkl=zeroSigv(Nf)


sigvijkl[0][0][4][4]=0.1
sigvijkl[1][1][4][4]=0.1
sigvijkl[2][2][4][4]=0.1
sigvijkl[3][3][4][4]=0.1

sigvijkl[0][1][2][3]=0.1
sigvijkl[0][2][1][3]=0.1
sigvijkl[0][3][1][2]=0.1


##gamma should be between 10^(-9) to 10^(-10) if nothing is producing it
##gammaijk[1][2][2] = 10**-8


##sigvijkl[1][1][5][5]=0.1
##sigvijkl[2][2][5][5]=1.0
##sigvijkl[3][3][5][5]=10.0
##sigvijkl[4][4][5][5]=100.0
##
##sigvijkl[1][3][5][5]=10.0
##sigvijkl[1][3][5][5]=5.0
##sigvijkl[2][3][5][5]=10.0
##sigvijkl[3][4][5][5]=10.0




##################################################
##################################################
##################### Checks #####################
##################################################
##################################################

if len(masses)!= Nx:
    print("ERROR: Mass list lenght inconsistent with number of dark matter fields")
    exit()
if len(sigvijkl)!= Nf:
    print("ERROR: sigvijkl list lenght inconsistent with number of fields")
    exit()
for i in range(Nf):
    if len(sigvijkl[i]) != Nf:
        print("ERROR: sigvijkl["+str(i)+"] list lenght inconsistent with number of fields")
        exit()
    else:
        for j in range(Nf):
            if len(sigvijkl[i][j]) != Nf:
                print("ERROR: sigvijkl["+str(i)+"]["+str(j)+"] list lenght inconsistent with number of fields")
                exit()
            else:
                for k in range(Nf):
                    if len(sigvijkl[i][j][k]) != Nf:
                        print("ERROR: sigvijkl["+str(i)+"]["+str(j)+"]["+str(k)+"] list lenght inconsistent with number of fields")
                        exit()
if len(sigvijklp)!= Nf:
    print("ERROR: sigvijklp list lenght inconsistent with number of fields")
    exit()
for i in range(Nf):
    if len(sigvijklp[i]) != Nf:
        print("ERROR: sigvijklp["+str(i)+"] list lenght inconsistent with number of fields")
        exit()
    else:
        for j in range(Nf):
            if len(sigvijklp[i][j]) != Nf:
                print("ERROR: sigvijklp["+str(i)+"]["+str(j)+"] list lenght inconsistent with number of fields")
                exit()
            else:
                for k in range(Nf):
                    if len(sigvijklp[i][j][k]) != Nf:
                        print("ERROR: sigvijklp["+str(i)+"]["+str(j)+"]["+str(k)+"] list lenght inconsistent with number of fields")
                        exit()

##################################################
##################################################
################# Write To File ##################
##################################################
##################################################                     

file = open(outputFile,'w')
file.write("##Use wave approximation to sigmav##\n" + str(wave))
file.write("\n\n##Number Of Dark Matter Fields##\n")
file.write(str(Nx) + "\n\n##Masses Of Dark Matter Fields##\n")
for i in range(Nx):
    file.write(str(masses[i])+" ")
file.write("\n\n##n_particle - n_antiparticle, positive if more, negative if less##\n")
for i in range(Nx):
    file.write(str(nDiff[i])+" ")




if wave==1:
    file.write("\n\n##Termally Averaged cross section matrices (i j -> k l) (s-wave)##\n")
    file.write("##i and j are the rows and columns of the BIG matrix##\n")
    file.write("##k and l are the rows and columns of the SMALLER matrices##\n\n")
    for i in range(Nf):
        for k in range(Nf):
            for j in range(Nf):
                for l in range(Nf):
                    file.write(str(sigvijkl[i][j][k][l])+"\t")
                file.write("\t")
            file.write("\n")
        file.write("\n")
        
    file.write("##same thing but for p-wave##\n\n")
    for i in range(Nf):
        for k in range(Nf):
            for j in range(Nf):
                for l in range(Nf):
                    file.write(str(sigvijklp[i][j][k][l])+"\t")
                file.write("\t")
            file.write("\n")
        file.write("\n")

    file.write("##decay widths##\n\n")
    for i in range(Nf):
        for k in range(Nf):
            for j in range(Nf):
                file.write(str(gammaijk[i][j][k])+"\t")
            file.write("\n")
        file.write("\n")
else:
    file.write("\n\n##Path to sigmav file##\n"+pathToSigmav+ "\n\n")


file.write("##Step Size for Numerical Procedure##\n"+str(stepSize1)+'\n' + str(stepSize2))
file.write("\n\n##smallest possible starting value for x##\n" +str(smallestX))
file.write("\n\n##Use Pre Calculated values for Yeq##\n" +str(usePreCalc))
file.write("\n\n##Use RK4##\n" +str(useRK4))
file.write("\n\n##g factors##\n")
for i in range(Nx):
    file.write(str(gFactors[i])+" ")
file.write("\n\n##Use W instead of Y##\n" +str(useW))
file.close()
