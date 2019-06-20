#Script for generating initial configurations for Ethylbenzene system: a single GB disk attached to two methyl beads

import sys
import math
import numpy as np
from scipy import stats
import copy
import glob

dataFile = open("input.data","w")

box_x = 500.0
box_y = 500.0
box_z = 500.0

numMoly = 609

nAtoms = numMoly*4

#num_spacing = 10

vec2 = np.array([1.0,0.0,0.0])
vec1 = np.array([0.333,-0.943,0.0])
dr = 1.51
rad = 1.4
current = np.zeros(3)
prev = np.zeros(3)
ID = 0
#Generating all the molecular configuration information
dataFile.write("LAMMPS Description\n\n")

dataFile.write("\t"+str(nAtoms)+"\tatoms\n")
dataFile.write("\t"+str(3*numMoly)+"\tbonds\n")
dataFile.write("\t"+str(3*numMoly)+"\tangles\n")
dataFile.write("\t"+str(0*numMoly)+"\tdihedrals\n")
dataFile.write("\t"+str(0*numMoly)+"\timpropers\n")
dataFile.write("\t"+str(numMoly)+"\tellipsoids\n")

dataFile.write("\n")

dataFile.write("\t"+"4\tatom types\n")
dataFile.write("\t"+"3\tbond types\n")
dataFile.write("\t"+"3\tangle types\n")
dataFile.write("\t"+"0\tdihedral types\n")
dataFile.write("\t"+"0\timproper types\n")

dataFile.write("\n")

dataFile.write(str(0.0)+" "+str(box_x)+" xlo xhi\n")
dataFile.write(str(0.0)+" "+str(box_y)+" ylo yhi\n")
dataFile.write(str(0.0)+" "+str(box_z)+" zlo zhi\n")

dataFile.write("\n")

dataFile.write("Masses\n\n")
dataFile.write("1\t14.0\n")
dataFile.write("2\t13.0\n")
dataFile.write("3\t77.0\n")
dataFile.write("4\t0.00001\n")

dataFile.write("\n")

dataFile.write("Atoms\n\n")

flag = 0
count = 1
bodyid = 1
for i in range(numMoly):
    x = box_x*np.random.uniform()
    y = box_x*np.random.uniform()
    z = box_x*np.random.uniform()

    #Main Ellipsoid
    dataFile.write(str(count)+" 3 "+str(x)+" "+str(y)+" "+str(z) + " " +str(bodyid) + " 0.0 1 66.55\n")

    count += 1
    #Rigid Site
    dataFile.write(str(count)+" 4 "+str(x+1.4)+" "+str(y)+" "+str(z) + " " +str(bodyid) + " 0.0 0 0.0001\n")
    count += 1

    #Bonded Beads
    dataFile.write(str(count)+" 2 "+str(x+1.4+1.52)+" "+str(y)+" "+str(z) + " " +str(bodyid) + " 0.0 0 13.0\n")
    count += 1
    dataFile.write(str(count)+" 2 "+str(x+1.4+1.52+0.637)+" "+str(y+1.39)+" "+str(z) + " " +str(bodyid) + " 0.0 0 14.0\n")
    count += 1

    bodyid += 1

    if (bodyid > numMoly):
        flag = 1
        break

dataFile.write("\n")

dataFile.write("Bonds\n\n")

count = 1
for i in range(numMoly):

    m = i*4+1

    dataFile.write(str(count)+" 1 "+str(m)+" "+str(m+1)+"\n")
    count += 1
    dataFile.write(str(count)+" 2 "+str(m+1)+" "+str(m+2)+"\n")
    count += 1
    dataFile.write(str(count)+" 3 "+str(m+2)+" "+str(m+3)+"\n")
    count += 1


dataFile.write("\n")

dataFile.write("Angles\n\n")

count = 1
for i in range(numMoly):
    m = i*4+1

    dataFile.write(str(count)+" 1 " + str(m) + " " +str(m+1) + " " + str(m+2) + "\n")
    count += 1
    dataFile.write(str(count)+" 2 " + str(m+1) + " " +str(m+2) + " " + str(m+3) + "\n")
    count += 1
    dataFile.write(str(count)+" 3 " + str(m) + " " +str(m+2) + " " + str(m+3) + "\n")
    count += 1

dataFile.write("\n")

dataFile.write("Ellipsoids\n\n")

for i in range(numMoly):
    m = i*4+1
    dataFile.write(str(m)+" 4.834 0.096 4.834 1 0 0 0\n")

dataFile.write("\n")

dataFile.write("CoreIDs\n\n")


count = 1
for i in range(numMoly):
    dataFile.write(str(count)+" "+str(i+1)+"\n")
    count += 1
    dataFile.write(str(count)+" "+str(i+1)+"\n")
    count += 1
    dataFile.write(str(count)+" "+str(0)+"\n")
    count += 1
    dataFile.write(str(count)+" "+str(0)+"\n")
    count += 1
