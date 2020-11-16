import numpy as np
import math as math
import matplotlib.pyplot as plt
f = open("0012.dat", "r")
fn = open("spar50.dat", "w+")

L = f.readlines()
s = len(L) #number of coords

xpoints = [1]*s
ypoints = [1]*s

for i in range(s):
    q = str.split(L[i])
    xpoints[i] = float(q[0])
    ypoints[i] = float(q[1])

#locate closest x chord to quarter chord on lower surface
half = math.floor(s/2) +1
index = half + abs(np.asarray(xpoints[half:s -1]) - 0.25).argmin() #index lower surf at quarter chord

#define the points for our hole
c = 0.15 #wing chord
r = 0.50*0.0074549/c #radius of our circle
diff = 0.04 #distance from edge of circle to surface of af
x_offset = xpoints[index] #xloc of circle
y_offset =  ypoints[index] + diff #yloc of circle
c_res = 16 #number of points along circle
xpc=[1]*c_res
ypc=[1]*c_res
for j in range(c_res):
    xpc[j] = x_offset - r*math.cos(2*math.pi/(c_res-1)*j + math.pi/2) #start from lower point of circle
    ypc[j] =  - r*math.sin(2*math.pi/(c_res-1)*j + math.pi/2)

#insert points before hole

for k in range(index+1):
    fn.writelines("  " + "{:.5f}".format((xpoints[k])) + "  " + "{:.5f}".format((ypoints[k])) +"\n")

fn.writelines("  " + "{:.5f}".format((xpoints[index])) + "  " + "{:.5f}".format((ypoints[index])) +"\n")
#insert circle
for l in range(len(xpc)):
    fn.writelines("  " + "{:.5f}".format((xpc[l])) + "  " + "{:.5f}".format((ypc[l])) +"\n")

#insert rest of foil
fn.writelines("  " + "{:.5f}".format((xpoints[index])) + "  " + "{:.5f}".format((ypoints[index])) +"\n")
for m in range(s-index):
    fn.writelines("  " + "{:.5f}".format((xpoints[m+index])) + "  " + "{:.5f}".format((ypoints[m+index])) +"\n")
    
