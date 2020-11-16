import numpy as np
import math as math
import matplotlib.pyplot as plt

fn = open("circle.dat", "w+")
c_res = 24 #number of points along circle
xpc=[1]*c_res
ypc=[1]*c_res
for j in range(c_res):
    xpc[j] = math.cos(2*math.pi/(c_res-1)*j + math.pi/2) #start from lower point of circle
    ypc[j] = math.sin(2*math.pi/(c_res-1)*j + math.pi/2)
for l in range(len(xpc)):
    fn.writelines("  " + "{:.5f}".format((xpc[l])) + "  " + "{:.5f}".format((ypc[l])) +"\n")

