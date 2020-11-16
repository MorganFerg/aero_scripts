import numpy as np
import math
from scipy import linalg
from scipy.linalg import toeplitz
from scipy import sparse
from scipy.sparse.linalg import spsolve


af_name = "combined.dat"
avg_density = 24.8 #density of pink foam
chord = 0.275 #chord of the foil, dummy
spar_rad = 0.016#radius of the spar
hspan = 0.762#wingspan of half of the wing
span = 1.524
spar_mass = 0.05#mass of the spar

N = 1
g = 9.81
z = np.linspace(0,hspan,100)
E = 2.28*10**11#young's modulus
q_inf = 551

def L(z,q_inf,span,chord):
    #insert lift equations
    
    L = 4*q_inf*span*(0.075*math.pi*chord)/(math.pi*chord + span*math.sqrt(2))*math.sin(math.acos(2*z/span))
    return L

def load(af_name): # loads airfoil datapoints into the workspace
    f = open(af_name, "r")
    L = f.readlines()
    size = len(L)
    x = np.ones(size)
    y = np.ones(size)
    #obtain all points
    for i in range(size):
        row = str.split(L[i])
        x[i] = chord*float(row[0])
        y[i] = chord*float(row[1])
    return x,y
def moi(af_name, avg_density, chord, spar_rad, spar_mass, hspan):
    x,y = load(af_name)
    #integrate?
    moi = 0 #initialize moi as 0, assuming equal mass distribution
    A = 0
    for j in range(len(x) - 1):
        l = abs(x[j] -x[j+1])
        h = abs(0.5*(y[j] -y[j+1]))
        dA = h*l
        A += dA
        moi += avg_density*hspan*dA*(x[j] + l/2)
    A -= math.pi*spar_rad**2
    d = chord*0.25 #location of center of spar
    moi = moi - (0.5*avg_density*hspan*math.pi*spar_rad**4) - (avg_density*hspan*math.pi*spar_rad**2*d**2)
    moi = moi +  spar_mass*spar_rad**2 + spar_mass
    dprime = ((0.0904*chord/2)**2+(0.25*chord)**2)
    moi = (A*hspan +spar_mass)*dprime**2
    return moi


def m(af_name, avg_density, chord, spar_rad, spar_mass, span):
    x,y = load(af_name)
    A = 0
    m = np.ones(100)
    for j in range(len(x) - 1):
        l = abs(x[j] -x[j+1])
        h = abs(0.5*(y[j] -y[j+1]))
        A += h*l
    A -= spar_rad**2 
    fm = A*span/100 +  spar_mass/100
    M = fm*m
    return M

#method to add masses to locations on wing
def addmass(m,z,zloc, mass):
    i = abs(np.asarray(z[0:len(z) -1]) - zloc).argmin()
    m[i] += mass

#generates the Q(z) function that governs the bending equations    
def q(z, q_inf, af_name, avg_density, chord, spar_rad, spar_mass, span,hpsan):
    q= np.zeros(100)
    m1 =m(af_name, avg_density,chord,spar_rad,spar_mass,hspan)
    for i in range(100):
        zp = z[i]
        mp = m1[i]
        q[i] = L(zp,q_inf, span,chord) - N*g*mp

    return q
#create values for Q and I
Q = q(z, q_inf, af_name, avg_density, chord, spar_rad, spar_mass, span, hspan)
I = moi(af_name, avg_density, chord, spar_rad, spar_mass, hspan)
def solve(Q,z,E,I):

    S = np.zeros(100)
    M = np.zeros(100)
    print("######## Solving for S")
    for i in range(98):
        n = 98-i
        S[n] = S[n+1] - (Q[n+1]+Q[n])/2*(z[n+1] - z[n])
    print("######## Solving for M")
    for i in range(98):
        n = 98-i
        M[n] = M[n+1] - (S[n+1]+S[n])/2*(z[n+1] - z[n])

    theta = np.zeros(100)
    omega = np.zeros(100)
    print("######## Solving for theta")
    for i in range(99):
        i = i+1
        theta[i] = theta[i-1] + 0.5*(M[i] + M[i-1])/(E*I)*(z[i] - z[i-1])
    print("######## Solving for w")
    for i in range(99):
        i = i+1
        omega[i] = omega[i-1] + 0.5*(theta[i] + theta[i-1])*(z[i] - z[i-1])

    return S, M, theta, omega
S, M , O, W = solve(Q,z,E,I)

def write_vec(Vector, file, name):
    file.write(name)
    file.write(",")
    for val in Vector:
      file.write(str(val))
      file.write(",")
    file.writelines("\n")
outfile = "bar_values.csv"
file = open(outfile, "w+")
write_vec(S, file, "S")
write_vec(M, file, "M")
write_vec(O, file, "θ")
write_vec(W, file, "ω")