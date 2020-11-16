e64 = open("e64.dat","r")
sd = open("sd7032.dat","r")
f = open("combined.dat", "w+")

L64 = e64.readlines()
Lsd = sd.readlines()

for i in range(61):

    q1 = str.split(L64[i])
    q2 = str.split(Lsd[i])
    v11 = float(q1[0])
    v12 = float(q1[1])
    v21 = float(q2[0])
    v22 = float(q2[1])
    co1 = (v11 + v21)/2
    co2 = (v12 + v22)/2
    cr1 = ("{:.5f}".format(co1))
    cr2 = ("{:.5f}".format(co2))

    
    f.writelines("  " + cr1 + "  " + cr2 +"\n")