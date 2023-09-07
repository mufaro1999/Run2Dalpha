import numpy as np
import matplotlib.pyplot as plt
def pees(charges, locs, r2, q2, d): # charges and locs are for sphere 1 and we are finding locations of charges and their values on sphere 2 with radius r2 and charge q2
    newcharges = []
    newlocs = []
    for i in range(len(charges)):   
            p = r2*r2/(d-locs[i])
            newcharges.append(-charges[i]*p/r2)
            newlocs.append(p)
    newcharges.append(-sum(newcharges) + q2)
    newlocs.append(0)
    return(newcharges, newlocs)

if __name__ == "__main__":
    charges=[1]
    locs = [0]
    isit1 =True
    count = 0
    force = []
    d = 2.5
    while(True):
        if(isit1):
            newcharges, newlocs = pees(charges, locs, 1, 1,d)
        if(not isit1):
            newcharges, newlocs = pees(charges, locs, 1, 1,d)
        #now make for loop to calculate force between charges, newcharges and locs newlocs
        emforce = 0
        for i in range(len(newcharges)):
             for j in range(len(charges)):
                  emforce = emforce + (charges[j]*newcharges[i])/((d-newlocs[i]-locs[j])**2)
        force.append(emforce)
        locs = newlocs
        charges = newcharges   
        isit1 = not isit1 
        count = count + 1
        if (count == 10):
             break
    #print(charges)
    #print(locs)
    print(force)
    plt.plot(force)
    plt.xlabel('iteration', fontsize=14)
    plt.ylabel(r'Force($\Frac{1}{4 \pi \Epsilon}}$)',fontsize=14)
    plt.show()
    plt.save('convergence.png')








'''
while(True): 
    if(isit1): 
        charges1 = []
        for i in range(count):
            pees1.append(r1*r1/(d-pees2[i]))
            charges1.append(-charges2[i]*pees1[i]/r1)
        charges1.append(q1+sum(-charges1))
    if(not isit1):
        charges2 = []
        for i in range(count):
            pees2.append(r2*r2/(d-pees1[i]))
            charges2.append(-charges1[i]*pees2[i]/r1)
        charges1.append(q1+sum(-charges1))

'''
   


 