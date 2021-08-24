import numpy as np
import matplotlib.pyplot as plt
import os
import sys

data=np.loadtxt(sys.argv[1])
def potential_integrate(x):
    return 4*(1/(x**10)-1./(x**4))

arg=np.asarray([potential_integrate(data[i,0]) for i in range(len(data[:,0]))])*2*np.pi*0.29
print(np.trapz(arg,data[:,1],data[:,0]))
plt.plot(data[:,0],data[:,1])
plt.show()
