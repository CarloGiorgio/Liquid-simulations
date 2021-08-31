import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def take_thermo(s):
    data=np.loadtxt(s)
    data_first=[]
    with open(s,"r") as f:
        line=f.readline()
        for t in line.split():
            for j in t.split(":"):
                try :
                    data_first.append(float(j))
                except ValueError:
                    pass
    return [ data[-1,5], data_first[3],data_first[1]]  
end_file=sys.argv[1]+'.txt'
os.chdir("data/NPT")
files=[f for f in os.listdir() if f.endswith(end_file) and f.startswith('NPT')]
data_all=[]
for f in files:
    data_all.append(take_thermo(f))
data_all=np.asarray(data_all)
order=np.argsort(data_all[:,-1])
for i in range(3):
    data_all[:,i]=data_all[order,i]
os.chdir("../EoS")
np.savetxt('eos_NPT_'+end_file,data_all)
plt.plot(data_all[:,0],data_all[:,2],'o',label="data")
if sys.argv[2]=="True":
    rho,p=np.loadtxt('eos_true_1.4.txt',unpack=True)
    plt.plot(rho,p,label='true')
plt.xlabel(r'$\rho^*$')
plt.ylabel(r'$P^*$')
plt.legend()
plt.show()