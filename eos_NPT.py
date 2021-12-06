import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import sys

#function for taking thermodynamics data 
#           [rho,T,P]
# last of data file(mean value), first line of txt filefor the other two

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
    return [ data[-1,5], data_first[3],data_first[1],data[-1,6]] 


 
end_file=sys.argv[1]+'.txt'
os.chdir("data/"+sys.argv[2])

files =glob.glob("*"+end_file)
#files=[f for f in os.listdir() if f.endswith(end_file) and f.startswith('NPT')]
data_all=[]
for f in files:
    data_all.append(take_thermo(f))

data_all=np.asarray(data_all)
order=np.argsort(data_all[:,0])
for i in range(4):
    data_all[:,i]=data_all[order,i]

order=(data_all[:,-1]!=data_all[:,-1])
_f=np.amin(data_all[~order,-1]/data_all[~order,0])
print(_f)


data_all[order,-1]=data_all[order,0]*_f
print(data_all[order,-1])


os.chdir("../EoS")
np.savetxt('eos_NPT_'+end_file,data_all)
plt.errorbar(data_all[:,0],data_all[:,2],xerr=data_all[:,-1],marker='o',linestyle=' ',label="data")

if sys.argv[3]=="True":
    rho,p=np.loadtxt('eos_true_1.4.txt',unpack=True)
    plt.plot(rho,p,label='true')

plt.xlabel(r'$\rho^*$')
plt.ylabel(r'$P^*$')
plt.legend()
plt.show()