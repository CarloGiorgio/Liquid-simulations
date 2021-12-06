import numpy as np
#from scipy.integrate import trapezoid as tp
import matplotlib.pyplot as plt
#import sys
import os
plt.style.use('classic')
plt.rcParams['figure.facecolor']='white'

def potential_integrate(x):
    return 4*(1/(x**10)-1./(x**4))

def force_integrate(x):
    return 48*(1./(x**10)-0.5/(x**4))

def g_r(s,ter,_ax,c,k,_plot=True):
    data=np.loadtxt(s)
    with open(s,'r') as f:
        l=f.readline()
        l=l.split()
        rho=float(l[1].split(':')[-1])
        T=float(l[3].split(':')[-1])
    print(T)
    arg_e=np.array([potential_integrate(data[i,0]) \
        for i in range(data.shape[0])])*2*np.pi*rho

    arg_p=np.array([force_integrate(data[i,0]) \
        for i in range(data.shape[0])])*2*np.pi*rho*rho/3

    ter.append([rho,np.trapz(arg_e*data[:,1],data[:,0]),
    (rho*T)+np.trapz(arg_p*data[:,1],data[:,0])])



    print(f' Density: {ter[-1][0]} Energy: {ter[-1][1]} \
        Pressure {ter[-1][2]}')
    
    if (rho<0.31 or rho>0.05) and _plot:
        _ax.plot(data[:,0],7-k+data[:,1],color=c,linewidth=2.5)
        _ax.text(data[-1,0]*0.8,(7-k+data[-1,1])+0.4,f'$\\rho:${rho:2.2}')
        _ax.axhline(7-k+1,color='green',linestyle='--',linewidth=0.5)


parent=os.getcwd()
dirr=os.path.join(os.getcwd(),'data','gr')
os.chdir(dirr)
ter=[]
fig=plt.figure(figsize=(6,8))
ax=fig.add_subplot(111)
color=plt.get_cmap('jet')
color= [color(i) for i in np.linspace(0, 1,len(os.listdir()))]
i=0
for f, c in zip(os.listdir(),color):
    g_r(f,ter,ax,c,i,True)
    i+=1
ter=np.array(ter)
ar=np.argsort(ter[:,0])
plt.show()
i=0
files=os.listdir()
files=[files[i] for i in ar]

#
os.chdir(parent)
ter=np.array(ter)
ar=np.argsort(ter[:,0])
ter=ter[ar,:]
np.savetxt('Plots/Data_ultimate/gr_res.txt',ter)
#

os.chdir(dirr)
ter=[]
fig,ax=plt.subplots(2,2,figsize=(10,10),sharey=True)
k=0
for f, c in zip(files,color):
    k=int((i//7)%2)
    j=int((i//7)//2)
    g_r(f,ter,ax[k,j],c,i%7,True)
    i+=1
ax[0,0].set_ylabel('$g(r)$',fontsize=18)
ax[1,0].set_ylabel('$g(r)$',fontsize=18)
ax[1,0].set_xlabel('$r$',fontsize=18)
ax[1,1].set_xlabel('$r$',fontsize=18)

"""for i in range(2):
    for j in range(2):
        ax[i,j].legend()
        """
plt.show()

os.chdir(parent)
ter=np.array(ter)
ar=np.argsort(ter[:,0])
ter=ter[ar,:]
np.savetxt('Plots/Data_ultimate/gr_res.txt',ter)
