#Project trials

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D



#import dataset
pdbset= open('4r52.pdb')

#loop to place all of necessary information from pdb into an array
array1=[]

for first in pdbset.readlines():
    if first[0:4]=='ATOM':
        line=first.split()
        array1.append(line)
#creation of dataframe
df=pd.DataFrame(data=array1)
#simplification of dataframe and organization
df2=df.iloc[:,[2,3,4,5,10]]
df2=df2.rename(columns={2:'Atom', 3:'Residue', 4:'Chain', 5:'Noresidue', 
                        10:'Bfactor'})
#make B-factors into numbers
df2=df2.apply(pd.to_numeric, errors='ignore', downcast='float')

#pivot table to categorize by Amino Acid number and Chain, creates a new dataframe    
newdf=pd.pivot_table(df2, index='Noresidue', columns='Chain', values='Bfactor',
                     aggfunc=np.mean)

#looping over the dataframe to plot all average functions
for frame in newdf:
    plt.figure()
    plt.xlabel('Amino Acid Number')
    plt.ylabel('Average B-factor')
    plt.title('Chain '+frame)
    plt.scatter(x=newdf.index, y= newdf[frame], c=newdf[frame], cmap='magma')
#main chain B-factors
mainchain=(df2['Atom']=='N') | (df2['Atom']=='C') | (df2['Atom']=='CA') | (df2['Atom']=='O')
df3=df2[mainchain]
maindf=pd.pivot_table(df3, index='Noresidue', columns='Chain', values='Bfactor', 
                      aggfunc=np.mean)

for chain in maindf:
    plt.figure()
    plt.xlabel('Amino Acid Number')
    plt.ylabel('Average Main Chain B-factor')
    plt.title('Chain '+chain)
    plt.scatter(x=maindf.index, y= maindf[chain], c=maindf[chain], cmap='magma')
    

#3D visualization Data preparation
df3d=df.iloc[:,[2,3,4,5,6,7,8,10]]
df3d=df3d.rename(columns={2:'Atom', 3:'Residue', 4:'Chain', 5:'Noresidue',
                          6:'x', 7:'y', 8:'z', 10:'Bfactor'})
dmainchain=(df3d['Atom']=='N') | (df3d['Atom']=='C') | (df3d['Atom']=='CA') | (df3d['Atom']=='O')
df3d=df3d[dmainchain]
df3d=df3d.apply(pd.to_numeric, errors='ignore', downcast='float')
x=df3d['x']
y=df3d['y']
z=df3d['z']
#to plot 3D protien
fig=plt.figure()
ax=fig.gca(projection='3d')
ax.scatter(x, y, z, c=df3d['Bfactor'], cmap='plasma', s=df3d['Bfactor'])
ax.plot(x, y, z, c='blue', alpha=0.3)
#Removes the plotting axis background
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

ax.set_xticks([]) 
ax.set_yticks([]) 
ax.set_zticks([])                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           

#normalize data and plot graphs overall and main-chain
#overall b-factors normalized
normnewdf=(newdf-newdf.min())/(newdf.max()-newdf.min())
for normchain in normnewdf:
    plt.figure()
    plt.xlabel('Amino Acid Number')
    plt.ylabel('Normalized Average B-factor')
    plt.title('Chain '+normchain)
    plt.scatter(x=normnewdf.index, y= normnewdf[normchain], c=normnewdf[normchain], cmap='magma')

#mainchain normalized   
normmainchain=(maindf-maindf.min())/(maindf.max()-maindf.min())

for chainnorm in normmainchain:
    plt.figure()
    plt.xlabel('Amino Acid Number')
    plt.ylabel('Main-Chain Normalized Average B-factor')
    plt.title('Chain '+chainnorm)
    plt.scatter(x=normmainchain.index, y= normmainchain[chainnorm], c=normmainchain[chainnorm], cmap='magma')
    
    
    



    



 



        
        