import matplotlib.pyplot as plt
from ase.io import read, write
from ase.io.trajectory import Trajectory
import numpy as np
from math import pi
import matplotlib.pyplot as plt
from os import remove

def traj_xdat(filename,frame):
    pictures=Trajectory(filename)
    pictures=pictures[-frame:]
    write('XDATCAR',pictures)
    return pictures

def extracter(filename):
    element = np.loadtxt('XDATCAR' ,skiprows= 5, max_rows=1,dtype= 'str')
    atom = np.loadtxt(filename ,skiprows= 6, max_rows=1)
    num_of_Atoms = 0
    name_tag = []
    
    for i in range(len(atom)):
        num_of_Atoms += atom[i]
        name_tag += [element[i] for j in range(int(atom[i]))]    
    
    num_of_Atoms = int(num_of_Atoms)
    cell_len = int(np.loadtxt(filename,skiprows= 1, max_rows=1))
    lattice_constant = np.loadtxt(filename,skiprows= 2, max_rows=1)[0]
    return num_of_Atoms ,cell_len,lattice_constant,name_tag

def distance (atom1, atom2, cell_len):
    dx = abs(atom1[0] - atom2[0])
    x = min(dx, abs(cell_len - dx))
    
    dy = abs(atom1[1] - atom2[1])
    y = min(dy, abs(cell_len - dy))
    
    dz = abs(atom1[2] - atom2[2])
    z = min(dz, abs(cell_len - dz))   
    
    return np.lib.scimath.sqrt(x**2 + y**2 + z**2)

def RDF_drawer(sys_name,filename,frame,bin,r_max,xlim = None,ylim=None,save_pic = False,del_X = True, pair = None):
    
    pictures = traj_xdat(filename,frame)
    num_of_Atoms ,cell_len,lattice_constant,name_tag = extracter('XDATCAR')
    name_tag*=len(pictures)
   
    xyz = np.loadtxt('XDATCAR',skiprows=8, max_rows = num_of_Atoms)
    r_list=[]
    
    for i in range(1,len(pictures)):
        temp_xyz = np.loadtxt('XDATCAR',skiprows=8+((num_of_Atoms + 1)*i), max_rows = num_of_Atoms)
        xyz = np.vstack([xyz,temp_xyz])
      
    for k in range(0,len(pictures)):    
        for i in range (0+k*num_of_Atoms,(num_of_Atoms-1)+k*num_of_Atoms):
            for j in range(i+1,num_of_Atoms+k*num_of_Atoms):
                if pair:
                    if name_tag[i] in pair and name_tag[j] in pair:
                        if  pair[0] != pair[1] and name_tag[i] != name_tag[j]:
                            r_list+=[distance (xyz[i], xyz[j],cell_len)]
                        elif pair[0] == name_tag[i]:
                            r_list+=[distance (xyz[i], xyz[j],cell_len)]
                else:
                    r_list+=[distance (xyz[i], xyz[j],cell_len)]                         
    r_list.sort()
    
    Na=0
    Nb=0
    if pair: 
        for i in range(int(len(name_tag)/len(pictures))):
            if name_tag[i] == pair[0]:
                Na += 1
            if name_tag[i] == pair[1]:
                Nb += 1    
    else:
        Na=num_of_Atoms
        Nb=num_of_Atoms
    
    r_list = np.array(r_list)*lattice_constant      
    bin=bin
    rmax=r_max
    interval=r_list[-1]/bin
    bins = np.arange(0,rmax,interval)
    hist, bins = np.histogram(r_list, bins)
    h = [i for i in hist]
    
    for i in range(len(hist)):
        N_F = (2.0*pi*(interval**3))/lattice_constant**3
        h[i] = h[i]/(N_F*(i-0.5)**2)
        h[i] = h[i]/(Na*Nb*len(pictures))
    
    x = np.arange(interval, rmax, interval)
    y = h
    plt.plot(x, y)
    if xlim == None:
        xlim = [0,rmax]
    plt.xlim(xlim)      
    plt.ylim(ylim)   
    plt.title(f'RDF of {sys_name}')
    plt.xlabel('distance')
    plt.ylabel('g(r)')
    
    if save_pic:
        plt.savefig(f'rdf_{sys_name}.png')
    if del_X:
        remove('XDATCAR')
