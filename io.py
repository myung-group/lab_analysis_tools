import numpy as np
from ase.io.trajectory import Trajectory
from ase.io import read,write
from ase.build import sort

def return_last_frame(filename='md.traj',format='vasp'):
    '''
    Parse the last time step structure in POSCAR format
    '''
    a=Trajectory(filename)
    if format == 'vasp':
        write('POSCAR_last',sort(a[-1]))
    else:
        write('init_last.xyz',sort(a[-1]),format='extxyz')

def return_min_structure(filename='structures.traj'):
    '''
    Output the minimum structure in POSCAR format
    '''

    #read global optimization traj
    a=Trajectory(filename)
    ene=[i.get_potential_energy() for i in a]
    minarg=a[np.argmin(ene)]

    write('POSCAR_cand',minarg)
