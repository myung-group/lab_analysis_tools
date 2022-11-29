import numpy as np
from ase.io.trajectory import Trajectory
from ase.io import read,write
from ase.build import sort
from ase import Atoms
from ase.constraints import FixBondLengths
from ase.calculators.tip3p import TIP3P, rOH, angleHOH
from ase.md import Langevin
import ase.units as units
from ase.io.trajectory import Trajectory
import numpy as np
from ase.io import write
from ase.build import sort
from ase.visualize import view
from ase import Atoms
import secrets
from ase.optimize import BFGS
from ase.calculators.emt import EMT

def return_frame(inum=0,filename='md.traj',format='vasp'):
    '''
    Parse the structure at inum time step
    '''
    a=Trajectory(filename)
    if format == 'vasp':
        write(f'POSCAR_t{inum}',sort(a[inum]))
    else:
        write(f'init_t{inum}.xyz',sort(a[inum]),format='extxyz')


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

def aqs_gen(pos_s=[[0.0,0.0,0.0]],mols='Cl',format='vasp',input=None):
    '''
    Generate aqueous solution input file
    mols: solvent molecule
    pos_s: position info of solvent molecule
    '''

    # Set up water box at 20 deg C density
    x = angleHOH * np.pi / 180 / 2
    pos = [[0, 0, 0],
           [0, rOH * np.cos(x), rOH * np.sin(x)],
           [0, rOH * np.cos(x), -rOH * np.sin(x)]]
    atoms = Atoms('OH2', positions=pos)
    
    vol = ((18.01528 / 6.022140857e23) / (0.9982 / 1e24))**(1 / 3.)
    atoms.set_cell((vol, vol, vol))
    atoms.center()
    
    atoms = atoms.repeat((3, 3, 3))
    atoms.set_pbc(True)
    del atoms[:3]
    
    #solvent molecule
    if input is not None:
        solvent=read(input)
    else:
        solvent = Atoms(mols, positions=pos_s)
    vol = ((18.01528 / 6.022140857e23) / (0.9982 / 1e24))**(1 / 3.)
    solvent.set_cell((vol, vol, vol))
    solvent.center()
    print(solvent.positions)
    
    atoms+=solvent
    
    atoms=sort(atoms)
    
    if format == 'vasp':
        write('POSCAR',atoms)
    else: 
        write('init.xyz',atoms,format='extxyz')

def random_str_generator(system=6*['Au'], d=3.0, print_xyz='init.xyz'):
    '''
    Random structure generator
     - Gas phase
     - Surface?
    '''
    
    emt_element=['H', 'C', 'N', 'O', 'Al', 'Ni', 'Cu', 'Pd', 'Ag', 'Pt', 'Au']

    for i,atom in enumerate(system): 
        if i>0:
            disp=np.random.rand(1,3)
            disp*=d/(np.linalg.norm(disp))
            iadd=secrets.randbelow(len(atoms))
            pos=atoms.get_positions()[iadd] + disp
            atoms+=Atoms(f'{atom}',positions=pos)
        else:
            atoms=Atoms(f'{atom}',positions=[[0,0,0]])
            
    if [*set(system)][0] in emt_element:
        atoms.calc=EMT()
        dyn = BFGS(atoms)
        dyn.run(fmax=0.01)
    
    write(f'{print_xyz}',atoms)
    
    return atoms
