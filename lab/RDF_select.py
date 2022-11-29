import matplotlib.pyplot as plt
import scipy as sp
from ase.io import read, write
from ase.io.trajectory import Trajectory
import numpy as np
import pandas as pd
import math
from os import remove

def distance(a, b, cell_length):
    dx = abs(a[0] - b[0])
    x = min(dx, abs(cell_length - dx))
    
    dy = abs(a[1] - b[1])
    y = min(dy, abs(cell_length - dy))
    
    dz = abs(a[2] - b[2])
    z = min(dz, abs(cell_length - dz))   
    
    return np.sqrt(x**2 + y**2 + z**2)

def atom_name(w):
    atoms = {
    'H':1, 'He':2, 'Li':3, 'Be':4, 'B':5, 'C':6, 'N':7, 'O':8, 'F':9, 'Ne':10,
    'Na':11, 'Mg':12, 'Al':13, 'Si':14, 'P':15, 'S':16, 'Cl':17, 'Ar':18, 'K':19, 'Ca':20,
    'Sc':21, 'Ti':22, 'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28,'Cu':29, 'Zn':30,
    'Ga':31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36, 'Rb':37, 'Sr':38, 'Y':39, 'Zr':40,
    'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48, 'In':49, 'Sn':50,
    'Sb':51, 'Te':52, 'I':53, 'Xe':54, 'Cs':55, 'Ba':56, 'La':57, 'Ce':58, 'Pr':59, 'Nd':60,
    'Pm':61, 'Sm':62, 'Eu':63, 'Gd':64, 'Tb':65, 'Dy':66, 'Ho':67, 'Er':68, 'Tm':69, 'Yb':70,
    'Lu':71, 'Hf':72, 'Ta':73, 'W':74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80,
    'Ti':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86, 'Fr':87, 'Ra':88, 'Ac':89, 'Th':90,
    'Pa':91, 'U':92, 'Np':93, 'Pu':94, 'Am':95, 'Cm':96, 'Bk':97, 'Cf':98, 'Es':99, 'Fm':100,
    'Md':101, 'No':102, 'Lr':103, 'Rf':104, 'Db':105, 'Sg':106, 'Bh':107, 'Hs':108, 'Mt':109, 'Ds':110,
    'Rg':111, 'Cn':112, 'Nh':113, 'Fl':114, 'Mc':115, 'Lv':116, 'Ts':117, 'Og':118
    }
    
    return atoms[w]

def Trajectory_to_XDATCAR(filename, steps):    # filename = traj 파일
    a=Trajectory(filename)
    XDAT=a[-steps:]
    write('XDATCAR', XDAT)
    return a, XDAT

def extracter(filename):    # filename = XDATCAR
    ### cell_parameter ### (1.0)
    cell_parameter = float(np.loadtxt(filename, skiprows= 1, max_rows=1))
    
    ### cell_length(lattice_constant) ### (15.60... 0 0)
    extract_cell_length = np.loadtxt(filename, skiprows=2, max_rows=3)
    A=extract_cell_length[0][0]
    B=extract_cell_length[1][1]
    C=extract_cell_length[2][2]
    cell_length = A
    
    ### num_of_atoms ### (각 원자 개수)
    atoms = np.loadtxt('XDATCAR' ,skiprows= 6, max_rows=1)    # filename = XDATCAR
    total_atoms=np.sum(atoms)
    num_of_Atoms = int(total_atoms)
    
    ### couple ###
    elements = np.loadtxt('XDATCAR' ,skiprows= 5, max_rows=1,dtype= 'str')
    atoms = np.loadtxt('XDATCAR' ,skiprows= 6, max_rows=1)
    couple = dict(zip(elements, atoms))
    
    
    return num_of_Atoms, cell_parameter, cell_length, A, B, C, couple

def make_rdf_fig(a, b):
    plt.plot(a, b)
    #plt.title('radial distribution function of liquid water')
    #plt.xlim([0, 8])    <- cutoff로 대체 
    #plt.ylim([0, 5])
    plt.xlabel('r')
    plt.ylabel('g(r)')
    #plt.show

def draw_RDF_select(filename, steps, cutoff, Resolution, sys_name, atom1, atom2,
                    save_fig=False, del_XDATCAR=True) :
    
    #함수 불러오기
    a, XDAT = Trajectory_to_XDATCAR(filename, steps)    
    num_of_Atoms, cell_parameter, cell_length, A, B, C, couple = extracter('XDATCAR')

    at1 = atom_name(atom1)
    at2 = atom_name(atom2)

    atom_position1 = []
    atom_position2 = []
    atom_position3 = [] 

    for i in range(len(XDAT)):
        for j in range(len(a[0].get_atomic_numbers())):
            if atom_name(atom1) == a[i].get_atomic_numbers()[j] :
                atom_position1.append(a[i].get_positions()[j])
            
            if atom_name(atom2) == a[i].get_atomic_numbers()[j] :
                atom_position2.append(a[i].get_positions()[j])
                       
        for k in range (2+i*int(couple[atom2]), (i+1)*int(couple[atom2]), 3):  # 이중결합 산소(3번째 'O') 추출
            atom_position3.append(atom_position2[k])

                
    atom_position = []
    aa=[]

    for i in range(len(XDAT)):
        # atom_position1    
        for j in range(i*int(len(atom_position1)/len(XDAT)), (i+1)*int(len(atom_position1)/len(XDAT))) :   
            aa.append(atom_position1[j])
        # atom_position3  
        for k in range(i*int(len(atom_position3)/len(XDAT)), (i+1)*int(len(atom_position3)/len(XDAT))) :
            aa.append(atom_position3[k])
    
    atom_position.append(aa)            
    
    # (atom1+atom2)개씩 분리
    atom_position_separation = []
    for i in range(len(XDAT)):
        b=[]
        for j in range (i*int(couple[atom1]+(couple[atom2]/3)), (i+1)*int(couple[atom1]+(couple[atom2]/3))):
            aa = atom_position[0][j]%cell_length
            b.append(aa)
        atom_position_separation.append(b)
        
    
    select_r_list = []
    for k in range(len(XDAT)):
        for i in range(0, (int(couple[atom1])+int(couple[atom2]/3)-1)) :
            for j in range(i+1, (int(couple[atom1])+int(couple[atom2]/3))) :
                r=distance(atom_position_separation[k][i], atom_position_separation[k][j], cell_length)
                if r <= cutoff :
                    select_r_list.append(r)
        
    R=np.array(select_r_list)
    
    bn=np.linspace(0, cutoff, cutoff*Resolution+1)
    nr=pd.cut(R, bn)    # Categorical형
    ns=pd.Series.value_counts(nr)    # 각 거리별 몇개 있는지, Series형
    ni=ns.sort_index()    # index 오름차순 정렬
    na=np.array(ni)    # 거리별 개수
    
    # 그래프 y축, g(r)
    RDF=[
        (A*B*C)*(na[i])/
        (2*math.pi*(int(couple[atom1])+int((couple[atom2]/3))**2)*len(XDAT)*((bn[1]-bn[0])**3)*((i+0.5)**2))
        for i in range(len(na))
        ]
    
    # 그래프 x축, r
    bx=np.linspace(1/Resolution, cutoff, cutoff*Resolution) 
    
    # 그래프 그리기
    make_rdf_fig(bx, RDF)
    
    # 사진 저장 유/무
    if save_fig:
        plt.savefig(f'rdf_{sys_name}.png')
    
    # XDATCAR 삭제 유/무
    if del_XDATCAR:
        remove('XDATCAR')

draw_RDF_select('md_310.traj', 1000, 7, 10, 'ec310_None_select', 'Li', 'O',
                del_XDATCAR=False)


