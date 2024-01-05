# This module is a part from XMECP program.
# This module is used for adding external electrostatic field in QM/MM calculations.
# External electrostatic field tracks one bond or plane when applied to optimization jobs.

# 2023-03-14
# xujiawei@fjirsm.ac.cn

import math,numpy
from parameters import *

def plane_normvec(A,B,C,x1,y1,z1): # Z=AX+BY+C
    # Project (x1,y1,z1) onto the plane.
    z0=(z1+(A*x1+B*y1+C))/(1+1/(A**2+B**2))
    x0=x1+A*(z1-z0)
    y0=y1+B*(z1-z0)
    # Normalize.
    normvec=numpy.array([x1,y1,z1])-numpy.array([x0,y0,z0])
    x,y,z=normvec/numpy.linalg.norm(normvec)
    return x,y,z # Norm vector of Z=AX+BY+C passing (x1,y1,z1); Direction: plane->(x1,y1,z1)

def fitplane(points,refpoint=[]):
    print(f'Total {len(points)} atoms are selected to define the plane.\n')

    if len(points)==3 and refpoint==[]:
        sys.exit('Error: You must define a reference atom to mark the direction of ESF when only three atoms given to define the plane!')
    elif len(points)>3 and refpoint==[]:
        print('Reference atom not defined.')
        print('The atom fastest from the fitted plane will be set as reference.')

    tmpA,tmpB=[],[]
    for i in points:
        tmpA.append([i[0],i[1],1])
        tmpB.append(i[2])
    
    A=numpy.matrix(tmpA)
    B=numpy.matrix(tmpB).T

    coeffA,coeffB,coeffC=numpy.matrix.tolist(((A.T*A).I*A.T*B).T)[0]

    print('The plane fitted by least square method is:')
    print('    A = '+f'{coeffA:.8f}'.rjust(16))
    print('    B = '+f'{coeffB:.8f}'.rjust(16))
    print('    C = '+f'{coeffC:.8f}'.rjust(16)+'  (Z=AX+BY+C)\n')

    rmax,MPP=0.,0.
    xmax,ymax,zmax=0.,0.,0.
    for i in range(len(points)):
        x,y,z=points[i]
        r=point2plane(x,y,z,coeffA,coeffB,coeffC)
        MPP+=r**2
        if r>rmax:
            rmax=r
            xmax,ymax,zmax=x,y,z
    MPP=math.sqrt(MPP/len(points))

    print(f'The atom fastest from the fitted plane locates at:')
    print('    x: '+f'{xmax:.8f}'.rjust(16))
    print('    y: '+f'{ymax:.8f}'.rjust(16))
    print('    z: '+f'{zmax:.8f}'.rjust(16))

    if refpoint==[]:
        print('    This is set as reference atom to mark the ESF direction.')
    else:
        print('User-defined reference atom to mark the ESF direction locates at:')
        print('    x: '+f'{refpoint[0]:.8f}'.rjust(16))
        print('    y: '+f'{refpoint[1]:.8f}'.rjust(16))
        print('    z: '+f'{refpoint[2]:.8f}'.rjust(16))
        xmax,ymax,zmax=refpoint

    print(f'The MPP indice of selected atom group is: {MPP:.8f}.')
    print('    Reference: Journal of Molecular Modeling, 27, 263 (2021)')

    # Now calculate the norm vector passing the reference point.
    vecx,vecy,vecz=plane_normvec(coeffA,coeffB,coeffC,xmax,ymax,zmax)

    print('The normalized ESF vector (orthogonal to Z=AX+BY+C):')
    print('    x: '+f'{vecx:.8f}'.rjust(16))
    print('    y: '+f'{vecy:.8f}'.rjust(16))
    print('    z: '+f'{vecz:.8f}'.rjust(16)+'\n')

    return coeffA,coeffB,coeffC,vecx,vecy,vecz

def median_plane(points):
    x1,y1,z1=points[0]
    x2,y2,z2=points[1]

    vecx,vecy,vecz=x2-x1,y2-y1,z2-z1

    coeffA=-(vecx)/(vecz)
    coeffB=-(vecy)/(vecz)
    coeffC=((x2**2+y2**2+z2**2)-(x1**2+y1**2+z1**2))/(2*vecz)

    print('The median plane of selected two atoms is:')
    print('    A = '+f'{coeffA:.8f}'.rjust(16))
    print('    B = '+f'{coeffB:.8f}'.rjust(16))
    print('    C = '+f'{coeffC:.8f}'.rjust(16)+'  (Z=AX+BY+C)')

    E=numpy.array([vecx,vecy,vecz])
    vecx,vecy,vecz=E/numpy.linalg.norm(E)

    print('The normalized ESF vector (orthogonal to Z=AX+BY+C):')
    print('    x: '+f'{vecx:.8f}'.rjust(16))
    print('    y: '+f'{vecy:.8f}'.rjust(16))
    print('    z: '+f'{vecz:.8f}'.rjust(16)+'\n')

    return coeffA,coeffB,coeffC,vecx,vecy,vecz

def point2plane(x,y,z,A,B,C): # Plane equation: Z=AX+BY+C.
    return abs(A*x+B*y-z+C)/math.sqrt(A**2+B**2+1)

# Format of .chg file:
#     Line =1: Number of charges.
#     Line >1: [Serial or element] [X] [Y] [Z] [charge]
def readchg(chgfile):
    chg,nchg=[],0
    with open(chgfile,'r') as read:
        for line in read.readlines():
            l=line.split()
            if len(chg)==nchg:
                break
            if len(l)==1:
                nchg=int(l[0])
                continue
            elif len(l)>1:
                chg.append(float(l[-1]))
    return chg

# You can also use .prmtop file as charge input:
def readchg_from_prmtop(prmtop):
    chg=[]
    with open(prmtop,'r') as read:
        flag=False
        for line in read.readlines():
            l=line.split()
            if '%FLAG CHARGE' in line:
                flag=True
                continue
            elif '%FLAG ATOMIC_NUMBER' in line:
                break
            if flag and len(l)!=0:
                import geom
                if geom.is_number(l[0]):
                    for i in l:
                        chg.append(float(i))
    return chg

# Define external electrostatic field:
#     1. Parallel to specific line: [intensity(au),[A1,A2]]
#        The field vector will be parallel to A1->A2.
#     2. Orthogonal to specific plane: [intensity(au),[A1,A2,A3,...],(A0)]
#        The field vector will be parallel to norm vector of fitted plane of [A1,A2,A3,...] passing A0.
#        If A0 not defined, when number of [A1,A2,A3,...] larger than 3, the atom fastest from the fitted plane will be A0.
def electrostatic_field(X,chg,fielddef,calceandg=True):
    print('\n================ >>> External Electrostatic Field <<< ================')

    # Get coordinates.
    fitps,refp=[],[]
    for i in fielddef[1]:
        fitps.append([X[0+3*i],X[1+3*i],X[2+3*i]])
    field=[fielddef[0],fitps]

    if len(fielddef)==3:
        i=fielddef[2]
        refp=[X[0+3*i],X[0+3*i],X[0+3*i]]
        field.append(refp)

    A,B,C=0.,0.,0.
    intensity=field[0]
    esf_type=None
    if len(field[1])==2:
        esf_type=1
        A,B,C,Ex,Ey,Ez=median_plane(field[1])
    elif len(field[1])>2:
        esf_type=2
        if len(field)==2:
            if len(field[1])==3:
                sys.exit('You must give direction when only three points used to define the plane!')
            else:
                A,B,C,Ex,Ey,Ez=fitplane(field[1],[])
        else:
            A,B,C,Ex,Ey,Ez=fitplane(field[1],field[2])
    else:
        sys.exit('Wrong field definition!')

    if calceandg==False:
        scale=math.sqrt(intensity)
        return scale*Ex,scale*Ey,scale*Ez

    esf_bq=[]
    dim=int(len(X)/3)
    global natms
    print(f'  Number of atomic charges: {len(chg)}')
    if len(chg)==dim:
        print('  Same as number of atoms in active region.')
        print('  Charges recognized as MM charges for active region.\n')
        esf_bq.extend(chg)
    elif 3*len(chg)==natms:
        print('  Same as number of atoms of the whole QM/MM system.')
        print('  Charges recognized as MM charges for the whole QM/MM system.\n')
        for i in geom.act_region():
            esf_bq.append(chg[i-1])
    else:
        sys.exit('Wrong with number of input charges. Please check!')

    # Calculate ESF additional energy and gradient.
    items=[]
    for i in range(dim):
        items.append([             intensity,     # ESF intensity
                      numpy.array([Ex,Ey,Ez]),    # ESF vector
                                  [A,B,C],        # ESF zero plane
                                  [X[0+3*i],   # Charge coordinates
                                   X[1+3*i],
                                   X[2+3*i]],
                                   esf_bq[i]])    # Charge

    global nproc
    import multiprocessing
    with multiprocessing.Pool(nproc) as CalcESF:
        ESF_results=CalcESF.map(ESF_eandg,items)

    ESF_ene,ESF_grad=0.,[]
    for i in ESF_results:
        ESF_ene+=i[0]
        ESF_grad.extend(i[1].tolist())
    ESF_grad=numpy.array(ESF_grad)

    print('  Additional electrostatic energy:  '+f'{ESF_ene:.8f}'.rjust(16)+' a.u.')
    print('      Note: XMECP does not add this energy into total energy.')
    print('            If you want to consider it, please add it by yourself.\n')

    print('  Additional electrostatic gradient:'+f'{numpy.linalg.norm(ESF_grad):.8f}'.rjust(16)+' a.u.\n')

    print('====== >>> Successfully Added External Electrostatic Field! <<< ======\n')
    return ESF_ene,ESF_grad

def ESF_eandg(item):
    # Calculate ESF energy of one point charge.
    A,B,C=item[2] # ESF zero plane defined as Z=AX+BY+C.
    x,y,z=item[3] # Position of point charge.
    # Distance to zero plane.
    global bohr
    r=point2plane(x,y,z,A,B,C)/bohr
    # Compare direction.
    normvec=item[1]
    pchgvec=numpy.array(plane_normvec(A,B,C,x,y,z))
    if numpy.dot(normvec,pchgvec)<0:
        r*=-1
    # ESF energy.
    charge=item[4]
    intensity=item[0]
    ESF_ene=-intensity*charge*r
    # ESF gradient.
    ESF_grad=intensity*charge*normvec
    return ESF_ene,ESF_grad # in a.u.
