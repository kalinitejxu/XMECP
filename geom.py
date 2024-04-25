bohr=0.52917721067
import math,numpy

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False

def all_numbers(numlist):
    for s in numlist:
        if is_number(s)==False:
            return False
    return True

def readxyzs(xyzfile):
    num_atm=0
    counter=0
    xyzs=[]
    atm,singlexyz=[],[]
    xyz=open(xyzfile,'r')
    for line in xyz.readlines():
        l=line.split()
        if len(l)!=0:
            if len(l)==1 and is_number(l[0]):
                num_atm=int(l[0])
                break
    flag,thefirst=False,True
    geom=0
    xyz.close()
    xyz=open(xyzfile,'r')
    for line in xyz.readlines():
        l=line.split()
        if len(l)!=0:
            if flag==False and len(l)==1 and int(l[0])==num_atm:
                flag=True
                geom+=1
                if geom>=2:
                    thefirst=False
                    xyzs.append(singlexyz)
                counter,singlexyz=0,[]
                continue
            if flag:
                flag=False
                continue
            if len(l)==4 and is_number(l[0])==False and is_number(l[1]) and is_number(l[2]) and is_number(l[3]):
                counter+=1
                if thefirst and len(atm)<num_atm:
                    atm.append(l[0])
                elif len(atm)>num_atm or counter>num_atm or (thefirst==False and l[0]!=atm[counter-1]):
                    print('Problems with geomtries! Exit program!')
                    exit()
                for i in range(3):
                    singlexyz.append(float(l[1+i]))
    xyzs.append(singlexyz)
    return atm,xyzs

def writexyz(atm,xyz,filename,comment='xyz'):
    xyzfile=open(filename,'w')
    xyzfile.write(str(len(atm))+'\n'+comment)
    for i in range(len(atm)):
        xyzline='\n'+atm[i]
        for j in range(3):
            xyzline+='    '+str(xyz[j+3*i])
        xyzfile.write(xyzline)
    xyzfile.write('\n')
    xyzfile.close()
    return

def writexyzs(atm,xyzs,filename,comment='xyz'):
    xyzsfile=open(filename,'w')
    for i in range(len(xyzs)):
        xyz=xyzs[i]
        com=comment+str(i)
        xyzsfile.write(str(len(atm))+'\n'+com)
        for i in range(len(atm)):
            xyzline='\n'+atm[i]
            for j in range(3):
                xyzline+='    '+str(xyz[j+3*i])
            xyzsfile.write(xyzline)
        xyzsfile.write('\n')
    xyzsfile.close()
    return

# Special for Q-Chem.
def qchemxyz(atm,xyz,chg=0,mult=1,filename='qchem.geom'):
    xyzfile=open(filename,'w')
    xyzfile.write(f'$molecule\n{chg} {mult}')
    for i in range(len(atm)):
        xyzline='\n'+atm[i]
        for j in range(3):
            xyzline+='    '+str(xyz[j+3*i])
        xyzfile.write(xyzline)
    xyzfile.write('\n$end\n')
    xyzfile.close()
    return

# Special for ChemShell.
def conn_from_prmtop(prmtop):
    flag1,flag2=False,False
    binch,bouth=[],[]
    with open(prmtop,'r') as prm:
        for line in prm.readlines():
            if '%FLAG BONDS_INC_HYDROGEN' in line:
                flag1=True
                continue
            elif '%FLAG BONDS_WITHOUT_HYDROGEN' in line:
                flag2=True
                flag1=False
                continue
            elif '%FLAG ANGLES_INC_HYDROGEN' in line:
                break
            if flag1:
                l=line.split()
                if all_numbers(l):
                    for i in l:
                        binch.append(int(i))
                continue
            elif flag2:
                l=line.split()
                if all_numbers(l):
                    for i in l:
                        bouth.append(int(i))
    conn=[]
    for i in range(int(len(binch)/3)):
        conn.append([1+int(binch[0+3*i]/3),1+int(binch[1+3*i]/3)])
    for i in range(int(len(bouth)/3)):
        conn.append([1+int(bouth[0+3*i]/3),1+int(bouth[1+3*i]/3)])
    return conn # List[[atm1,atm2],etc.]

def chemshc(atm,xyz,prmtop,filename='chemsh.c'):
    global bohr
    conn=conn_from_prmtop(prmtop)
    with open(filename,'w') as out:
        out.write('block = fragment records = 0\n')
        out.write('block = title records = 1\n')
        out.write('Generated by XMECP\n')
        out.write(f'block = coordinates records = {len(atm)}\n')
        for i in range(len(atm)):
            out.write(f'{atm[i].lower()} {xyz[0+3*i]/bohr:.14e} {xyz[1+3*i]/bohr:.14e} {xyz[2+3*i]/bohr:.14e}\n')        
        out.write(f'block = connectivity records = {len(conn)}\n')
        for i in range(len(conn)):
            out.write(f'{conn[i][0]} {conn[i][1]}\n')
    return

def active_region(tag='active_atom_numbers'):
    from parameters import chemsh_template
    act_region=[]
    flag=False
    with open(chemsh_template,'r') as read:
        for line in read.readlines():
            l=line.split()
            if 'set' in line and tag in line:
                flag=True
            if flag:
                for i in l:
                    if is_number(i):
                        act_region.append(int(i))
                    elif all_numbers(i.split('-')):
                        for j in range(int(i.split('-')[0]),1+int(i.split('-')[1])):
                            act_region.append(j)
                if line[-2:]!='\\\n':
                    flag=False
    return act_region

def active_xyz(xyz):
    act_region=active_region()
    actxyz=[]
    for i in act_region:
        for j in range(3):
            actxyz.append(xyz[j+3*i-3])
    return actxyz

def active_gradient(gradient):
    act_region=active_region()
    actgrad=[]
    for i in range(int(len(gradient)/3)):
        if 1+i in act_region:
            for j in range(3):
                actgrad.append(gradient[j+3*i])
    return actgrad # Gradients of atoms not in the active region are deleted.

def upd_active_xyz(Xact):
    from parameters import geom0
    act_region=active_region()
    act,Xupd=0,[]
    for i in range(int(len(geom0)/3)):
        if 1+i in act_region:
            for j in range(3):
                Xupd.append(Xact[j+3*act])
            act+=1
        else:
            for j in range(3):
                Xupd.append(geom0[j+3*i])
    return Xupd

def qm_region(tag='qm_atom_numbers'):
    from parameters import chemsh_template
    qm=[]
    flag=False
    with open(chemsh_template,'r') as read:
        for line in read.readlines():
            l=line.split()
            if 'set' in line and tag in line:
                flag=True
            if flag:
                for i in l:
                    if is_number(i):
                        qm.append(int(i))
                    elif all_numbers(i.split('-')):
                        for j in range(int(i.split('-')[0]),1+int(i.split('-')[1])):
                            qm.append(1+j)
                if line[-2:]!='\\\n':
                    flag=False
    return qm

def qm_xyz(xyz):
    qm=qm_region()
    qmxyz=[]
    for i in qm:
        for j in range(3):
            qmxyz.append(xyz[j+3*i-3])
    return qmxyz

def qm_gradient(gradient):
    qm=qm_region()
    qmgrad=[]
    for i in qm:
        for j in range(3):
            qmgrad.append(gradient[j+3*i-3])
    return qmgrad

def upd_qm_xyz(Xqm):
    from parameters import geom0
    qm=qm_region()
    n,Xupd=0,[]
    for i in range(int(len(geom0)/3)):
        if 1+i in qm:
            for j in range(3):
                Xupd.append(Xqm[j+3*n])
            n+=1
        else:
            for j in range(3):
                Xupd.append(geom0[j+3*i])
    return Xupd

def qm_from_active(XorG):
    XorG=upd_active_xyz(XorG)
    return qm_xyz(XorG)
    # Note: QM region should be included in active region and usually it is this case in QM/MM calculations. If needed, add temperary variable for storage.

# Special for linear insert.
def bond(a,b):
    return numpy.linalg.norm(a-b)

def angle(a,b,c,radunit=False):
    ab=a-b
    cb=c-b
    cos=numpy.dot(ab,cb)/(numpy.linalg.norm(ab)*numpy.linalg.norm(cb))
    abc=numpy.arccos(cos)
    if radunit==False:
        abc=180.0*abc/math.pi
    return abc

def dotpreptoline(a,b,c):
    fa=(a-b)-numpy.dot(a-b,c-b)/numpy.linalg.norm(c-b)*(c-b)
    h=numpy.linalg.norm(fa)
    f=a-fa
    return h,f

def dihedral(a,b,c,d,radunit=False):
    h1,f1=dotpreptoline(a,b,c)
    h2,f2=dotpreptoline(d,b,c)
    cosf=numpy.dot(a-f1,d-f2)/(h1*h2)
    m=numpy.cross(a-f1,d-f2)
    sind=numpy.linalg.norm(m)/(h1*h2)
    abcd=(cosf/abs(cosf))*numpy.arcsin(sind)+(sind/abs(sind))*(math.pi/2)*(1-cosf/abs(cosf))
    if radunit==False:
        abcd=180.0*abcd/math.pi
    return abcd

def cart2inte(cart,sub=[[],[],[]],radunit=True,mat=False):
    global bohr
    xyzs,Bs,As,Ds=[],[],[],[]
    redint=[]
    N=int(len(cart)/3)
    redmat=numpy.zeros((len(sub[0])+len(sub[1])+len(sub[2]),3*N))
    if mat:
        import tensorflow as tf
        icm=numpy.zeros((3*N-6,3*N))
    for i in range(N):
        xyzs.append(numpy.array([cart[0+3*i],cart[1+3*i],cart[2+3*i]]))
 #   for i in range(3*N):
 #       cart[i]=float(cart[i])/bohr
    for i in range(N-1):
        coords=[]
        # Bonds.
        for j in range(6):
            if mat:
                coords.append(tf.Variable(initial_value=cart[j+3*i]))
            else:
                coords.append(cart[j+3*i])
        a=numpy.array([coords[0],coords[1],coords[2]])
        b=numpy.array([coords[3],coords[4],coords[5]])
        if mat:
            with tf.GradientTape(persistent=True) as calcB:
                B12=bond(a,b)
            for j in range(6):
                icm[i,j+3*i]=calcB.gradient(B12,coords[j])
        Bs.append(bond(a,b))
        # Bonds added by user.
        if len(sub[0])!=0:
            add=[]
            for j in range(len(sub[0])):
                for k in range(2):
                    for h in range(3):
                        if mat:
                            add.append(tf.Variable(initial_value=cart[3*sub[0][j][k]-3+h]))
                        else:
                            add.append(cart[3*sub[0][j][k]-3+h])
                a=numpy.array([add[0],add[1],add[2]])
                b=numpy.array([add[3],add[4],add[5]])
                if mat:
                    with tf.GradientTape(persistent=True) as calcaddB:
                        B12=bond(a,b)
                    for k in range(2):
                        for h in range(3):
                            icm[N-1-len(sub[0])+j,3*sub[0][j][k]-3+h]=calcaddB.gradient(B12,add[3*k+h])
                Bs[N-1-len(sub[0])+j]=bond(a,b)
        # Angles.
        if i<N-2:
            for j in range(3):
                if mat:
                    coords.append(tf.Variable(initial_value=cart[j+6+3*i]))
                else:
                    coords.append(cart[j+6+3*i])
            c=numpy.array([coords[6],coords[7],coords[8]])
            if mat:
                with tf.GradientTape(persistent=True) as calcA:
                    A123=angle(a,b,c,radunit)
                for j in range(9):
                    icm[N-1+i,j+3*i]=calcA.gradient(A123,coords[j])
            As.append(angle(a,b,c,radunit))
            # Angles added by user.
            if len(sub[1])!=0:
                add=[]
                for j in range(len(sub[1])):
                    for k in range(3):
                        for h in range(3):
                            if mat:
                                add.append(tf.Variable(initial_value=cart[3*sub[1][j][k]-3+h]))
                            else:
                                add.append(cart[3*sub[1][j][k]-3+h])
                    a=numpy.array([add[0],add[1],add[2]])
                    b=numpy.array([add[3],add[4],add[5]])
                    c=numpy.array([add[6],add[7],add[8]])
                    if mat:
                        with tf.GradientTape(persistent=True) as calcaddA:
                            A123=angle(a,b,c,radunit)
                        for k in range(3):
                            for h in range(3):
                                icm[2*N-3-len(sub[1])+j,3*sub[1][j][k]-3+h]=calcaddA.gradient(A123,add[3*k+h])
                    As[2*N-3-len(sub[1])+j]=angle(a,b,c,radunit)
        # Dihedrals.
        if i<N-3:
            for j in range(3):
                if mat:
                    coords.append(tf.Variable(initial_value=cart[j+9+3*i]))
                else:
                    coords.append(cart[j+9+3*i])
            d=numpy.array([coords[9],coords[10],coords[11]])
            if mat:
                with tf.GradientTape(persistent=True) as calcD:
                    D1234=dihedral(a,b,c,d,radunit)
                for j in range(12):
                    icm[2*N-3+i,j+3*i]=calcD.gradient(D1234,coords[j])
            Ds.append(dihedral(a,b,c,d,radunit))
            # Dihedrals added by user.
            if len(sub[2])!=0:
                add=[]
                for j in range(len(sub[2])):
                    for k in range(4):
                        for h in range(3):
                            if mat:
                                add.append(tf.Variable(initial_value=cart[3*sub[2][j][k]-3+h]))
                            else:
                                add.append(cart[3*sub[2][j][k]-3+h])
                    a=numpy.array([add[0],add[1],add[2]])
                    b=numpy.array([add[3],add[4],add[5]])
                    c=numpy.array([add[6],add[7],add[8]])
                    d=numpy.array([add[9],add[10],add[11]])
                    if mat:
                        with tf.GradientTape(persistent=True) as calcaddD:
                            D1234=dihedral(a,b,c,d,radunit)
                        for k in range(4):
                            for h in range(3):
                                icm[3*N-6-len(sub[2])+j,3*sub[2][j][k]-3+h]=calcaddD.gradient(D1234,add[3*k+h])
                    Ds[3*N-6-len(sub[2])+j]=dihedral(a,b,c,radunit)
    inte=Bs+As+Ds
    if mat:
 #       icm=numpy.delete(icm,[0,1,2,3,4,6],1)
        return inte,icm
    else:
        return inte

def abc2d(Ps,Bs,As,D1):
    # Solve P4 with (P1,P2,P3), (B1,B2,B3), (A1,A2) and D1 given.
    B1,B2,B3=Bs
    A1,A2=As
    # Set P2(0.0,0.0,0.0) for clarity.
    fP1=numpy.array([0.0,B1*numpy.sin(A1),B1*numpy.cos(A1)])
    fP2=numpy.array([0.0,0.0,0.0])
    fP3=numpy.array([0.0,0.0,B2])
    # Add one to avoid singular matrix problem.
    fPs=numpy.mat([fP1,fP2,fP3])+numpy.mat(numpy.ones((3,3)))
    # fP4 is:
    m=numpy.array([0.0,B3*numpy.sin(A2),-B3*numpy.cos(A2)])
    cD1,sD1=numpy.cos(D1),numpy.sin(D1)
    RzD=numpy.mat([[ cD1,sD1,0.0],
                   [-sD1,cD1,0.0],
                   [ 0.0,0.0,1.0]])
    fP4=fP3+m*RzD
    fP4+=1.0
    # Solve coordinate system rotation matrix.
    iPs=numpy.mat([Ps]-Ps[1])+numpy.mat(numpy.ones((3,3)))
    T=numpy.linalg.solve(iPs,fPs)
    Tinv=numpy.linalg.inv(T)
    # Then iP4 is:
    iP4=numpy.dot(fP4,Tinv).A[0]-1.0+Ps[1]
    return iP4

def inte2cart(inte,sub=[[],[],[]]):
    cart=[]
    N=int((len(inte)+6)/3)
    Bs,As,Ds=[],[],[]
    for i in range(len(inte)):
        if len(Bs)<N-1:
            Bs.append(inte[i])
        elif len(As)<N-2:
            As.append(inte[i])
        elif len(Ds)<N-3:
            Ds.append(inte[i])
    zp2=Bs[0]
    yp3=Bs[1]*numpy.sin(As[0])
    zp3=Bs[0]-Bs[1]*numpy.cos(As[0])
    p1=[0.0,0.0,0.0]
    p2=[0.0,0.0,zp2]
    p3=[0.0,yp3,zp3]
    cart.extend(p1)
    cart.extend(p2)
    cart.extend(p3)
    xyzs=[numpy.array(p1),
          numpy.array(p2),
          numpy.array(p3)]
    for i in range(N-3):
        a=xyzs[-3]
        b=xyzs[-2]
        c=xyzs[-1]
        d=abc2d([a,b,c],
                [Bs[i],Bs[i+1],Bs[i+2]],
                [As[i],As[i+1]],
                 Ds[i])
        xyzs.append(d)
        cart.extend(d.tolist())
    return cart

def geniX(xyzs,n=1):
    Xm=numpy.zeros((numpy.shape([xyzs[0]])))[0]
    for i in range(len(xyzs)):
        Xm+=numpy.array(xyzs[i])
    Xm=Xm/len(xyzs)
    iX=[Xm.tolist()]
    for i in range(len(xyzs)):
        Xi=numpy.array(xyzs[i])
        # iX.append(Xi.tolist()) # Whether to add starting points.
        dX=Xm-Xi
        if n>1:
            for j in range(n):
                Xj=Xi+j*dX/(n+1)
                iX.append(Xj.tolist())
    return iX

def linear_insert(xyzfile,n=1,internal=False):
    atm,xyzs=readxyzs(xyzfile)
    if len(xyzs)<2:
        print('Error! At least two structures should be provided for linear insert.')
        print('Exit program.')
        exit()
    if internal:
        from parameters import add_ints
        for i in range(len(xyzs)):
            xyzs[i]=cart2inte(xyzs[i],add_ints)
    iX=geniX(xyzs,n)
    print(f'Inserted structures: {n}')
    if internal:
        from parameters import add_ints
        for i in range(len(iX)):
            iX[i]=inte2cart(iX[i],add_ints)
    monitor=linear_search(iX)
    return iX[monitor]

def linear_search(iX):
    from parameters import layer,qm_software
    if layer=='qm':
        if qm_software=='g16':
            from rung16 import mecpstep
        elif qm_software=='orca':
            from runorca import mecpstep
        elif qm_software=='qchem':
            from runqchem import mecpstep
        elif qm_software=='pyscf':
            from runpyscf import mecpstep
    elif layer=='qmmm':
        if qm_software=='orca':
            from runchemsh_orca import mecpstep
        if qm_software=='pyscf':
            from runchemsh_pyscf import mecpstep
    ed=0.0
    monitor=0
    for i in range(len(iX)):
        e1,g1,e2,g2=mecpstep(iX(i))
        if ed<abs(e1-e2):
            ed=abs(e1-e2)
            monitor=i
    return monitor

