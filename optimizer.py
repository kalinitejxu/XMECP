import os,math,numpy
from parameters import *

if 'mecp' in algorithm:
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
elif algorithm=='opt':
    if layer=='qm':
        if qm_software=='pyscf':
            from runpyscf import optstep
    if layer=='qmmm':
        if qm_software=='orca':
            from runchemsh_orca import optstep
        if qm_software=='pyscf':
            from runchemsh_pyscf import optstep

### >>> Calculate force for MECP search >>> ###
def fgforce(e1,g1,e2,g2):
    ed=e2-e1
    gd=numpy.array(g2)-numpy.array(g1)
    norm_gd=gd/numpy.linalg.norm(gd)
    fvec=ed*norm_gd
    gvec=(g2-numpy.dot(norm_gd,g2)*norm_gd)

    mecpfrc=fvec+gvec
    return mecpfrc
    # Reference: Theor. Chem. Acc., 99, 95 (1998); Chem. Phys. Lett. 223, 269 (1994)

def bpupd(xs,ys,g1,g2):
    gd=g2-g1
    gm=(g1+g2)/2

    global gdiffs,gmeans
    gdiffs.append(gd)
    gmeans.append(gm)

    xk=gd/numpy.linalg.norm(gd)
    xs.append(xk)

    if ys==[]:
        yk=gm-(numpy.dot(xk,gm)/numpy.linalg.norm(xk)**2)*xk
        yk=yk/numpy.linalg.norm(yk)
        ys.append(yk)
    else:
        yk=(numpy.dot(ys[-1],xk)*xs[-2]-numpy.dot(xs[-2],xk)*ys[-1])/\
            math.sqrt(numpy.dot(ys[-1],xk)**2+numpy.dot(xs[-2],xk)**2)
        ys.append(yk)

    with open('bpupdvec.txt','w') as gen:
        dim=int(len(xk)/3)
        gen.write(f'Branching plane and gradient vectors. Dim={dim}.\n')
        gen.write('  x_x             x_y             x_z\n')
        for i in range(dim):
            for j in range(3):
                gen.write(f'{xk[j+3*i]:.8f}'.rjust(16))
            gen.write('\n')
        gen.write('  y_x             y_y             z_z\n')
        for i in range(dim):
            for j in range(3):
                gen.write(f'{yk[j+3*i]:.8f}'.rjust(16))
            gen.write('\n')
        gen.write(' g1_x            g1_y            g1_z\n')
        for i in range(dim):
            for j in range(3):
                gen.write(f'{g1[j+3*i]:.8f}'.rjust(16))
            gen.write('\n')
        gen.write(' g2_x            g2_y            g2_z\n')
        for i in range(dim):
            for j in range(3):
                gen.write(f'{g2[j+3*i]:.8f}'.rjust(16))
            gen.write('\n')

    return xs,ys
    # Reference: J. Comput. Theory Chem., 6, 1538 (2010)

def bpupdforce(e1,g1,e2,g2):
    g1,g2=numpy.array(g1),numpy.array(g2)
    global xs,ys

    xs,ys=bpupd(xs,ys,g1,g2)
    x,y=xs[-1],ys[-1]

    a,b=1.0,0.0
    if len(xs)>1:
        xp,yp=xs[-2],ys[-2]
        a=numpy.dot(y,xp)
        b=numpy.dot(y,yp)

    E=numpy.identity(numpy.shape(x)[0])
    xm=numpy.mat(x)
    ym=numpy.mat(y)
    P=E-xm.T*xm-ym.T*ym

    gd=numpy.mat(g2-g1).T
    gm=numpy.mat((g1+g2)/2).T

    # gm=gm/numpy.linalg.norm(gm)
    vec1=2*(e2-e1)*xm.T
    vec2=numpy.dot(P,gm)
    # vec2=gm-(float(numpy.dot(gm.T,xm.T))/numpy.linalg.norm(xm.T)**2)*xm.T-(float(numpy.dot(gm.T,ym.T))/numpy.linalg.norm(ym.T)**2)*ym.T

    bpupdfrc=vec1+vec2
    # print(f'debug {numpy.linalg.norm(vec1):5f} {numpy.linalg.norm(vec2):5f}')

    # Reference: J. Comput. Theory Chem., 6, 1538 (2010)
    return numpy.array(bpupdfrc.T)[0]

def pfforce(e1,g1,e2,g2):
    g1=numpy.array(g1)
    g2=numpy.array(g2)

    ed=e2-e1
    em=(e1+e2)/2
    gd=g2-g1
    gm=(g1+g2)/2

    global a,s
    a=5.0*627.511525
    s=5.0/627.511525
    W=1.0+ed**2/a**2
    G=em+s*a**2*math.log(W)

    gscale=2*s*ed/W
    grad=(0.5+gscale)*g2+(0.5-gscale)*g1

    global pfF
    pfF.append(G)

    return grad

# def pfforce(e1,g1,e2,g2):
#     g1=numpy.array(g1)
#     g2=numpy.array(g2)
# 
#     ed=e2-e1
#     em=(e1+e2)/2
#     gd=g2-g1
#     gm=(g1+g2)/2
# 
#     global a,s
#     vec1=gm
#     vec2=s*gd*(ed**2+2*a*ed)/(ed+a)**2
#     grad=vec1+vec2
# 
#     global pfF
#     pfF.append(em+s*ed**2/(ed+a))
#     return grad

def mecpforce(e1,g1,e2,g2):
    global e1s,g1s,e2s,g2s
    e1s.append(e1)
    g1s.append(g1)
    e2s.append(e2)
    g2s.append(g2)

    global algorithm
    if algorithm=='mecp_bpupd':
        frc0=bpupdforce(e1,g1,e2,g2)
    elif algorithm=='mecp_pf':
        frc0=pfforce(e1,g1,e2,g2)
    else:
        frc0=fgforce(e1,g1,e2,g2)
    return numpy.array(frc0)
### <<< Calculate force for MECP search <<< ###

def stepsize(X0,X1,factor=1):
    global maxstep
    dX=X1-X0
    dX*=factor
    stepsize_norm=numpy.linalg.norm(dX)
    if stepsize_norm>maxstep:
        print(f'Current stepsize: {stepsize_norm} is reduced to max stepsize.\n')
        dX=dX*maxstep/numpy.linalg.norm(dX)
    return X0+dX

def hessian_upd(B0,dF,dX,method='psb'):
    global algorithm
    if algorithm=='mecp_bpupd1':
        global gmeans
        dX=dX.T
        dh=numpy.mat(gmeans[-1]-gmeans[-2]).T
        B1=B0+(dh*dh.T)/(dh.T*dX)-(B0*dX)*(dX.T*Bk)/(dX.T*B0*dX)
        global xs,ys
        x,y=xs[-1],ys[-1]
        E=numpy.identity(numpy.shape(x)[0])
        xm=numpy.mat(x)
        ym=numpy.mat(y)
        P=E-xm.T*xm-ym.T*ym
        B1=P*B1*P
    else:
#       if algorithm=='mecp_pf':
#           method='bofill'
        if method=='psb':
            B1=B0+((dF.T-B0*dX.T)*dX+dX.T*(dF.T-B0*dX.T).T)/(dX*dX.T)-\
                  (float(dX*(dF.T-B0*dX.T))*dX.T*dX)/float(dX*dX.T)**2
        elif method=='bfgs':
            B1=B0+float(1.0/(dX*dF.T))*\
                 (float(1.0+(dF*B0*dF.T)/(dX*dF.T))*(dX.T*dX)-\
                                     (dX.T*dF*B0+B0*dF.T*dX))
        elif method=='bofill':
            B1=B0+((dF.T-B0*dX.T)*dX+dX.T*(dF.T-B0*dX.T).T)/(dX*dX.T)-\
                  (float(dX*(dF.T-B0*dX.T))*dX.T*dX)/float(dX*dX.T)**2
            B2=B0+(dF.T-B0*dX.T)*(dF.T-B0*dX.T).T/float((dF.T-B0*dX.T).T*dX.T)
            phi=float((dF.T-B0*dX.T).T*dX.T)**2/\
                float((dF.T-B0*dX.T).T*(dF.T-B0*dX.T))/float(dX*dX.T)
            B1=phi*B2+(1-phi)*B1
    return B1

def prop_bfgs(X0,F0,B0):
    scaledX=15.0
    global algorithm
    if 'mecp' in algorithm:
        scaledX=15.0
    elif algorithm=='opt':
        scaledX=0.01
    dX=-numpy.linalg.solve(B0,F0)
    if numpy.linalg.norm(dX)>0.1:
        dX=dX*0.1/numpy.linalg.norm(dX)
    X1=X0+scaledX*dX
    return X1
def bfgs(X0,F0,B0,step):
    print(f'Entering BFGS step {step}\n')
    X1=prop_bfgs(X0,F0,B0)
    X1=stepsize(X0,X1)
    dX=X1-X0
    if 'mecp' in algorithm:
        E1,G1,E2,G2=mecpstep(X1,step)
        F1=mecpforce(E1,G1,E2,G2)
        dF=numpy.mat(F1-F0)
        B1=hessian_upd(B0,dF,dX)
        return X1,F1,B1,E1,E2
    elif algorithm=='opt':
        E1,G1=optstep(X1,step)
        F1=numpy.array(G1)
        dF=numpy.mat(F1-F0)
        B1=hessian_upd(B0,dF,dX)
        return X1,F1,B1,E1,0.0

def prop_gdiis(Xs,Fs,Hs):
    # Produce a new geometry based on the GDIIS algorithm, see https://manual.q-chem.com/5.3/A1.S7.html
    nX=len(Xs)
    nA=numpy.shape(Xs[-1])[1]
    if len(Hs)!=nX:
        raise Exception('Runtime exception: H and X dimensions are different.')
    EMat=numpy.mat(numpy.zeros(shape=(nX,nA)))
    Hm=sum(Hs)/nX
    for i in range(nX):
        EMat[i]=(Hm.I*numpy.mat(Fs[i]).T).flatten()
    L1=numpy.mat(numpy.ones(nX))
    BMat=numpy.block([[EMat*EMat.T,L1.T],
                      [L1,numpy.mat([0])]])
    y=numpy.append(numpy.zeros(nX),1)
    c=numpy.linalg.solve(BMat,y)
    c=numpy.delete(c,-1)
    X2=numpy.mat(numpy.zeros(nA))
    F2=numpy.mat(numpy.zeros(nA))
    H2=numpy.mat(numpy.zeros(shape=(nA,nA)))
    for i in range(nX):
        X2+=Xs[i]*c[i]
        F2+=Fs[i]*c[i]
        H2+=Hs[i]*c[i]
    X1=X2-(Hm.I*numpy.mat(F2).T).flatten()
    return X1
def prop_gediis(Xs,Fs,Hs,Es):
    # Produce a new geometry based on the GEDIIS algorithm (J. Chem. Theory Comput. 2006, 2, 835-839)
    # Note that the energy to be minimized should not be E1 or E2, but produced from the pritimive function of G (not implemented yet)
    nX=len(Xs)
    nA=numpy.shape(Xs[-1])[1]
    if len(Hs)!=nX:
        raise Exception('Runtime exception: H and X numbers are different.')
    EMat=numpy.mat(numpy.zeros(shape=(nX,nX)))
    for i in range(nX):
        EMat[i,i]=0
        for j in range(i+1,nX):
            EMat[i,j]=-1*numpy.dot(Fs[i]-Fs[j],
                                   numpy.array((Xs[i]-Xs[j])).flatten())  # Gs is force rather than gradient
            EMat[j,i]=EMat[i,j]
    L1=numpy.mat(numpy.ones(nX))
    BMat=numpy.block([[EMat,L1.T],\
                      [L1,numpy.mat([0])]])
    y=numpy.append(-1*numpy.array(Es),1)
    c=numpy.linalg.solve(BMat,y)
    c=numpy.delete(c,-1)  # delete the last element, lambda, in c vector
    X2=numpy.mat(numpy.zeros(nA))
    F2=numpy.mat(numpy.zeros(nA))
    for i in range(nX):
        X2+=Xs[i]*c[i]
        F2+=Fs[i]*c[i]
    X1=X2+F1
    return X1
def gdiis(Xs,Fs,Bs,Es,step,gediis=False):
    global conver
    thresh_rms_g=conver[4]
    reduced_factor=0.5
    
    X1=prop_gdiis(Xs,Fs,Bs)
    if gediis:
        X2=prop_gediis(Xs,Fs,Bs,Es)
        X1=X1*0.5+X2*0.5
        print(f'\nEntering GEDIIS Step {step}\n')
    else:
        print(f'\nEntering GDIIS Step {step}\n')

    factor=1.0
    if numpy.linalg.norm(Fs)<thresh_rms_g*10:
        factor=reduced_factor
    X1=stepsize(Xs[-1],X1,factor)

    dX=X1-Xs[-1]
    if 'mecp' in algorithm:
        E1,G1,E2,G2=mecpstep(X1,step)
        F1=mecpforce(E1,G1,E2,G2)
        dF=numpy.mat(F1-Fs[-1])
        B1=hessian_upd(Bs[-1],dF,dX)
        return X1,F1,B1,E1,E2
    elif algorithm=='opt':
        E1,G1=optstep(X1,step)
        F1=numpy.array(G1)
        dF=numpy.mat(F1-Fs[-1])
        B1=hessian_upd(Bs[-1],dF,dX)
        return X1,F1,B1,E1,0.0

def geom_upd(X1,cfile,chemsh_template):
    numatm,atm,xyz=readc(cfile)
    act_region=active_region(chemsh_template)
    for i in range(len(act_region)):
        for j in range(3):
            xyz[3*act_region[i]-3+j]=X1[0,3*i+j]
    return atm,xyz

def opt(X0):
    step=1
    global opt_method,maxcycle,restart,restart_from
    global algorithm,layer

    # Cut out active region for optimization.
    if layer=='qmmm':
        X0=geom.active_xyz(X0)

    # Run optimization.
    X0=numpy.mat(X0)
    B0=numpy.eye(numpy.shape(X0)[1])
    X1,F1,B1,E1,E2=[None,None,None,0,0]
    Xs,Bs,Fs,Es=[[],[],[],[]]
    
    F0=None
    # Run first step to get F0.
    if 'mecp' in algorithm:
        E1,G1,E2,G2=mecpstep(X0,step)
        print('Initial structure')
        print(f'Energy of state 1:{f"{E1:.8f}".rjust(16)} a.u.')
        print(f'Energy of state 2:{f"{E2:.8f}".rjust(16)} a.u.')
        print(f'Energy difference:{f"{(E2-E1):.8f}".rjust(16)} a.u.')
        F0=mecpforce(E1,G1,E2,G2)
    elif algorithm=='opt':
        E1,G1=optstep(X0,step)
        print('Initial structure')
        print(f'Energy:{f"{E1:.8f}".rjust(16)}')
        F0=numpy.array(G1)

    # Entering iteration.
    while True:
        step+=1
        if opt_method=='bfgs' or step<3:
            X1,F1,B1,E1,E2=bfgs(X0,F0,B0,step)
        elif opt_method=='gdiis':
            X1,F1,B1,E1,E2=gdiis(Xs,Fs,Bs,Es,step)
        elif opt_method=='gediis':
            X1,F1,B1,E1,E2=gdiis(Xs,Fs,Bs,Es,step,gediis=True)
        if len(Xs)>3:
            Xs.pop(0)
            Bs.pop(0)
            Fs.pop(0)
            Es.pop(0)
        Xs.append(X1)
        Bs.append(B1)
        Fs.append(F1)
        Es.append(E1)

        if isconver(E1,E2,X0,X1,F1):
            break

        if step==maxcycle:
            print('Maximium number of steps exceeded.')
            break
        X0=X1
        F0=F1
        B0=B1

def isconver(E1,E2,X0,X1,G1):
    global algorithm,layer
    if layer=='qm':
        if algorithm=='mecp_pf':
            return isconver_pf(E1,E2,X0,X1,G1)
        else:
            return isconver_standard(E1,E2,X0,X1,G1)
    if layer=='qmmm':
        import geom
        qmG1=geom.qm_from_active(G1)
        qmX0=geom.qm_from_active(numpy.matrix.tolist(X0)[0])
        qmX1=geom.qm_from_active(numpy.matrix.tolist(X1)[0])
        qmX0=numpy.mat(qmX0)
        qmX1=numpy.mat(qmX1)
        # Check QM region first.
        qmConverged=False
        save_stdout=sys.stdout
        with open('trash.txt','w') as sys.stdout:
            if algorithm=='mecp_pf':
                qmConverged=isconver_pf(E1,E2,qmX0,qmX1,qmG1)
            else:
                qmConverged=isconver_standard(E1,E2,qmX0,qmX1,qmG1)
        sys.stdout=save_stdout
        if qmConverged:
            print('QM region has already converged.')
        else:
            print('QM region has not converged yet.')
        # Check the whole system and print convergence infomation.
        if algorithm=='mecp_pf':
            return isconver_pf(E1,E2,X0,X1,G1)
        else:
            return isconver_standard(E1,E2,X0,X1,G1)

def isconver_standard(E1,E2,X0,X1,G1):
    dim=len(G1)
    global thresh_de,thresh_rms_g,thresh_max_g,thresh_rms,thresh_max_dis
    # Calculate convergence criterions.
    de=E2-E1
    rms_g=numpy.linalg.norm(G1)/math.sqrt(dim)
    max_g=0.0

    for i in range(0,dim,3):
        g=abs(G1[i])
        if g>max_g:
            max_g=g
    rms=numpy.linalg.norm(X0-X1)/math.sqrt(dim)
    max_dis=0.0
    for i in range(int(dim/3)):
        dis=0.0
        for j in range(3):
            dis+=(X0[0,3*i+j]-X1[0,3*i+j])**2
        dis=math.sqrt(dis)
        if dis>max_dis:
            max_dis=dis

    # Compare with convergence criterions.
    de_conv,rms_conv,max_dis_conv,max_g_conv,rms_g_conv=False,False,False,False,False
    yes1,yes2,yes3,yes4,yes5='NO','NO','NO','NO','NO'
    if abs(de)<thresh_de:
        de_conv=True
        yes1='YES'
    if abs(rms_g)<thresh_rms_g:
        rms_g_conv=True
        yes2='YES'
    if abs(max_g)<thresh_max_g:
        max_g_conv=True
        yes3='YES'
    if abs(rms)<thresh_rms:
        rms_conv=True
        yes4='YES'
    if abs(max_dis)<thresh_max_dis:
        max_dis_conv=True
        yes5='YES'

    # Print convergence infomation of current cycle.
    print('------------------------------------------------------')
    print('           >>> Summary for Current Step <<<           ')
    print('------------------------------------------------------')
    if 'mecp' in algorithm:
        print(f'E1 = {E1:.8f}\nE2 = {E2:.8f}')
    else:
        print(f'Energy = {E1:.8f}')
    print('------------------------------------------------------')
    if 'mecp' in algorithm:
        print(f'           Energy Gap     {abs(de):.6f}     {thresh_de:.6f}     {yes1}')
    print(f'         RMS Gradient     {rms_g:.6f}     {thresh_rms_g:.6f}     {yes2}')
    print(f'    Maximium Gradient     {max_g:.6f}     {thresh_max_g:.6f}     {yes3}')
    print(f'     RMS Displacement     {rms:.6f}     {thresh_rms:.6f}     {yes4}')
    print(f'Maximium Displacement     {max_dis:.6f}     {thresh_max_dis:.6f}     {yes5}')
    print('------------------------------------------------------')
    if rms_conv and max_dis_conv and max_g_conv and rms_g_conv:
        if ('mecp' in algorithm and de_conv) or algorithm=='opt':
            print('Optimization Completed!')
            return True
    else:
        return False

def isconver_pf(E1,E2,X0,X1,G1):
    global natms
    global thresh_de,pf_thresh_step,pf_thresh_grad

    # Calculate convergence criterions.
    de=E2-E1
    global a,s
    global g1s,g2s
    g1=numpy.array(g1s[-1])
    g2=numpy.array(g2s[-1])
    gd=g2-g1
    gm=(g1+g2)/2
    Z=numpy.dot(gd,gm)/numpy.dot(gm,gm)
    pfE=None
    try:
        pfE=a*(math.sqrt((s*Z)/(1+s*Z))-1)
    except:
        pfE='******'
    global pfF
    func1=pfF[-2]-pfF[-1]
    
 #    dG=gd*((de**2+2*a*de)/(de+a)**2)
 #    u=dG/numpy.linalg.norm(dG)
 #    func2=numpy.dot(-G1,u)/s
 #    func3=numpy.linalg.norm(-G1-numpy.dot(-G1,u)*u)

    a=5.0*627.511525
    s=5.0/627.511525
    W=1.0+de**2/a**2
    G=(E1+E2)/2+s*a**2*math.log(W)

    gscale=2*s*de/W
    dG=(0.5+gscale)*g2+(0.5-gscale)*g1
    u=dG/numpy.linalg.norm(dG)

    func2=numpy.dot(-G1,u)
    func3=numpy.linalg.norm(-G1-numpy.dot(-G1,u)*u)

    rms=numpy.linalg.norm(X0-X1)/math.sqrt(natms*3)
    max_dis=0.0
    for i in range(natms):
        dis=0.0
        for j in range(3):
            dis+=(X0[0,3*i+j]-X1[0,3*i+j])**2
        dis=math.sqrt(dis)
        if dis>max_dis:
            max_dis=dis

    # Compare with convergence criterions.
    func1_conv,func2_conv,func3_conv,rms_conv,max_dis_conv=False,False,False,False,False
    yes1,yes2,yes3,yes4,yes5='NO','NO','NO','NO','NO'
    if abs(func1)<pf_thresh_step:
        func1_conv=True
        yes1='YES'
    if abs(func2)<pf_thresh_grad:
        func2_conv=True
        yes2='YES'
    if abs(func3)<pf_thresh_grad:
        func3_conv=True
        yes3='YES'
    if abs(rms)<thresh_rms:
        rms_conv=True
        yes4='YES'
    if abs(max_dis)<thresh_max_dis:
        max_dis_conv=True
        yes5='YES'

    # Print convergence infomation of current cycle.
    print('------------------------------------------------------')
    print('          Penalty Function MECI Optimization          ')
    print('           >>> Summary for Current Step <<<           ')
    print('------------------------------------------------------')
    print(f'E1 = {E1}\nE2 = {E2}')
    print(f'Energy Gap: {de:5f}')
    try:
        print(f'pfE: {pfE:5f}')
    except:
        print(f'pfE: {pfE}')
    print('------------------------------------------------------')
    print(f'         PF Function1     {func1:5f}     {pf_thresh_step:5f}     {yes1}')
    print(f'         PF Function2     {func2:5f}     {pf_thresh_grad:5f}     {yes2}')
    print(f'         PF Function3     {func3:5f}     {pf_thresh_grad:5f}     {yes3}')
 #    print(f'     RMS Displacement     {rms:5f}     {thresh_rms:5f}     {yes4}')
 #    print(f'Maximium Displacement     {max_dis:5f}     {thresh_max_dis:5f}     {yes5}')
    print('------------------------------------------------------')
    if func1_conv and func2_conv and func3_conv:
        if 1==1:
 #        if abs(de)<thresh_de:
            print('Optimization Completed!')
            return True
        else:
            s*=2.0
            print(f'Increse sigma to {s:2f} for following steps.')
            return False
    else:
        return False

