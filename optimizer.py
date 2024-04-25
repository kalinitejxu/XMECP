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

def stepsize(oldX,newX,factor=1):
    global maxstep
    dX=newX-oldX
    dX*=factor
    stepsize_norm=numpy.linalg.norm(dX)
    if stepsize_norm>maxstep:
        print(f'Reduce current stepsize {stepsize_norm:.4f} to max stepsize {maxstep:.4f}.\n')
        dX=dX*maxstep/numpy.linalg.norm(dX)
    return oldX+dX

def hessian_upd(old_hess,dF,dX,method='psb'):
    global algorithm
    if algorithm=='mecp_bpupd1':
        global gmeans
        dX=dX.T
        dh=numpy.mat(gmeans[-1]-gmeans[-2]).T
        new_hess=old_hess+(dh*dh.T)/(dh.T*dX)-(old_hess*dX)*(dX.T*Bk)/(dX.T*old_hess*dX)
        global xs,ys
        x,y=xs[-1],ys[-1]
        E=numpy.identity(numpy.shape(x)[0])
        xm=numpy.mat(x)
        ym=numpy.mat(y)
        P=E-xm.T*xm-ym.T*ym
        new_hess=P*new_hess*P
    else:
#       if algorithm=='mecp_pf':
#           method='bofill'
        if method=='psb':
            new_hess=old_hess+((dF.T-old_hess*dX.T)*dX+dX.T*(dF.T-old_hess*dX.T).T)/(dX*dX.T)-\
                  (float(dX*(dF.T-old_hess*dX.T))*dX.T*dX)/float(dX*dX.T)**2
        elif method=='bfgs':
            new_hess=old_hess+float(1.0/(dX*dF.T))*\
                 (float(1.0+(dF*old_hess*dF.T)/(dX*dF.T))*(dX.T*dX)-\
                                     (dX.T*dF*old_hess+old_hess*dF.T*dX))
        elif method=='bofill':
            new_hess=old_hess+((dF.T-old_hess*dX.T)*dX+dX.T*(dF.T-old_hess*dX.T).T)/(dX*dX.T)-\
                  (float(dX*(dF.T-old_hess*dX.T))*dX.T*dX)/float(dX*dX.T)**2
            tmp_hess=old_hess+(dF.T-old_hess*dX.T)*(dF.T-old_hess*dX.T).T/float((dF.T-old_hess*dX.T).T*dX.T)
            phi=float((dF.T-old_hess*dX.T).T*dX.T)**2/\
                float((dF.T-old_hess*dX.T).T*(dF.T-old_hess*dX.T))/float(dX*dX.T)
            new_hess=phi*tmp_hess+(1-phi)*new_hess
    return new_hess

def bfgs(oldX,oldF,oldH,step):
    print(f'Entering BFGS step {step}\n')
    scaledX=15.0
    global algorithm
    if 'mecp' in algorithm:
        scaledX=15.0
    elif algorithm=='opt':
        scaledX=0.01
    dX=-numpy.linalg.solve(oldH,oldF)
    if numpy.linalg.norm(dX)>0.1:
        dX=dX*0.1/numpy.linalg.norm(dX)
    newX=stepsize(oldX,oldX+scaledX*dX)
    dX=newX-oldX

    if 'mecp' in algorithm:
        E1,G1,E2,G2=mecpstep(newX,step)
        newF=mecpforce(E1,G1,E2,G2)
        dF=numpy.mat(newF-oldF)
        newH=hessian_upd(oldH,dF,dX)
        return newX,newF,newH,E1,E2
    elif algorithm=='opt':
        E1,G1=optstep(newX,step)
        newF=numpy.array(G1)
        dF=numpy.mat(newF-oldF)
        newH=hessian_upd(oldH,dF,dX)
        return newX,newF,newH,E1,0.0

def gdiis(Xhist,Fhist,Hhist,Ehist,step,gediis=False):
    global conver
    thresh_rms_g=conver[4]

    # generate newX for gdiis
    nX=len(Xhist)
    nA=numpy.shape(Xhist[-1])[1]
    matrix_ene=numpy.mat(numpy.zeros(shape=(nX,nA)))
    Hm=sum(Hhist)/nX
    for i in range(nX):
        matrix_ene[i]=(Hm.I*numpy.mat(Fhist[i]).T).flatten()
    matrix_1=numpy.mat(numpy.ones(nX))
    matrix_hess=numpy.block([[matrix_ene*matrix_ene.T,matrix_1.T],
                      [matrix_1,numpy.mat([0])]])
    y=numpy.append(numpy.zeros(nX),1)
    c=numpy.linalg.solve(matrix_hess,y)
    c=numpy.delete(c,-1)
    tmpX=numpy.mat(numpy.zeros(nA))
    tmpF=numpy.mat(numpy.zeros(nA))
    tmpH=numpy.mat(numpy.zeros(shape=(nA,nA)))
    for i in range(nX):
        tmpX+=Xhist[i]*c[i]
        tmpF+=Fhist[i]*c[i]
        tmpH+=Hhist[i]*c[i]
    newX=tmpX-(Hm.I*numpy.mat(tmpF).T).flatten()
    # update newX for gediis
    if gediis:
        nX=len(Xhist)
        nA=numpy.shape(Xhist[-1])[1]
        matrix_ene=numpy.mat(numpy.zeros(shape=(nX,nX)))
        for i in range(nX):
            matrix_ene[i,i]=0
            for j in range(i+1,nX):
                matrix_ene[i,j]=-1*numpy.dot(Fhist[i]-Fhist[j],
                                       numpy.array((Xhist[i]-Xhist[j])).flatten())
                matrix_ene[j,i]=matrix_ene[i,j]
        matrix_1=numpy.mat(numpy.ones(nX))
        matrix_hess=numpy.block([[matrix_ene,matrix_1.T],\
                          [matrix_1,numpy.mat([0])]])
        y=numpy.append(-1*numpy.array(Ehist),1)
        c=numpy.linalg.solve(matrix_hess,y)
        c=numpy.delete(c,-1)
        tmpX=numpy.mat(numpy.zeros(nA))
        tmpF=numpy.mat(numpy.zeros(nA))
        for i in range(nX):
            tmpX+=Xhist[i]*c[i]
            tmpF+=Fhist[i]*c[i]
        newX=(newX+tmpX+newF)/2
        print(f'\nEntering GEDIIS Step {step}\n')
    else:
        print(f'\nEntering GDIIS Step {step}\n')

    factor=1.0
    if numpy.linalg.norm(Fhist)<thresh_rms_g*10:
        factor=0.5
    newX=stepsize(Xhist[-1],newX,factor)

    dX=newX-Xhist[-1]
    if 'mecp' in algorithm:
        E1,G1,E2,G2=mecpstep(newX,step)
        newF=mecpforce(E1,G1,E2,G2)
        dF=numpy.mat(newF-Fhist[-1])
        newH=hessian_upd(Hhist[-1],dF,dX)
        return newX,newF,newH,E1,E2
    elif algorithm=='opt':
        E1,G1=optstep(newX,step)
        newF=numpy.array(G1)
        dF=numpy.mat(newF-Fhist[-1])
        newH=hessian_upd(Hhist[-1],dF,dX)
        return newX,newF,newH,E1,0.0

def geom_upd(newX,cfile,chemsh_template):
    numatm,atm,xyz=readc(cfile)
    act_region=active_region(chemsh_template)
    for i in range(len(act_region)):
        for j in range(3):
            xyz[3*act_region[i]-3+j]=newX[0,3*i+j]
    return atm,xyz

def opt(oldX):
    step=1
    global opt_method,maxcycle,restart,restart_from
    global algorithm,layer

    # Cut out active region for optimization.
    if layer=='qmmm':
        oldX=geom.active_xyz(oldX)

    # Run optimization.
    oldX=numpy.mat(oldX)
    oldH=numpy.eye(numpy.shape(oldX)[1])
    newX,newF,newH,E1,E2=[None,None,None,0,0]
    Xhist,Bhist,Fhist,Ehist=[[],[],[],[]]
    
    oldF=None
    # Run first step to get oldF.
    if 'mecp' in algorithm:
        E1,G1,E2,G2=mecpstep(oldX,step)
        print('Initial structure')
        print(f'Energy of state 1:{f"{E1:.8f}".rjust(16)} a.u.')
        print(f'Energy of state 2:{f"{E2:.8f}".rjust(16)} a.u.')
        print(f'Energy difference:{f"{(E2-E1):.8f}".rjust(16)} a.u.')
        oldF=mecpforce(E1,G1,E2,G2)
    elif algorithm=='opt':
        E1,G1=optstep(oldX,step)
        print('Initial structure')
        print(f'Energy:{f"{E1:.8f}".rjust(16)}')
        oldF=numpy.array(G1)

    # Entering iteration.
    while True:
        step+=1
        if opt_method=='bfgs' or step<3:
            newX,newF,newH,E1,E2=bfgs(oldX,oldF,oldH,step)
        elif opt_method=='gdiis':
            newX,newF,newH,E1,E2=gdiis(Xhist,Fhist,Bhist,Ehist,step)
        elif opt_method=='gediis':
            newX,newF,newH,E1,E2=gdiis(Xhist,Fhist,Bhist,Ehist,step,gediis=True)
        if len(Xhist)>3:
            Xhist.pop(0)
            Bhist.pop(0)
            Fhist.pop(0)
            Ehist.pop(0)
        Xhist.append(newX)
        Bhist.append(newH)
        Fhist.append(newF)
        Ehist.append(E1)

        if isconver(E1,E2,oldX,newX,newF):
            break

        if step==maxcycle:
            print('Maximium number of steps exceeded.')
            break
        oldX=newX
        oldF=newF
        oldH=newH

def isconver(E1,E2,oldX,newX,G1):
    global algorithm,layer
    if layer=='qm':
        if algorithm=='mecp_pf':
            return isconver_pf(E1,E2,oldX,newX,G1)
        else:
            return isconver_standard(E1,E2,oldX,newX,G1)
    if layer=='qmmm':
        import geom
        qmG1=geom.qm_from_active(G1)
        qmoldX=geom.qm_from_active(numpy.matrix.tolist(oldX)[0])
        qmnewX=geom.qm_from_active(numpy.matrix.tolist(newX)[0])
        qmoldX=numpy.mat(qmoldX)
        qmnewX=numpy.mat(qmnewX)
        # Check QM region first.
        qmConverged=False
        save_stdout=sys.stdout
        with open('trash.txt','w') as sys.stdout:
            if algorithm=='mecp_pf':
                qmConverged=isconver_pf(E1,E2,qmoldX,qmnewX,qmG1)
            else:
                qmConverged=isconver_standard(E1,E2,qmoldX,qmnewX,qmG1)
        sys.stdout=save_stdout
        if qmConverged:
            print('QM region has already converged.')
        else:
            print('QM region has not converged yet.')
        # Check the whole system and print convergence infomation.
        if algorithm=='mecp_pf':
            return isconver_pf(E1,E2,oldX,newX,G1)
        else:
            return isconver_standard(E1,E2,oldX,newX,G1)

def isconver_standard(E1,E2,oldX,newX,G1):
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
    rms=numpy.linalg.norm(oldX-newX)/math.sqrt(dim)
    max_dis=0.0
    for i in range(int(dim/3)):
        dis=0.0
        for j in range(3):
            dis+=(oldX[0,3*i+j]-newX[0,3*i+j])**2
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

def isconver_pf(E1,E2,oldX,newX,G1):
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

    rms=numpy.linalg.norm(oldX-newX)/math.sqrt(natms*3)
    max_dis=0.0
    for i in range(natms):
        dis=0.0
        for j in range(3):
            dis+=(oldX[0,3*i+j]-newX[0,3*i+j])**2
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

