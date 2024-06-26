import os,numpy,geom
from parameters import *

### >>> Notice >>> ###

# Set current directory as QCSCRATCH into your job script, for example:
# QCSCRATCH=/home/jwxu/xmecp (where you run XMECP with Q-Chem)
# export QCSCRATCH=/home/jwxu/xmecp

### <<< Notice <<< ###

### >>> Special for Q-Chem spin-flip framework >>> ###
def spin_flip(xyzfile,step,incnroots=0):
    global job_name,system_name,nproc,qm_charge,qm_mult,qm_template
    global mecp_wordlist
    global sf_ref_mult,sf_thresh
    global atm

    print(f'Spin-flip (flip-down) is turned on. Multiplicity of reference state is {sf_ref_mult}.')
    # Generate a spin-flip TD-DFT single point task to decide iroot and jroot.
    inp=system_name+'_step'+str(step)+'_sf.inp'
    out=system_name+'_step'+str(step)+'_sf.out'

    qc_inp4mecp=open(inp,'w')
    qc_inp4mecp.write('# Q-Chem input file generated by XMECP\n')
    qc_inp4mecp.write(f'$molecule\n  read {xyzfile}\n$end\n')

    # $rem section
    qc_inp4mecp.write('\n$rem')
    qc_inp4mecp.write('\n  Spin_Flip True')
    qc_inp4mecp.write('\n  CIS_Max_Cycles 128')
    nroots=3+5*incnroots+max(state1_root,state2_root)
    print(f'Search for required states in the lowest {nroots} states.')
    qc_inp4mecp.write(f'\n  CIS_N_Roots {nroots}\n')

    flag=False
    extrakeys=open(qm_template,'r')
    for line in extrakeys.readlines():
        L=line.split()
        if len(L)!=0:
            if L[0]=='$rem':
                flag=True
                continue
            if flag and L[0]!='$end':
                qc_inp4mecp.write(line)
            if L[0]=='$end':
                flag=False
    extrakeys.close()
    qc_inp4mecp.write('$end\n')

    # Other sections required by user.
    flag=False
    extrakeys=open(qm_template,'r')
    for line in extrakeys.readlines():
        L=line.split()
        if len(L)!=0:
            probe=line[0]
            if probe[0]=='$' and L[0]!='$end' and L[0]!='$rem':
                flag=True
            if flag:
                qc_inp4mecp.write(line)
            if L[0]=='$end':
                flag=False
    extrakeys.close()
    qc_inp4mecp.close()

    # Now run spin-flip TD-DFT single point task.
    global nproc
    with os.popen(f'qchem -nt {nproc} {inp} {out} {scratch_path}/save_step{step}_sf') as run:
        null=run.read()

    # Search for required states.
    flag=False
    s2val,ijroots=[],[]
    sfsp=open(out,'r')
    for line in sfsp.readlines():
        L=line.split()
        if len(L)!=0:
            if L[0]=='SF-DFT' and L[1]=='Excitation' and L[2]=='Energies':
                flag=True
                continue
            if flag:
                if L[0]=='<S**2>':
                    s2val.append(float(L[-1]))
                    if len(s2val)==nroots:
                        break
                    continue
    flags=[False,False]
    for i in range(2):
        count=0
        s=(mecp_wordlist[0+i]-1)/2
        s2=s*(s+1)
        for j in range(len(s2val)):
            print (f'debug: {s2-(2-sf_thresh)} {s2val[j]} {s2+sf_thresh}')
            if s2-(2-sf_thresh) <= s2val[j] <= s2+sf_thresh:
                count+=1
                if count==mecp_wordlist[2+i]:
                    ijroots.append(1+j)
                    flags[i]==True
                    print(f'Found target state {1+i} as spin-flip excited state {1+j}.')
                    print(f'Expected value for spin-square operator is: {s2val[j]:3f}')
                    break
                elif j==len(s2val)-1:
                    print(f'Target state {1+i} was not found in current CIS space composed of {len(s2val)} states.')
    return ijroots
### <<< Special for ORCA spin-flip framework <<< ###

### >>> Generate ORCA input file for two-states MECP optimization >>> ###
def inp4mecp(X0,step=1):
    global job_name,system_name,nproc,qm_charge,qm_mult,uks,qm_template
    global atm
    global sf,sf_ref_mult

    xyz=numpy.matrix.tolist(X0)[0]
    xyzfile=system_name+'_step'+str(step)+'.geom'
    if sf:
        geom.qchemxyz(atm,xyz,qm_charge,sf_ref_mult,xyzfile)

    ijroots=[]
    incnroots=0
    if sf:
        while len(ijroots)!=2:
            ijroots=spin_flip(xyzfile,step,incnroots)
            if len(ijroots)!=2:
                incnroots+=1
                print('Increase CIS space for another 5 states.')
            if incnroots==5:
                print('Target state was not found with 20 more states expanding the CIS space.')
                print('Your system may not be suitable for spin-flip calculation. Exit program!')
                exit()

    global mecp_wordlist
    for i in range(2):
        # Q-Chem input files are generated in current folder and filenames are set to ${system_name}_stepN_state1.gjf and ${system_name}_stepN_state2.gjf.
        if sf==False:
            xyzfile=system_name+'_step'+str(step)+'_state'+str(i+1)+'.geom'
            geom.qchemxyz(atm,xyz,qm_charge,qm_mult,xyzfile)

        filename=system_name+'_step'+str(step)+'_state'+str(i+1)+'.inp'
        qc_inp4mecp=open(filename,'w')
        qc_inp4mecp.write('! Q-Chem input file generated by XMECP\n')
        qc_inp4mecp.write(f'$molecule\n  read {xyzfile}\n$end\n')

        # $rem section
        qc_inp4mecp.write('\n$rem')
        qc_inp4mecp.write('\n  JobType Force')

        tddft=False
        if sf or\
           (mecp_wordlist[0+i]==qm_mult and mecp_wordlist[2+i]!=0) or\
           (mecp_wordlist[0+i]!=qm_mult and mecp_wordlist[2+i]==1 and mecp_wordlist[4+i]==True) or\
           (mecp_wordlist[0+i]!=qm_mult and mecp_wordlist[2+i]>=2):
            tddft=True
            iroot=mecp_wordlist[2+i]
            qc_inp4mecp.write('\n  CIS_Max_Cycles 128')
            if sf:
                qc_inp4mecp.write('\n  Spin_Flip True')
                iroot=ijroots[i]
            else:
                if qm_mult==1 and mecp_wordlist[0+i]!=3:
                    qc_inp4mecp.write('\n  CIS_Triplets False')
            qc_inp4mecp.write(f'\n  CIS_State_Deriv {iroot}')
            qc_inp4mecp.write(f'\n  CIS_N_Roots {3+iroot}')

        if sf:
            qc_inp4mecp.write('\n  SCF_Guess Read')

        qc_inp4mecp.write('\n')
        flag=False
        extrakeys=open(qm_template,'r')
        for line in extrakeys.readlines():
            L=line.split()
            if len(L)!=0:
                if L[0]=='$rem':
                    flag=True
                    continue
                if flag and L[0]!='$end':
                    qc_inp4mecp.write(line)
                if L[0]=='$end':
                    flag=False
        extrakeys.close()
        qc_inp4mecp.write('$end\n')

        # Other sections required by user.
        flag=False
        extrakeys=open(qm_template,'r')
        for line in extrakeys.readlines():
            L=line.split()
            if len(L)!=0:
                probe=line[0]
                if probe[0]=='$' and L[0]!='$end' and L[0]!='$rem':
                    flag=True
                if flag:
                    qc_inp4mecp.write(line)
                if L[0]=='$end':
                    flag=False
        extrakeys.close()
    return
### <<< Generate Q-Chem input file for two-states MECP optimization <<< ###

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

def energy(energyfile,td=False,iroot=1):
    ene=0.0
    with open(energyfile,'r') as readene:
        for line in readene.readlines():
            L=line.split()
            if len(L)!=0:
                if td==False and 'Total energy in the final basis set' in line:
                    ene=float(L[-1])
                    break
                if td and 'Total energy for state' in line:
                    if L[4]==str(iroot)+':':
                        if sasf:
                            ene=float(L[-1])
                        else:
                            ene=float(L[-2])
                        break
    return ene

def gradient(gradfile,td=False):
    global sasf
    if sasf:
        return sasf_grad(gradfile)

    global natms

    flag=False
    grad,gradx,grady,gradz=[],[],[],[]

    counter=0
    dim=3*natms
    mark=int(natms/6)
    res=natms-mark*6

    with open(gradfile,'r') as readgrad:
        for line in readgrad.readlines():
            l=line.split()
            if len(l)!=0:
                if td and 'Gradient of the state energy' in line:
                    flag=True
                elif td==False and 'Gradient of SCF Energy' in line:
                    flag=True
                if flag:
                    if len(gradz)==natms:
                        break
                    if counter<mark:
                        if len(l)==7:
                            if l[0]=='1':
                                for i in range(6):
                                    gradx.append(float(l[1+i]))
                                continue
                            elif l[0]=='2':
                                for i in range(6):
                                    grady.append(float(l[1+i]))
                                continue
                            elif l[0]=='3':
                                for i in range(6):
                                    gradz.append(float(l[1+i]))
                                counter+=1
                                continue
                    elif counter==mark and res!=0:
                        if len(l)==res+1:
                            if l[0]=='1':
                                for i in range(res):
                                    gradx.append(float(l[1+i]))
                                continue
                            elif l[0]=='2':
                                for i in range(res):
                                    grady.append(float(l[1+i]))
                                continue
                            elif l[0]=='3':
                                for i in range(res):
                                    gradz.append(float(l[1+i]))
                                continue
    for i in range(natms):
        grad.extend([gradx[i],grady[i],gradz[i]])
    return grad

def sasf_grad(gradfile):
    grad=[]
    flag=False
    with open(gradfile,'r') as readgrad:
        for line in readgrad.readlines():
            L=line.split()
            if len(L)>=2:
                if 'SA-SF-DFT Gradient' in line:
                    flag=True
                    continue
                if flag:
                    if is_number(L[-3]) and is_number(L[-2]) and is_number(L[-1]):
                        for i in range(3):
                            grad.append(-float(L[-(3-i)]))
                    elif L[0]=='SA-SF-DFT':
                        break
    return grad

def mecpstep(XNew,step):
    print('Optimization steps: '+str(step))

    global scratch_path,xyz_backup,ene_backup
    if step==1:
        os.system('mkdir '+scratch_path)
        os.system('rm -f '+xyz_backup+' '+ene_backup)
        with open(ene_backup,'w') as gen:
            gen.write('MECP Optimization\nStep  E_state1          E_state2          deltaE\n')

    global atm
    inp4mecp(XNew,step)
    xyz=numpy.matrix.tolist(XNew)[0]

    global mecp_wordlist
    energylist,gradlist=[],[]
    for i in range(2):
        print(f'Run Q-Chem force calculation for state {1+i}.')
        print(f'State {1+i}: {mecp_wordlist[i]} {mecp_wordlist[2+i]}')

        inp=f'{system_name}_step{step}_state{1+i}.inp'
        out=f'{system_name}_step{step}_state{1+i}.out'
        
        tddft=False
        iroot=0
        for line in open(inp,'r').readlines():
            l=line.split()
            for i in range(len(l)):
                if 'CIS_State_Deriv' in line:
                    tddft=True
                    iroot=int(l[-1])
                    break

        if sf:
            tddft=True
            with os.popen(f'qchem -nt {nproc} {inp} {out} {scratch_path}/save_step{step}_sf') as run:
                null=run.read()
        else:
            with os.popen(f'qchem -nt {nproc} {inp} {out}') as run:
                null=run.read()

        energylist.append(energy(out,tddft,iroot))
        gradlist.append(gradient(out,tddft))

    e1,e2=energylist[0],energylist[1]
    g1,g2=gradlist[0],gradlist[1]
    ed=e2-e1

    # Backup energy.
    with open(ene_backup,'a') as backup:
        backup.write(f'{f"{step}".rjust(4)}  {f"{e1:.8f}".rjust(16)}  {f"{e2:.8f}".rjust(16)}  {f"{ed:.8f}".rjust(12)}\n')

    # Backup geometry.
    global atm
    n=len(atm)
    with open(xyz_backup,'a') as gen:
        gen.write(f'{n}\nMECP Step {step}\n')
        for i in range(n):
            gen.write(f'{atm[i]}'.ljust(2))
            for j in range(3):
                gen.write(f'{xyz[j+3*i]:.8f}'.rjust(16))
            gen.write('\n')

    os.system(f'mv *step*.inp scratch-{job_name}')
    os.system(f'mv *step*.out scratch-{job_name}')
    os.system(f'{system_name}_step{step}.geom')

    return e1,g1,e2,g2
