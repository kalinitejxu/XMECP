import os,numpy,geom
from parameters import *

### >>> Special for ORCA spin-flip framework >>> ###
def spin_flip(xyzfile,step,incnroots=0):
    global job_name,system_name,nproc,qm_charge,qm_mult,qm_template,qm_path
    global mecp_wordlist
    global sf_ref_mult

    print(f'Spin-flip (flip-down) is turned on. Multiplicity of reference state is {sf_ref_mult}.')
    # Generate a spin-flip TD-DFT single point task to decide iroot and jroot.
    inp=f'{system_name}_step{step}_sf.inp'
    out=f'{system_name}_step{step}_sf.out'

    with open(inp,'w') as gen:
        gen.write('# ORCA input file generated by XMECP\n')
        gen.write(open(qm_template,'r').readlines()[0])
        gen.write(f'\n%pal nproc {nproc} end')

        gen.write('\n%tddft')
        gen.write('\n  sf true')
        nroots=3+5*incnroots+max(state1_root,state2_root)
        print(f'Search for required states in the lowest {nroots} states.')
        gen.write(f'\n  nroots {nroots}')
        gen.write('\nend')

        gen.write('\n')
        flag=False
        extrakeys=open(qm_template,'r')
        for line in extrakeys.readlines():
            L=line.split()
            if len(L)!=0:
                probe=line[0]
                if probe[0]=='%':
                    flag=True
                if flag:
                    gen.write('\n'+line)
                if L[-1]=='end':
                    flag=False
        extrakeys.close()
        gen.write(f'\n* xyzfile {qm_charge} {sf_ref_mult} {xyzfile}\n')

    # Now run spin-flip TD-DFT single point task.
    with os.popen(f'{qm_path} {inp} > {out}') as run:
        null=run.read()

    # Search for required states.
    flag=False
    s2val,ijroots=[],[]
    
    with open(out,'r') as read:
        for line in read.readlines():
            L=line.split()
            if len(L)!=0:
                if (L[0]=='SF-CIS' or L[0]=='SF-TDA') and L[1]=='EXCITED' and L[2]=='STATES':
                    flag=True
                    continue
                if flag:
                    if L[0]=='STATE' and L[-3]=='<S**2>':
                        s2val.append(float(L[-1]))
                        if len(s2val)==nroots:
                            break
    flags=[False,False]
    for i in range(2):
        count=0
        s=(mecp_wordlist[0+i]-1)/2
        s2=s*(s+1)
        for j in range(len(s2val)):
            if s2-0.8 <= s2val[j] <= s2+1.2:
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
    global sf
    xyz=numpy.matrix.tolist(X0)[0]
    xyzfile=system_name+'_step'+str(step)+'.xyz'
    geom.writexyz(atm,xyz,xyzfile,'step'+str(step))

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

    for i in range(2):
        # ORCA input files are generated in current folder and filenames are set to ${system_name}_stepN_state1.gjf and ${system_name}_stepN_state2.gjf.
        filename=f'{system_name}_step{step}_state{1+i}.inp'

        with open(filename,'w') as gen:
            gen.write('# ORCA input file generated by XMECP\n')
            gen.write(open(qm_template,'r').readlines()[0].strip('\n'))
            gen.write(f' engrad\n%pal nproc {nproc} end')

            if sf or step>1:
                gen.write('\n%scf guess moread end')
                if sf:
                    gen.write(f'\n%moinp "{system_name}_step{step}_sf.gbw"')
                else:
                    gen.write(f'\n%moinp "{system_name}_step{step}_state{i+1}.gbw"')

            global mecp_wordlist
            tddft=False
            if sf or\
               (mecp_wordlist[0+i]==qm_mult and mecp_wordlist[2+i]!=0) or\
               (mecp_wordlist[0+i]!=qm_mult and mecp_wordlist[2+i]==1 and mecp_wordlist[4+i]==True) or\
               (mecp_wordlist[0+i]!=qm_mult and mecp_wordlist[2+i]>=2):
                tddft=True
                gen.write('\n%tddft')
                iroot=mecp_wordlist[2+i]
                if sf:
                    gen.write('\n  sf true')
                    iroot=ijroots[i]
                else:
                    if qm_mult==1 and mecp_wordlist[0+i]==3:
                        gen.write('\n  triplets true) ')
                        gen.write('\n  tda false')
                gen.write(f'\n  iroot {iroot}')
                gen.write(f'\n  nroots {3+iroot}')
                gen.write('\nend')

            gen.write('\n')
            extrakeys=open(qm_template,'r')
            for line in extrakeys.readlines()[1:]:
                L=line.split()
                if len(L)!=0:
                    gen.write(line)
                else:
                    break
            extrakeys.close()

            if sf:
                gen.write(f'\n* xyzfile {qm_charge} {sf_ref_mult} {xyzfile}\n')
            else:
                gen.write(f'\n* xyzfile {qm_charge} {qm_mult} {xyzfile}\n')
    return
### <<< Generate ORCA input file for two-states MECP optimization <<< ###

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

def energy(energyfile):
    energy=0.0
    flag=False
    readenergy=open(energyfile,'r')
    for line in readenergy.readlines():
        L=line.split()
        if len(L)!=0:
            if len(L)>4:
                if L[1]=='The' and L[2]=='current' and L[3]=='total' and L[4]=='energy':
                    flag=True
            if flag and is_number(L[0]):
                energy=float(L[0])
                break
    readenergy.close()
    return energy

def gradient(gradfile):
    natms,grad=0,[]
    flagatm,flag=False,False
    readgrad=open(gradfile,'r')
    for line in readgrad.readlines():
        L=line.split()
        if len(L)!=0:
            if len(L)>3:
                if L[1]=='Number' and L[2]=='of' and L[3]=='atoms':
                    flagatm=True
                    continue
                if L[1]=='The' and L[2]=='current' and L[3]=='gradient':
                    flag=True
                    continue
            if flagatm and is_number(L[0]):
                natms=int(L[0])
                flagatm=False
                continue
            if flag and is_number(L[0]):
                grad.append(float(L[0])/0.52917721067)
                continue
            if flag and len(grad)==3*natms:
                break
    readgrad.close()
    return grad

def mecpstep(XNew,step):
    print(f'Optimization steps: {step}')
    
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
        print(f'Run ORCA engrad calculation for state {1+i}.')
        print(f'State {1+i}: {mecp_wordlist[i]} {mecp_wordlist[2+i]}')

        inp=f'{system_name}_step{step}_state{1+i}.inp'
        out=f'{system_name}_step{step}_state{1+i}.out'
        engrad=f'{system_name}_step{step}_state{1+i}.engrad'
        
        with os.popen(f'{qm_path} {inp} > {out}') as run:
            null=run.read()

        energylist.append(energy(engrad))
        gradlist.append(gradient(engrad))
    
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
    os.system('rm -f *.cis *.densities *_property.txt')

    for i in range(2):
        gbw=f'{system_name}_step{step}_state{1+i}.gbw'
        nextgbw=f'{system_name}_step{step+1}_state{1+i}.gbw'
        os.system(f'mv {gbw} {nextgbw}')
    
    return e1,g1,e2,g2
