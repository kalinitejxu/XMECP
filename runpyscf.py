import os,sys,numpy
from parameters import *

def elemid(elem):
    id=0
    Element=['H' ,                                                                                                                                                      'He',
             'Li','Be',                                                                                                                        'B' ,'C' ,'N' ,'O' ,'F' ,'Ne',
             'Na','Mg',                                                                                                                        'Al','Si','P' ,'S' ,'Cl','Ar',
             'K' ,'Ca',                                                                      'Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
             'Rb','Sr',                                                                      'Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I' ,'Xe',
             'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',
             'Fr','Ra','Ac','Th','Pa','U' ,'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og']
    for i in range(len(Element)):
        if elem.capitalize()==Element[i]:
            id=1+i
            break
    return id

# QM template for PySCF:
# JobType CASSCF/def2SVP
# MOKIT_Key MOKIT{}
# ReadNO "guess.fch"
# State_Specific False
# Fragment False
#     frag 1.MOLA.1 "1-15 20 22"      1  2
#     frag 1.MOLB.2 "16-19 21 23-25" -1 -2
#     frag 2.MOLA.1 1-15,20,22        0  1
#     frag 2.MOLA.2 16-19,21,23-25    0  1

def keywords(template='qm.tpl'):
    # Acceptable parameters.
    runtype='CASSCF/def2SVP'
    mokitkeys='mokit{}'
    readno='guess.fch'
    state_specific=False
    fragment=False
    # Read in parameters from QM template file.
    with open(template) as read:
        for line in read.readlines():
            l=line.lower().split()
            if len(l)!=0:
                if 'jobtype' in l:
                    jobtype=l[1]
                    continue
                if 'mokit_key' in l:
                    mokitkeys=l[1]
                    continue
                if 'readno' in l:
                    readno=l[1]
                    continue
                if 'state_specific' in l:
                    if l[1]=='true':
                        state_specific=True
                        continue
                if 'fragment' in l:
                    if l[1]=='true':
                        fragment=True
                        continue
    frags=[]
    if fragment:
        frags=parse_fragment(template)
    with open('keywords_tmp.py','w') as gen:
        gen.write(f'\
runtype="{runtype}"\n\
mokitkeys="{mokitkeys}"\n\
readno="{readno}"\n\
state_specific={state_specific}\n\
')
    return runtype,mokitkeys,readno,state_specific,fragment,frags

def parse_id(inpid):
    l,ids=[],[]
    tmp=inpid.strip('"').split()
    for i in tmp:
        l.extend(i.split(','))
    for i in l:
        if '-' in i:
            a,b=i.split('-')
            for j in range(int(a),int(b)+1):
                ids.append(j)
        else:
            ids.append(int(i))
    return ids

def parse_fragment(template='qm.tpl'):
    with open(template,'r') as read:
        for line in read.readlines():
            l=line.lower.split()
            if l[0]==frag:
                frags.append([l[1],parse_id(l[2]),int(l[3]),int(l[4])])
    return frags

# Generate MOKIT input.
def inp4mokit(X0,step=1,state=1):
    global nproc,memory,system_name
    filename=f'{system_name}_mokit_step{step}_state{state}'
    with open(f'{filename}.gjf','w') as gen:
        # Job settings and keywords.
        global qm_mult,mecp_wordlist
        mult=mecp_wordlist[state-1]
        gen.write(f'%nprocshared={nproc}\n')
        gen.write(f'%mem={memory}GB\n')
        with open(qm_template,'r') as read:
            flag=False
            for line in read.readlines():
                if '#' in line and f'state{state}' in line.lower(): # e.g. state1 #p CASSCF/def2SVP guess(fragment=2)
                    gen.write(line.lower().lstrip().lstrip(f'state{state}').lstrip())
                    gen.write('\n')
                    flag=True
                    continue
                if flag:
                    if 'mokit' in line.lower():
                        if 'force' not in line.lower():
                            line=line.replace('}',',force}')
                        if 'root' not in line and ((mecp_wordlist[state-1]==qm_mult and mecp_wordlist[state+1]>=1) or\
                                                   (mecp_wordlist[state-1]!=qm_mult and mecp_wordlist[state+1]>=2) or\
                                                   (mecp_wordlist[state-1]!=qm_mult and mecp_wordlist[state+1]==1 and mecp_wordlist[state+3])):
                            line=line.replace('}',f',root={mecp_wordlist[state+1]}}}')
                            mult=qm_mult
                            if mecp_wordlist[state-1]!=qm_mult:
                                line=line.replace('}',f',Xmult={mecp_wordlist[state-1]}}}')
                        gen.write(f'{line}\n')
                        break
        # Charge and multiplicity.
        global qm_charge
        gen.write(f'{qm_charge} {mult}')
        global fragment,frags # frags: list[fragname,ids,fragchg,fragmult]; p.s. fragname: 1.MOL.1 (state.name.id)
                              # e.g. for fragment definition:
                              #      1.MOLA.1 "1-15 20 22"      1  2
                              #      1.MOLB.2 "16-19 21 23-25" -1 -2
                              #      2.MOLA.1 1-15,20,22        0  1
                              #      2.MOLA.2 16-19,21,23-25    0  1
        if fragment:
            nfrag=0
            for frag in frags:
                if state==int(frag[0].split('.')[0]):
                    nfrag+=1
            for i in range(nfrag):
                if state==int(frag[0].split('.')[0]) and 1+i==int(frag[0].split('.')[2]):
                    gen.write(f' {frag[2]} {frag[3]}')
        gen.write('\n')
        # Structure.
        global atm
        xyz=numpy.matrix.tolist(X0)[0]
        for i in range(len(atm)):
            if fragment:
                for frag in frags:
                    if state==int(frag[0].split('.')[0]):
                        if 1+i in frag[1]:
                            gen.write(f'{atm[i]}'.rjust(4)+f'(fragment={frag[0].split(".")[2]})'.ljust(20))
                            break
            else:
                gen.write(f'{atm[i]}'.ljust(4))
            for j in range(3):
                gen.write(f'{xyz[j+3*i]:.10f}'.rjust(20))
            gen.write('\n')
        gen.write('\n')
    return

# Pay attention!
# MO energies should be corresponding to occupation numbers.
def readnofch(fch):
    alpha,beta=0,0
    nocc,norb=0,0
    with open(fch,'r') as read:
        flag=False
        for line in read.readlines():
            l=line.split()
            if 'Number of alpha electrons' in line:
                alpha=int(l[-1])
                continue
            if 'Number of beta electrons' in line:
                beta=int(l[-1])
                continue
            if 'Alpha Orbital Energies' in line:
                flag=True
                continue
            if 'Alpha MO coefficients' in line:
                break
            if flag:
                for i in l:
                    if float(i)==2.0:
                        nocc+=1
                    elif float(i)==0.0:
                        break
                    else:
                        norb+=1
    nele=alpha+beta-2*nocc
    return norb,nele

def moproj(X0,fchin,fchout,norb,nele):
    # Get wavefunction type (RHF or ROHF).
    import linecache
    wfntype=linecache.getline(fchin,2)[10:14].rstrip()
    import keywords_tmp as key
    # Generate Gaussian16 input file.
    global atm,qm_charge,qm_mult
    xyz=numpy.matrix.tolist(X0)[0]
    with open('tmp.gjf','w') as gen:
        gen.write('%chk=tmp.chk\n')
        gen.write(f'#p {wfntype}/{key.runtype.split("/")[1]} 5d 7f nosymm int=nobasistransform guess(only,save)\n')
        gen.write('\n')
        gen.write('Temporary job generated by runpyscf.py\n')
        gen.write(f'\n')
        gen.write(f'{qm_charge} {qm_mult}\n')
        for i in range(len(atm)):
            gen.write(f'{atm[i]} {xyz[0+3*i]} {xyz[1+3*i]} {xyz[2+3*i]}\n')
        gen.write('\n')
    # Call Gaussian16 to generate guess wavefunction .fch format.
    with os.popen('g16 < tmp.gjf > tmp.log') as run:
        null=run.read()
    with os.popen(f'formchk tmp.chk tmp.fch') as formchk:
        null=formchk.read()
    # Natural orbital projection.
    Z=0
    for i in atm:
        Z+=elemid(i)
    nmo=norb+int((Z-qm_charge-nele)/2)
    from mokit.lib.gaussian import make_orb_resemble
    make_orb_resemble('tmp.fch',fchin,nmo)
    os.system(f'mv tmp.fch {fchout}')
    os.system(f'rm -rf fort.7 tmp.* {fchin}')
    return

def inp4pyscf(X0,step=1):
    # Get keywords.
    global qm_template
    import keywords_tmp as key
    # Get size of CI space and initial guess NO wavefunction.
    norb,nele,nofch=0,0,''
    global system_name,scratch_path
    if step==1:
        nofch=key.readno
    elif step>=2:
        nofch=f'{system_name}_step{step-1}.fch'
        os.system(f'cp {scratch_path}/{nofch} .')
    norb,nele=readnofch(nofch)
    # MO projection.
    fch=f'{system_name}_step{step}_guess.fch'
    moproj(X0,nofch,fch,norb,nele)
    # Generate PySCF input.
    with os.popen(f'bas_fch2py {fch}') as proj: # Generate PySCF input *guess.py based on *guess.fch.
        null=proj.read()
    with open(f'{system_name}_step{step}.py','w') as gen: # Add CASSCF settings.
        gen.write(f'# PySCF input created by mokit_interface.py\n')
        gen.write('from pyscf import gto,scf,mcscf,lib\n')
        gen.write('from mokit.lib.fch2py import fch2py\n')
        gen.write('from mokit.lib.py2fch import py2fch\n')
        gen.write('from mokit.lib.ortho import check_orthonormal\n')
        gen.write('from shutil import copyfile\n')
        gen.write('import numpy as np\n')
        gen.write(f'lib.num_threads({nproc})\n')
        with open(f'{system_name}_step{step}_guess.py','r') as read:
            for line in read.readlines():
                if 'import' not in line:
                    gen.write(line)
        # CASSCF job section.
        global memory
        mem=200*memory
        alpha=int((qm_mult-1+nele)/2)
        beta=nele-alpha
#       # User defined active space.
#       if addocc!=0 or addvir!=0:
#           nele+=2*addocc
#           norb+=addocc+addvir
#           alpha+=addocc
#           beta+=addocc
        if key.runtype.split('/')[0][:6].upper()=='CASSCF':
    # ---------------------------------------------------------------------------------------------------------------------------------- #
            force,mecp=False,False
            global algorithm
            if algorithm!='sp':
                force=True
            if 'mecp' in algorithm:
                mecp=True
    # -------------------------------------------- PySCF: MECP State-Average CASSCF Section -------------------------------------------- #
            if mecp:
                global mecp_wordlist
                imult,jmult,iroot,jroot=mecp_wordlist[0:4]
                if os.path.exists('pyscf_tmp.py')==False:
                    os.system(f'echo nroots={3+max(iroot,jroot)} > pyscf_tmp.py')
                from pyscf_tmp import nroots
                gen.write(f'\
\n\
### Two-State MECP Optimization ###\n\
iroot,imult={iroot},{imult}\n\
jroot,jmult={jroot},{jmult}\n\
\n\
is2=(imult**2-1)/4\n\
js2=(jmult**2-1)/4\n\
\n\
nroots={nroots}\n\
weights=[1/nroots]*nroots\n\
\n\
# Run SA-CASSCF calculation to get target states.\n\
mc=mcscf.CASSCF(mf,{norb},({alpha},{beta})).state_average_(weights)\n\
mc.fcisolver.max_memory={mem} # MB\n\
mc.max_memory={mem*4} # MB\n\
mc.max_cycle=200\n\
mc.fcisolver.max_cycle=200\n\
mc.natorb=True\n\
mc.verbose=5\n\
mc.kernel()\n\
\n\
indexi,indexj=0,1\n\
findi,findj=False,False\n\
\n\
counti,countj=0,0\n\
if imult==mf.spin_square()[1]:\n\
    counti=-1\n\
if jmult==mf.spin_square()[1]:\n\
    countj=-1\n\
\n\
for i in range(len(mc.weights)):\n\
    s2=mcscf.CASSCF(mf,{norb},({alpha},{beta})).state_specific_(i).fcisolver.spin_square(mc.ci[i],{norb},{alpha+beta})[0]\n\
    if abs(s2-is2)<1e-4:\n\
        counti+=1\n\
    if abs(s2-js2)<1e-4:\n\
        countj+=1\n\
    if findi==False and counti==iroot:\n\
        indexi=i\n\
        findi=True\n\
        print(f"Found target state 1 as CASSCF state {{i}} with <S^2>={{s2:.5f}}.")\n\
        print(f"Target state 1 CASSCF energy: {{mc.e_states[i]}} a.u.")\n\
    if findj==False and countj==jroot:\n\
        indexj=i\n\
        findj=True\n\
        print(f"Found target state 2 as CASSCF state {{i}} with <S^2>={{s2:.5f}}.")\n\
        print(f"Target state 2 CASSCF energy: {{mc.e_states[i]}} a.u.")\n\
    if findi and findj:\n\
        break\n\
\n\
numinc=0\n\
while (findi and findj)==False:\n\
    numinc+=1\n\
    print(f"Fail to find target state(s) in current {{nroots}} roots. Try increase to {{nroots+3*numinc}} roots.")\n\
\n\
    nroots+=3*numinc\n\
    weights=[1/nroots]*nroots\n\
\n\
    mc=mcscf.CASSCF(mf,{norb},({alpha},{beta})).state_average_(weights)\n\
    mc.fcisolver.max_memory={mem} # MB\n\
    mc.max_memory={mem*4} # MB\n\
    mc.max_cycle=200\n\
    mc.fcisolver.max_cycle=200\n\
    mc.natorb=True\n\
    mc.verbose=5\n\
    mc.kernel()\n\
\n\
    indexi,indexj=0,1\n\
    findi,findj=False,False\n\
\n\
    counti,countj=0,0\n\
    if imult==mf.spin_square()[1]:\n\
        counti=-1\n\
    if jmult==mf.spin_square()[1]:\n\
        countj=-1\n\
\n\
    for i in range(len(mc.weights)):\n\
        s2=mcscf.CASSCF(mf,{norb},({alpha},{beta})).state_specific_(i).fcisolver.spin_square(mc.ci[i],{norb},{alpha+beta})[0]\n\
        if abs(s2-is2)<1e-4:\n\
            counti+=1\n\
        if abs(s2-js2)<1e-4:\n\
            countj+=1\n\
        if findi==False and counti==iroot:\n\
            indexi=i\n\
            findi=True\n\
            print(f"Found target state 1 as CASSCF state {{i}} with <S^2>={{s2:.5f}}.")\n\
            print(f"Target state 1 CASSCF energy: {{mc.e_states[i]}} a.u.")\n\
        if findj==False and countj==jroot:\n\
            indexj=i\n\
            findj=True\n\
            print(f"Found target state 2 as CASSCF state {{i}} with <S^2>={{s2:.5f}}.")\n\
            print(f"Target state 2 CASSCF energy: {{mc.e_states[i]}} a.u.")\n\
        if findi and findj:\n\
            break\n\
\n\
    if findi and findj:\n\
        import os\n\
        os.system(f"echo nroots={{nroots}} > pyscf_tmp.py")\n\
        break\n\
\n')
    # ---------------------------------------------------------------------------------------------------------------------------------- #
            else:
                imult,iroot=mecp_wordlist[0],mecp_wordlist[2]
                if os.path.exists('pyscf_tmp.py')==False:
                    os.system(f'echo nroots={3+iroot} > pyscf_tmp.py')
                from pyscf_tmp import nroots
    # ---------------------------------------------- PySCF: State-Specific CASSCF Section ---------------------------------------------- #
                if iroot!=0 and key.state_specific:
                    gen.write(f'\
\n\
imult,iroot={imult},{iroot}\n\
nroots={nroots}\n\
\n\
for i in range(250):\n\
    print("ITER=",1+i)\n\
    mc=mcscf.CASCI(mf,{norb},({alpha},{beta}))\n\
    mc.max_memory={4*mem} # MB\n\
    mc.fcisolver.max_memory={mem} # MB\n\
    mc.fcisolver.nroots=nroots\n\
    mc.fcisolver.spin={mult-1}\n\
    mc.fcisolver.max_cycle=200\n\
    mc.verbose=5\n\
    mc.kernel()\n\
\n\
    counti=0\n\
    if imult==mf.spin_square()[1]:\n\
        counti=-1\n\
\n\
    for j in range(nroots):\n\
        ss=mc.fcisolver.spin_square(mc.ci[j],mc.ncas,mc.nelecas)\n\
        if (abs(ss[0]-{(mult**2-1)/4})<1e-4):\n\
            counti+=1\n\
        if counti==iroot:\n\
            break\n\
    e=mc.e_tot[j]\n\
    nroots=3+j\n\
\n\
    mc=mcscf.CASSCF(mf,{norb},({alpha},{beta})).state_specific_(j)\n\
    mc.max_memory={4*mem} # MB\n\
    mc.fcisolver.max_memory={mem} # MB\n\
    mc.fcisolver.nroots=nroots\n\
    mc.fcisolver.spin={mult-1}\n\
    mc.fcisolver.max_cycle=200\n\
    mc.verbose=5\n\
    mc.max_cycle=1\n\
    mc.kernel()\n\
    mf.mo_coeff=mc.mo_coeff.copy()\n\
    if (abs(e-mc.e_tot)<1e-7):\n\
        print(f"State-Specific Tracking Iteration {{1+i}}: {{e}} {{mc.e_tot}}")\n\
        print(f"Energy Diff.: {{abs(e-mc.e_tot)}}")\n\
        print("Success in State-Specific Tracking.")\n\
        break\n\
\n\
mc=mcscf.CASSCF(mf,{norb},({alpha},{beta})).state_specific_(j)\n\
mc.fcisolver.max_memory={mem} # MB\n\
mc.max_memory={mem*4} # MB\n\
mc.fcisolver.nroots=nroots\n\
mc.fcisolver.spin={mult-1}\n\
mc.fcisolver.max_cycle=200\n\
mc.natorb=True\n\
mc.verbose=5\n\
mc.kernel()\n\
\n')
    # ---------------------------------------------------------------------------------------------------------------------------------- #
    # ---------------------------------------------- PySCF: State-Average CASSCF Section ----------------------------------------------- #
                elif iroot!=0:
                    gen.write(f'\
\n\
iroot,imult={iroot},{imult}\n\
is2=(imult**2-1)/4\n\
\n\
nroots={nroots}\n\
weights=[1/nroots]*nroots\n\
\n\
mc=mcscf.CASSCF(mf,{norb},({alpha},{beta})).state_average_(weights)\n\
mc.fcisolver.max_memory={mem} # MB\n\
mc.max_memory={mem*4} # MB\n\
mc.max_cycle=200\n\
mc.fcisolver.max_cycle=200\n\
mc.natorb=True\n\
mc.verbose=5\n\
mc.kernel()\n\
\n\
find=False\n\
counter=0\n\
if imult==mf.spin_square()[1]:\n\
    counter=-1\n\
\n\
for i in range(len(mc.weights)):\n\
    ss=mcscf.CASSCF(mf,{norb},({alpha},{beta})).state_specific_(i).fcisolver.spin_square(mc.ci[i],{norb},{alpha+beta})[0]\n\
    if abs(ss-is2)<1e-4:\n\
        counter+=1\n\
    if counter==iroot:\n\
        find=True\n\
        print(f"Found target state as CASSCF state {{i}} with <S^2>={{ss:.5f}}.")\n\
        print(f"Target root CASSCF energy: {{mc.e_states[i]}} a.u.")\n\
        break\n\
\n\
numinc=0\n\
while find==False:\n\
    numinc+=1\n\
    print(f"Fail to find target state(s) in current {{nroots}} roots. Try increase to {{nroots+3*numinc}} roots.")\n\
\n\
    nroots+=3*numinc\n\
    weights=[1/nroots]*nroots\n\
\n\
    mc=mcscf.CASSCF(mf,{norb},({alpha},{beta})).state_average_(weights)\n\
    mc.fcisolver.max_memory={mem} # MB\n\
    mc.max_memory={mem*4} # MB\n\
    mc.max_cycle=200\n\
    mc.fcisolver.max_cycle=200\n\
    mc.natorb=True\n\
    mc.verbose=5\n\
    mc.kernel()\n\
\n\
    counter=0\n\
    if imult==mf.spin_square()[1]:\n\
        counter=-1\n\
    for i in range(len(mc.weights)):\n\
        ss=mcscf.CASSCF(mf,{norb},({alpha},{beta})).state_specific_(i).fcisolver.spin_square(mc.ci[i],{norb},{alpha+beta})[0]\n\
        if abs(ss-is2)<1e-4:\n\
            counter+=1\n\
        if counter==iroot:\n\
            find=True\n\
            print(f"Found target state as CASSCF state {{i}} with <S^2>={{ss:.5f}}.")\n\
            print(f"Target root CASSCF energy: {{mc.e_states[i]}} a.u.")\n\
            break\n\
    if find:\n\
        import os\n\
        os.system(f"echo nroots={{nroots}} > pyscf_tmp.py")\n\
        break\n\
\n')
    # ---------------------------------------------------------------------------------------------------------------------------------- #
    # ----------------------------------------------- PySCF: Ground State CASSCF Section ----------------------------------------------- #
                elif iroot==0:
                    gen.write(f'\
\n\
mc=mcscf.CASSCF(mf,{norb},({alpha},{beta}))\n\
mc.fcisolver.max_memory={mem} # MB\n\
mc.max_memory={mem*4} # MB\n\
mc.max_cycle=200\n\
mc.fcisolver.spin={mult-1}\n\
mc.fcisolver.max_cycle=200\n\
mc.natorb=True\n\
mc.verbose=5\n\
mc.kernel()\n\
\n')
    # ---------------------------------------------------------------------------------------------------------------------------------- #
            gen.write('# save NOs into .fch file\n')
            gen.write(f'copyfile("{fch}", "{system_name}_step{step}.fch")\n')
            gen.write(f"py2fch('{system_name}_step{step}.fch',nbf,nif,mc.mo_coeff,'a',mc.mo_occ,True)\n")
    # --------------------------------------------------- CASSCF Analytical Gradient --------------------------------------------------- #
            if force:
                gen.write('from pyscf import grad\n')
                if mecp:                                        # MECP state-average CASSCF gradient.
                    gen.write(f'\
\n\
print("State-average CASSCF gradient for MECP state 1:")\n\
mcg=grad.sacasscf.Gradients(mc,state=indexi)\n\
mcg.max_memory={mem*5} # MB\n\
mcg.kernel()\n\
\n\
print("State-average CASSCF gradient for MECP state 2:")\n\
mcg=grad.sacasscf.Gradients(mc,state=indexj)\n\
mcg.max_memory={mem*5} # MB\n\
mcg.kernel()\n\
\n')
                elif iroot!=0 and key.state_specific==False:        # State-average CASSCF gradient.
                    gen.write(f'\
\n\
mcg=grad.sacasscf.Gradients(mc,state=i)\n\
mcg.max_memory={mem*5} # MB\n\
mcg.kernel()\n\
\n')
                else:                                           # State-specific CASSCF gradient.
                    gen.write(f'\
\n\
mcg=mc.Gradients()\n\
mcg.max_memory={mem*5} # MB\n\
mcg.kernel()\n\
\n')
    os.system(f'rm {system_name}_step{step}_guess.py')
    return

def pyscf_results(out,mecp=False):
    flag=False
    e1,e2=0.,0.
    g1,g2=[],[]
    with open(out,'r') as read:
        flag=False
        if mecp:
            for line in read.readlines():
                l=line.split()
                if 'Target state 1 CASSCF energy:' in line:
                    e1=float(l[-2])
                    continue
                if 'Target state 2 CASSCF energy:' in line:
                    e2=float(l[-2])
                    continue
                if 'State-average CASSCF gradient for MECP state 1:' in line:
                    flag=True
                    continue
                if 'State-average CASSCF gradient for MECP state 2:' in line:
                    flag=True
                    continue
                if flag and len(g1)<3*natms:
                    if len(l)==5:
                        for i in range(3):
                            g1.append(float(l[2+i]))
                    if len(g1)==3*natms:
                        flag=False
                elif flag:
                    if len(l)==5:
                        for i in range(3):
                            g2.append(float(l[2+i]))
                    if len(g2)==3*natms:
                        break
            return e1,g1,e2,g2
        else:
            import keywords_tmp as key
            for line in read.readlines():
                l=line.split()
                if iroot==0 or (iroot!=0 and key.state_specific):
                    if 'CASSCF energy' in line:
                        e1=float(l[-2])
                        continue
                    elif 'CASSCF gradients' in line:
                        flag=True
                        continue
                else:
                    if 'Target root CASSCF energy:' in line:
                        e1=float(l[-2])
                        continue
                    elif 'StateAverageMCSCF gradients' in line:
                        flag=True
                        continue
                if flag:
                    if len(l)==5:
                        for i in range(3):
                            g1.append(float(l[2+i]))
                if len(g1)==3*natms:
                    break
            return e1,g1

# Get electrostatic interaction energy that is to be excluded from final energy of MOKIT output.
def bqene(mokitout):
    with open(mokitout,'r') as read:
        for line in read.readlines():
            if 'Coulomb interaction energy of background point charges:' in line:
                return float(line.split()[-2])
    return 0.

def mokit_energy(out):
    ene=0.
    with open(out,'r') as read:
        for line in read.readlines():
            l=line.split()
            if 'E(CASSCF)' in line:
                ene=float(l[2])
                ene-=bqene(out) # Extract background charges electrostatic interaction energy.
                return ene

def mokit_gradient(out):
    grad=[]
    with open(out,'r') as read:
        flag=False
        for line in read.readlines():
            l=line.split()
            if 'Cartesian gradient (HARTREE/BOHR):' in line:
                flag=True
                continue
            elif 'Leave subroutine do_cas' in line:
                flag=False
            if flag:
                for i in range(len(l)):
                    grad.append(float(l[i]))
        if len(grad)!=3*len(atm):
            print('Dimension of cartesian gradient != 3*Natoms, please check!')
            sys.exit()
    return grad

def optstep(X0,step=1):
    return

def mecpstep(X0,step=1):
    print(f'Optimization steps: {step}')

    global scratch_path,xyz_backup,ene_backup
    if step==1:
        os.system(f'rm -rf {scratch_path} {ene_backup} {xyz_backup}')
        os.system(f'mkdir {scratch_path}')
        with open(ene_backup,'w') as gen:
            gen.write('MECP Optimization\nStep  E_state1          E_state2          deltaE\n')

    global mecp_wordlist
    e1,e2,g1,g2=0.,0.,[],[]

    import keywords_tmp as key
    if step==1 and key.readno=='':
        inp4mokit(X0,step,1+state)
        print(f'Run MOKIT AutoMR job for state {1+state}.')
        print(f'State {1+state}: {mecp_wordlist[state]} {mecp_wordlist[2+state]}')
        # Run MOKIT AutoMR job.
        filename=f'{system_name}_mokit_step{step}_state{1+state}'
        with os.popen(f'automr {filename}.gjf > {filename}.out') as run:
            null=run.read()
        energylist.append(mokit_energy(f'{filename}.out'))
        gradlist.append(mokit_gradient(f'{filename}.out'))
        # Backup files.
        os.system(f'mv {filename}.gjf {scratch_path}')
        os.system(f'mv {filename}.out {scratch_path}')
        os.system(f'mv {filename}*CASSCF_NO.fch {scratch_path}/{system_name}_step{step}_state{1+state}_CASSCF_NO.fch')
        os.system('rm *_mokit_*')

    else:
        inp4pyscf(X0,step)
        print(f'Run PySCF SA-CASSCF job.')
        # Run PySCF job.
        filename=f'{system_name}_step{step}'
        with os.popen(f'python -u {filename}.py > {filename}.out') as run:
            null=run.read()
        e1,g1,e2,g2=pyscf_results(f'{filename}.out',mecp=True)
        # Backup files.
        os.system('rm *guess*')
        os.system(f'mv {filename}* {scratch_path}')

    ed=e2-e1

    # Backup energy.
    with open(ene_backup,'a') as backup:
        backup.write(f'{f"{step}".rjust(4)}  {f"{e1:.8f}".rjust(16)}  {f"{e2:.8f}".rjust(16)}  {f"{ed:.8f}".rjust(12)}\n')

    # Backup geometry.
    global atm
    n=len(atm)
    xyz=numpy.matrix.tolist(X0)[0]
    with open(xyz_backup,'a') as gen:
        gen.write(f'{n}\nMECP Step {step}\n')
        for i in range(n):
            gen.write(f'{atm[i]}'.ljust(2))
            for j in range(3):
                gen.write(f'{xyz[j+3*i]:.8f}'.rjust(16))
            gen.write('\n')

    return e1,g1,e2,g2
