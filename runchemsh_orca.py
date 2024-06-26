# 2023-03-14 updated
# External electrostatic field.

# Combined with ORCA-ChemShell interface in ChemShell package.
# 2023-02-12
# xujiawei@fjirsm.ac.cn

import os,sys,numpy,geom
from parameters import *

# Region specify format in ChemShell template:
# set qm_atom_numbers { 1-10 20-30 35 }
# set active_atom_numbers { 1-1000 2000-3000 }

### >>> Generate ChemShell input file >>> ###
# Spin-flip single point check for MECP optimization.
def spin_flip(step,incnroots=0):
    global system_name,nproc,qm_charge,qm_mult,qm_template
    global mecp_wordlist
    global sf_ref_mult
    global atm
    
    print(f'Spin-flip (flip-down) is turned on. Multiplicity of reference state is {sf_ref_mult}.')
    # Generate a spin-flip TD-DFT single point task to decide iroot and jroot.
    chm=f'{system_name}_step{step}_sf.chm'
    out=f'{system_name}_step{step}_sf.out'

    with open(chm,'w') as gen:
        gen.write(f'# Tcl-ChemShell input file generated by XMECP\n')
        gen.write(f'# Spin-flip TD-DFT single point calculation. Step {step}.\n')

        # Write basci settings.
        gen.write(f'set sysname {system_name}\n')
        gen.write(f'set qmcharge {qm_charge}\n')
        gen.write(f'set qmmult {sf_ref_mult}\n')

        # Read ORCA simple input for QM template file.
        gen.write('\nset orcasimpleinput "')
        with open(qm_template,'r') as qmtpl:
            gen.write(qmtpl.readlines()[0].rstrip('\n'))

        if sf_ref_mult!=1 or (sf_ref_mult==1 and uks):
            gen.write(' uks')
        gen.write('"\n')

        # Write ORCA blocks.
        gen.write(f'set orcablocks {{\n%pal nprocs {nproc} end\n')
        if step>1:
            gen.write('%scf\n  guess moread\nend\n')
            gen.write(f'%moinp "./{scratch_path}/{system_name}_orca_step{step-1}_sf.gbw"\n')

        nroots=3+5*incnroots+max(state1_root,state2_root)
        # Write TD-DFT block.
        gen.write('%tddft sf true\n')
        gen.write(f'  nroots {nroots}\n')
        gen.write('end\n')

        # Write other blocks required by QM template file.
        with open(qm_template,'r') as qmtpl:
            flag=False
            for line in qmtpl.readlines():
                if len(line)!=0:
                    probe=line[0]
                    if probe[0]=='%':
                        flag=True
                    if flag:
                        gen.write('\n'+line)
                    if line[-1]=='end':
                        flag=False
        gen.write('}\n\n')

        # Write ORCA settings.
        with open(chemsh_template,'r') as tpl:
            for line in tpl.readlines():
                if 'gradient' in line and '${sysname}' in line:
                    continue
                if 'eandg' in line:
                    line=line.replace('eandg','energy')
                gen.write(line)

    # Now run spin-flip TD-DFT single point task.
    with os.popen(f'chemsh < {chm} > {out}') as run:
        null=run.read()

    orca_inp=f'{system_name}_orca.inp'
    orca_out=f'{system_name}_orca.out'
    orca_gbw=f'{system_name}_orca.gbw'

    # Search for required states.
    flag=False
    s2val,ijroots=[],[]
    with open(orca_out,'r') as sfsp:
        for line in sfsp.readlines():
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
                        continue
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
    
    # Delete useless files.
    os.system(f'rm -f {system_name}.e FIELD CONFIG CONTROL hybrid* pointcharges.xyz')
    os.system(f'rm -f {system_name}_orca.densities {system_name}_orca_property.txt')
    # Backup ChemShell files
    os.system(f'mv {chm} {scratch_path}/{chm}')
    os.system(f'mv {out} {scratch_path}/{out}')
    # Backup ORCA files.
    if os.path.exists(orca_inp):
        os.system(f'mv {orca_inp} {scratch_path}/{system_name}_orca_step{step}_sf.inp')
    if os.path.exists(orca_out):
        os.system(f'mv {orca_out} {scratch_path}/{system_name}_orca_step{step}_sf.out')
    if os.path.exists(orca_gbw):
        os.system(f'mv {orca_gbw} {scratch_path}/{system_name}_orca_step{step}_sf.gbw')

    return ijroots

# Note for ChemShell template file:
# ORCA keyword "UKS" and blocks %TDDFT, %MOInp, %Pal and %SCF(guess) will be automatically created.
def inp4mecp(X0,step=1):
    global system_name,nproc,qm_charge,qm_mult,uks,chemsh_template
    global atm

    # Write coordinates input.
    xyzfile=f'{system_name}.c'
    prmfile=f'{system_name}.prmtop'
    geom.chemshc(atm,X0,prmfile,xyzfile)

    ijroots=[]
    incnroots=0

    if sf:
        while len(ijroots)!=2:
            ijroots=spin_flip(step,incnroots)
            if len(ijroots)!=2:
                incnroots+=1
                print('Increase CIS space for another 5 states.')
            if incnroots==5:
                print('Target state was not found with 20 more states expanding the CIS space.')
                print('Your system may not be suitable for spin-flip calculation. Exit program!')
                sys.exit()

    # Write force inputs.
    global mecp_wordlist
    for state in range(2):
        chm=f'{system_name}_step{step}_state{1+state}.chm'
        with open(chm,'w') as gen:
            gen.write(f'# Tcl-ChemShell input file generated by XMECP\n')
            gen.write(f'# Searching for two-state MECP: step {step} / state {1+state}\n')

            tddft=False
            if sf or\
               (mecp_wordlist[0+state]==qm_mult and mecp_wordlist[2+state]!=0) or\
               (mecp_wordlist[0+state]!=qm_mult and mecp_wordlist[2+state]==1 and mecp_wordlist[4+state]==True) or\
               (mecp_wordlist[0+state]!=qm_mult and mecp_wordlist[2+state]>=2):
                tddft=True

            statemult=mecp_wordlist[state]
            if sf:
                statemult=sf_ref_mult
            elif tddft:
                statemult=qm_mult

            # Write basci settings.
            gen.write(f'set sysname {system_name}\n')
            gen.write(f'set qmcharge {qm_charge}\n')
            gen.write(f'set qmmult {statemult}\n')

            # Read ORCA simple input for QM template file.
            gen.write('\nset orcasimpleinput "')
            with open(qm_template,'r') as qmtpl:
                gen.write(qmtpl.readlines()[0].rstrip('\n'))

            if (tddft and qm_mult!=1) or\
               (tddft==False and mecp_wordlist[0+state]!=1) or\
               (tddft==False and mecp_wordlist[0+state]==1 and uks):
                gen.write(' uks')
            gen.write('"\n')

            # Write ORCA blocks.
            gen.write(f'set orcablocks {{\n%pal nprocs {nproc} end\n')

            global esf
            Ex,Ey,Ez=0.,0.,0.
            # Add external electrostatic field.
            if esf:
                global chg,fielddef
                save_stdout=sys.stdout
                with open('trash.txt','w') as sys.stdout:
                    from electrostatic_field import electrostatic_field
                    Ex,Ey,Ez=electrostatic_field(X0,chg,fielddef,False)
                sys.stdout=save_stdout
                gen.write(f'%scf\n  efield {Ex:.6f},{Ey:.6f},{Ez:.6f}\n  end\n')
            if step>1:
                gen.write('%scf\n  guess moread\n  end\n')
                if sf:
                    gen.write(f'%moinp "./{scratch_path}/{system_name}_orca_step{step-1}_sf.gbw"\n')
                else:
                    gen.write(f'%moinp "./{scratch_path}/{system_name}_orca_step{step-1}_state{1+state}.gbw"\n')

            # Write TD-DFT block.
            if tddft:
                gen.write('%tddft\n')
                iroot=mecp_wordlist[2+state]
                if sf:
                    iroot=ijroots[state]
                else:
                    if mecp_wordlist[0+state]==1:
                        gen.write('  triplets false\n')                
                gen.write(f'  iroot {iroot}\n')
                gen.write(f'  nroots {3+iroot}\n')
                gen.write('end\n')
            
            # Write other blocks required by QM template file.
            with open(qm_template,'r') as qmtpl:
                flag=False
                for line in qmtpl.readlines():
                    if len(line)!=0:
                        probe=line[0]
                        if probe[0]=='%':
                            flag=True
                        if flag:
                            gen.write(line)
                        if line[-1]=='end':
                            flag=False
            gen.write('}\n\n')
        os.system(f'cat {chemsh_template} >> {chm}')
    return

def inp4opt(X0,step=1):
    global system_name,nproc,qm_charge,qm_mult,chemsh_template
    global atm

    # Write coordinates input.
    xyzfile=f'{system_name}.c'
    prmfile=f'{system_name}.prmtop'
    geom.chemshc(atm,X0,prmfile,xyzfile)

    chm=f'{system_name}_step{step}.chm'
    with open(chm,'w') as gen:
        gen.write(f'# Tcl-ChemShell input file generated by XMECP\n')
        gen.write(f'# Optimization: step {step}\n')

        # Write basci settings.
        gen.write(f'set sysname {system_name}\n')
        gen.write(f'set qmcharge {qm_charge}\n')
        gen.write(f'set qmmult {qm_mult}\n')

        # Read ORCA simple input for QM template file.
        gen.write('\nset orcasimpleinput "')
        with open(qm_template,'r') as qmtpl:
            gen.write(qmtpl.readlines()[0].rstrip('\n'))
        gen.write('"\n')

        # Write ORCA blocks.
        gen.write(f'set orcablocks {{\n%pal nprocs {nproc} end\n')
        
        global esf
        Ex,Ey,Ez=0.,0.,0.
        # Add external electrostatic field.
        if esf:
            global chg,fielddef
            save_stdout=sys.stdout
            with open('trash.txt','w') as sys.stdout:
                from electrostatic_field import electrostatic_field
                Ex,Ey,Ez=electrostatic_field(X0,chg,fielddef,False)
            sys.stdout=save_stdout
            gen.write(f'%scf\n  efield {Ex:.6f},{Ey:.6f},{Ez:.6f}\n  end\n')

        if step>1:
            gen.write('%scf\n  guess moread\n  end\n')
            gen.write(f'%moinp "./{scratch_path}/{system_name}_orca_step{step-1}.gbw"\n')

        # Write other blocks required by QM template file.
        with open(qm_template,'r') as qmtpl:
            flag=False
            for line in qmtpl.readlines():
                if len(line)!=0:
                    probe=line[0]
                    if probe[0]=='%':
                        flag=True
                    if flag:
                        gen.write(line)
                    if line[-1]=='end':
                        flag=False
        gen.write('}\n\n')
    os.system(f'cat {chemsh_template} >> {chm}')
    return
### <<< Generate ChemShell input file <<< ###

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

### >>> Read energies and gradients >>> ###
def energy(energyfile):
    flag=False
    ene=0.0
    with open(energyfile,'r') as readene:
        for line in readene.readlines():
            L=line.split()
            if flag:
                ene=float(L[0])
            if len(L)==4 and L[2]=='dimensions=1':
                flag=True
    return ene
def gradient(gradfile):
    flag=False
    grad=[]
    readgrad=open(gradfile,'r')
    for line in readgrad.readlines():
        L=line.split()
        if flag:
            grad.append(float(L[0]))
        if len(L)==4 and L[2]=='dimensions=3':
            flag=True
    return grad
### <<< Read energies and gradients <<< ###

def mecpstep(XNew,step):
    print(f'Optimization steps: {step}')

    global scratch_path,xyz_backup,ene_backup
    if step==1:
        os.system(f'rm -rf {xyz_backup} {ene_backup} {scratch_path}; mkdir {scratch_path}')
        with open(ene_backup,'w') as gen:
            gen.write('MECP Optimization\nStep  E_state1          E_state2          deltaE\n')

    xyz=numpy.matrix.tolist(XNew)[0]
    # Recover coordinates of the whole QM/MM system.
    xyz=geom.upd_active_xyz(xyz)
    inp4mecp(xyz,step)

    global system_name,mecp_wordlist
    enelist,gradlist=[],[]
    for state in range(2):
        print(f'Run ChemShell eandg calculation for state {1+state}.')
        print(f'State {1+state}: {mecp_wordlist[state]} {mecp_wordlist[2+state]}')

        chm=f'{system_name}_step{step}_state{1+state}.chm'
        out=f'{system_name}_step{step}_state{1+state}.out'

        os.system(f'echo "set step {step} # Note: for PySCF interface" > xmecp.tmp\n')
        os.system(f'echo "set state {1+state} # Note: for PySCF interface" >> xmecp.tmp\n')
        os.system(f'echo "set target_mult {mecp_wordlist[0+state]} # Note: for PySCF SA-CASSCF interface" >> xmecp.tmp\n')
        os.system(f'echo "set target_root {mecp_wordlist[2+state]} # Note: for PySCF SA-CASSCF interface" >> xmecp.tmp\n')

        with os.popen(f'chemsh < {chm} > {out}') as run:
            null=run.read()

        e=f'{system_name}.e'
        g=f'{system_name}.g'

        orca_inp=f'{system_name}_orca.inp'
        orca_out=f'{system_name}_orca.out'
        orca_gbw=f'{system_name}_orca.gbw'

        enelist.append(energy(e))
        gradlist.append(gradient(g))

        # Delete useless files.
        os.system(f'rm -f {system_name}.e {system_name}.g FIELD CONFIG CONTROL hybrid* pointcharges.xyz')
        os.system(f'rm -f {system_name}_orca.pcgrad {system_name}_orca.engrad {system_name}_orca.cis {system_name}_orca.densities {system_name}_orca_property.txt')
        os.system(f'rm -f {system_name}_step{step}_guess_NO.fch')

        # Backup ChemShell files
        os.system(f'mv {chm} {scratch_path}/{system_name}_step{step}_state{1+state}.chm')
        os.system(f'mv {out} {scratch_path}/{system_name}_step{step}_state{1+state}.out')
        # Backup ORCA files. (if there are)
        if os.path.exists(orca_inp):
            os.system(f'mv {orca_inp} {scratch_path}/{system_name}_orca_step{step}_state{1+state}.inp')
        if os.path.exists(orca_out):
            os.system(f'mv {orca_out} {scratch_path}/{system_name}_orca_step{step}_state{1+state}.out')
        if os.path.exists(orca_gbw):
            os.system(f'mv {orca_gbw} {scratch_path}/{system_name}_orca_step{step}_state{1+state}.gbw')

        # Backup MOKIT and PySCF files. (if there are)
        mokit_gjf=f'{system_name}_mokit.gjf'
        mokit_out=f'{system_name}_mokit.out'
        pyscf_inp=f'{system_name}_step{step}_state1_CASSCF.py'
        pyscf_out=f'{system_name}_step{step}_state1_CASSCF.out'
        if step==1 and os.path.exists(mokit_gjf):
            os.system(f'mv {mokit_gjf} {scratch_path}/{system_name}_mokit_step{step}_state{1+state}.py')
        if step==1 and os.path.exists(mokit_out):
            os.system(f'mv {mokit_out} {scratch_path}/{system_name}_mokit_step{step}_state{1+state}.out')
        if step>1 and os.path.exists(pyscf_inp):
            os.system(f'mv {pyscf_inp} {scratch_path}/{pyscf_inp}')
        if step>1 and os.path.exists(pyscf_out):
            os.system(f'mv {pyscf_out} {scratch_path}/{pyscf_out}')

    e1,e2=enelist
    g1,g2=gradlist
    ed=e2-e1

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

    # Backup QM geometry.
    qm=geom.qm_region()
    n=len(qm)
    qmxyz=geom.qm_xyz(xyz)
    with open(xyz_backup_qm,'a') as gen:
        gen.write(f'{n}\nMECP Step {step}\n')
        for i in range(n):
            gen.write(f'{atm[qm[i]-1]}'.ljust(2))
            for j in range(3):
                gen.write(f'{qmxyz[j+3*i]:.8f}'.rjust(16))
            gen.write('\n')

    # Add external electrostatic field.
    global esf
    if esf:
        global chg,fielddef
        from electrostatic_field import electrostatic_field
        esfene,esfgrad=electrostatic_field(xyz,chg,fielddef)

        # Delete ESF gradient for QM atoms.
        qmids=geom.qm_region()
        for i in range(len(esfgrad)):
            if 1+i in qmids:
                esfgrad[i]=0.

        g1=(numpy.array(g1)+numpy.array(esfgrad)).tolist()
        g2=(numpy.array(g2)+numpy.array(esfgrad)).tolist()

    # Delete gradient of atoms not in active region.
    g1=geom.active_gradient(g1)
    g2=geom.active_gradient(g2)

    return e1,g1,e2,g2

def optstep(XNew,step):
    print(f'Optimization steps: {step}')

    global scratch_path,xyz_backup,ene_backup
    if step==1:
        os.system(f'rm -rf {xyz_backup} {ene_backup} {scratch_path}; mkdir {scratch_path}')
        with open(ene_backup,'w') as gen:
            gen.write('Optimization\nStep  Energy\n')

    xyz=numpy.matrix.tolist(XNew)[0]
    # Recover coordinates of the whole QM/MM system.
    xyz=geom.upd_active_xyz(xyz)
    inp4opt(xyz,step)

    print('Run ChemShell eandg calculation.')

    chm=f'{system_name}_step{step}.chm'
    out=f'{system_name}_step{step}.out'

    os.system(f'echo "set step {step} # Note: for MOKIT interface" > xmecp.tmp\n')

    # If you want to optimize specific target state in CASSCF job, please set state1 {mult} {root} in #section mecp.
    global mecp_wordlist
    os.system(f'echo "set target_mult {mecp_wordlist[0]} # Note: for PySCF SA-CASSCF interface" >> xmecp.tmp\n')
    os.system(f'echo "set target_root {mecp_wordlist[2]} # Note: for PySCF SA-CASSCF interface" >> xmecp.tmp\n')

    with os.popen(f'chemsh < {chm} > {out}') as run:
        null=run.read()

    e=f'{system_name}.e'
    g=f'{system_name}.g'

    e1=energy(e)
    g1=gradient(g)

    orca_inp=f'{system_name}_orca.inp'
    orca_out=f'{system_name}_orca.out'
    orca_gbw=f'{system_name}_orca.gbw'

    # Delete useless files.
    os.system(f'rm -f {system_name}.e {system_name}.g FIELD CONFIG CONTROL hybrid* pointcharges.xyz')
    os.system(f'rm -f {system_name}_orca.pcgrad {system_name}_orca.engrad {system_name}_orca.cis {system_name}_orca.densities {system_name}_orca_property.txt')
    # Backup ChemShell files
    os.system(f'mv {chm} {scratch_path}/{system_name}_step{step}.chm')
    os.system(f'mv {out} {scratch_path}/{system_name}_step{step}.out')

    # Backup ORCA files. (if there are)
    if os.path.exists(orca_inp):
        os.system(f'mv {orca_inp} {scratch_path}/{system_name}_orca_step{step}.inp')
    if os.path.exists(orca_out):
        os.system(f'mv {orca_out} {scratch_path}/{system_name}_orca_step{step}.out')
    if os.path.exists(orca_gbw):
        os.system(f'mv {orca_gbw} {scratch_path}/{system_name}_orca_step{step}.gbw')

    mokit_gjf=f'{system_name}_mokit.gjf'
    mokit_out=f'{system_name}_mokit.out'
    pyscf_inp=f'{system_name}_step{step}_state1_CASSCF.py'
    pyscf_out=f'{system_name}_step{step}_state1_CASSCF.out'
    # Backup MOKIT and PySCF files. (if there are)
    if step==1 and os.path.exists(mokit_gjf):
        os.system(f'mv {mokit_gjf} {scratch_path}/{system_name}_mokit_step{step}_state1.py')
    if step==1 and os.path.exists(mokit_out):
        os.system(f'mv {mokit_out} {scratch_path}/{system_name}_mokit_step{step}_state1.out')
    if step>1 and os.path.exists(pyscf_inp):
        os.system(f'mv {pyscf_inp} {scratch_path}/{pyscf_inp}')
    if step>1 and os.path.exists(pyscf_out):
        os.system(f'mv {pyscf_out} {scratch_path}/{pyscf_out}')

    with open(ene_backup,'a') as backup:
        backup.write(f'{f"{step}".rjust(4)}  {f"{e1:.8f}".rjust(16)}\n')

    # Backup geometry.
    global atm
    n=len(atm)
    with open(xyz_backup,'a') as gen:
        gen.write(f'{n}\nOptimization Step {step}\n')
        for i in range(n):
            gen.write(f'{atm[i]}'.ljust(2))
            for j in range(3):
                gen.write(f'{xyz[j+3*i]:.8f}'.rjust(16))
            gen.write('\n')

    # Backup QM geometry.
    qm=geom.qm_region()
    n=len(qm)
    qmxyz=geom.qm_xyz(xyz)
    with open(xyz_backup_qm,'a') as gen:
        gen.write(f'{n}\nOptimization Step {step}\n')
        for i in range(n):
            gen.write(f'{atm[qm[i]-1]}'.ljust(2))
            for j in range(3):
                gen.write(f'{qmxyz[j+3*i]:.8f}'.rjust(16))
            gen.write('\n')

    # Add external electrostatic field.
    global esf
    if esf:
        global chg,fielddef
        from electrostatic_field import electrostatic_field
        esfene,esfgrad=electrostatic_field(xyz,chg,fielddef)

        # Delete ESF gradient for QM atoms.
        qmids=geom.qm_region()
        for i in range(len(esfgrad)):
            if 1+i in qmids:
                esfgrad[i]=0.

        g1=(numpy.array(g1)+numpy.array(esfgrad)).tolist()

    # Delete gradient of atoms not in active region.
    g1=geom.active_gradient(g1)

    return e1,g1

