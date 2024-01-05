def jobtype(inputfile):
    job=['mecp','qmmm']

    avail_job0=['mecp_harvey','mecp_harvay_noed','mecp_bpupd','mecp_pf','opt','mecp']
    avail_job1=['qm','qmmm']
    
    # Available: mecp_harvey, mecp_harvay_noed, mecp_bpupd / qm, qmmm
    #   mecp_harvey: Harvey's MECP locating method.
    #   mecp_harvay_noed: Harvey's MECP locating method, but ignoring state energy difference as convergence criterion.
    #       Note: Be very careful with your result in case using this method for MECP between states with same symmetry.
    #   mecp_bpupd: Branching plane updating method for gradient projection algorithm.

    readinput=open(inputfile,'r')
    for line in readinput.readlines():
        L=line.split()
        if len(L)!=0:
            if L[0]=='#jobtype':
                try:
                    job[0]=L[1]
                    job[1]=L[2]
                except:
                    pass
                break
    readinput.close()
    if (job[0] in avail_job0) and (job[1] in avail_job1):
        pass
    else:
        print('Unsupported job type! Exit program!')
        exit()
    if job[0]=='mecp':
        job[0]='mecp_bpupd'
    return job

def section_control(inputfile):
    control=False

    ### >>> Parameters expected in control section >>> ###
    job_name='mecp'                 # Name of the job. This is used as name of Gaussian input or scratch path.
    system_name='example_system'    # Name of the system. This should be held the same for .prmtop, .pdb and .c files.
    nproc=32                        # Number of parallel processors.
    memory=32                       # GB.
    control_wordlist=[job_name,system_name,nproc,memory]
    ### <<< Parameters expected in input file <<< ###

    readinput=open(inputfile,'r')
    for line in readinput.readlines():
        L=line.split()
        if len(L)!=0:
            if L[0]=='#section' and L[1]=='control':
                control=True
                continue
            elif L[0]=='#end':
                control=False
        if control and len(L)!=0:
            if L[0]=='job_name':
                job_name=L[1]
            elif L[0]=='system_name':
                system_name=L[1]
            elif L[0]=='nproc':
                nproc=int(L[1])
            elif L[0]=='memory':
                memory=int(L[1])
    readinput.close()
    control_wordlist=[job_name,system_name,nproc,memory]
    return control_wordlist

def section_qm(inputfile):
    qm=False

    ### >>> Parameters expected in qm section >>> ###
    qm_software='orca'              # Software for QM calculation.
    qm_path='/home/jwxu/orca5/orca' # Path of QM software executable file.
    qm_charge=0                     # Charge of QM region.
    qm_mult=1                       # Multiplicity of QM region.
    qm_wordlist=[qm_software,qm_path,
                 qm_charge,qm_mult]
    ### <<< Parameters expected in qm section <<< ###

    readinput=open(inputfile,'r')
    for line in readinput.readlines():
        L=line.split()
        if len(L)!=0:
            if L[0]=='#section' and L[1]=='qm':
                qm=True
                continue
            elif L[0]=='#end':
                qm=False
        if qm and len(L)!=0:
            if L[0]=='qm_software':
                qm_software=L[1]
            elif L[0]=='qm_path':
                qm_path=L[1]
            elif L[0]=='qm_charge':
                qm_charge=int(L[1])
            elif L[0]=='qm_mult':
                qm_mult=int(L[1])
    readinput.close()
    qm_wordlist=[qm_software,qm_path,
                 qm_charge,qm_mult]
    return qm_wordlist

def section_mecp(inputfile):
    mecp=False

    ### >>> Parameters expected in mecp section >>> ###
    state1_mult,state1_root=1,0     # Multiplicity and root of the first state to be involved in MECP search.
                                    # Note: <state1_root=0> stands for ground state and can only be used for <state1_mult=qm_mult>.
    state2_mult,state2_root=3,1     # Multiplicity and root of the second state to be involved in MECP search.
    if_td1,if_td2=False,False       # Whether to treat the lowest one state with different multiplicity to the ground state as an excited state.
                                    # Note: When if_td1 or if_td2 is set to False, unrestricted Kohn-Sham method is applied, which is by-default and recommended.
    uks=False                       # Whether to calculation unrestricted singlet ground state.
    sf=False                        # Spin-flip TD-DFT.
    dynamic=True                    # Use dynamic gradient projection algorithm.
    pf_alpha=0.02                   # Alpha for penalty function algorithm.
    pf_sigma=3.50                   # Sigma for penalty function algorithm.
    pf_thresh=[0.000001,0.00500]    # Thresholds for penalty function algorithm.
    mecp_wordlist=[state1_mult,state2_mult,
                   state1_root,state2_root,
                   if_td1,if_td2,uks,
                   sf,
                   dynamic,
                   pf_alpha,pf_sigma,pf_thresh]
    ### <<< Parameters expected in mecp section <<< ###

    qm_wordlist=section_qm(inputfile)
    qm_mult=qm_wordlist[3]

    global conver

    readinput=open(inputfile,'r')
    for line in readinput.readlines():
        L=line.split()
        if len(L)!=0:
            if L[0]=='#section' and L[1]=='mecp':
                mecp=True
                continue
            elif L[0]=='#end':
                mecp=False
        if mecp and len(L)!=0:
            if L[0]=='state1':
                state1_mult,state1_root=int(L[1]),int(L[2])
                if state1_mult!=qm_mult and state1_root==1:
                    if len(L)==4 and (L[3]=='true' or L[3]=='True'):
                        if_td1=True
            elif L[0]=='state2':
                state2_mult,state2_root=int(L[1]),int(L[2])
                if state2_mult!=qm_mult and state2_root==1:
                    if len(L)==4 and (L[3]=='true' or L[3]=='True'):
                        if_td2=True
            elif L[0]=='u' or L[0]=='uks' or L[0]=='unrestricted':
                uks=True
            elif L[0]=='sf' or L[0]=='spin-flip':
                if L[1]=='true':
                    sf=True
            elif L[0]=='dynamic':
                if L[1]=='false':
                    dynamic=False
            elif L[0]=='pf_alpha':
                pf_alpha=float(L[1])
            elif L[0]=='pf_sigma':
                pf_sigma=float(L[1])
            elif L[0]=='pf_thresh':
                for i in range(len(L)-1):
                    conver[i]=float(L[1+i])
    readinput.close()
    mecp_wordlist=[state1_mult,state2_mult,
                   state1_root,state2_root,
                   if_td1,if_td2,uks,
                   sf,
                   dynamic,
                   pf_alpha,pf_sigma,pf_thresh]
    return mecp_wordlist

def section_templates(inputfile):
    templates=False

    ### >>> Parameters expected in templates section >>> ###
    qm_template='qm.tpl'            # Template for QM input.
    chemsh_template='chemsh.tpl'    # Template for ChemShell input.
                                    # Note: Contents in template files will be directly added into input files for QM softwares and ChemShell.
    templates_wordlist=[qm_template,chemsh_template]
    ### <<< Parameters expected in templates section <<< ###

    readinput=open(inputfile,'r')
    for line in readinput.readlines():
        L=line.split()
        if len(L)!=0:
            if L[0]=='#section' and L[1]=='templates':
                templates=True
                continue
            elif L[0]=='#end':
                templates=False
        if templates and len(L)!=0:
            if L[0]=='qm_template':
                qm_template=L[1]
            elif L[0]=='chemsh_template':
                chemsh_template=L[1]
    readinput.close()
    templates_wordlist=[qm_template,chemsh_template]
    return templates_wordlist

def section_opt(inputfile):
    opt=False
    addB,addA,addD=[],[],[]

    ### >>> Parameters expected in optimization section >>> ###
    method='gdiis'                 # Optimization algorithm. Available: BFGS, GDIIS and GEDIIS.
    maxcycle=128                    # Max optimization cycles.
    maxstep=0.10                    # Max optimization stepsize.
    restart=False                   # Whether to restart from a structure in xyz_backup file.
    restart_from=0                  # Restart from which structure in xyz_backup file.
    conver=[0.000050,               # Convergence criterions for Energy difference (for MECP).
            0.002500,               # Convergence criterions for RMS displacement.
            0.004000,               # Convergence criterions for max displacement.
            0.000700,               # Convergence criterions for max gradient.
            0.000500]               # Convergence criterions for RMS gradient.
    insert=1                        # Number of points for linear insert method.
    internal=False                  # Use internal coordinates for linear insert method.
    add_ints=[addB,addA,addD]       # Add internal coordinates that must be included in linear insert.
    insert_list=[insert,internal,add_ints]
    opt_wordlist=[method,maxcycle,maxstep,restart,restart_from,conver,insert_list]
    ### <<< Parameters expected in optimization section <<< ###

    if jobtype(inputfile)[0] in ['opt','ts']:
        conver=[0.000050,
                0.001200,
                0.001800,
                0.000450,
                0.000300] # When optimization, set to Gaussian default (opt=Tight).

    readinput=open(inputfile,'r')
    for line in readinput.readlines():
        L=line.split()
        if len(L)!=0:
            if L[0]=='#section' and (L[1]=='opt' or L[1]=='optimization'):
                opt=True
                continue
            elif L[0]=='#end':
                opt=False
        if opt and len(L)!=0:
            if L[0]=='method':
                method=L[1]
            elif L[0]=='maxcycle':
                maxcycle=int(L[1])
            elif L[0]=='maxstep':
                maxstep=float(L[1])
            elif L[0]=='restart':
                restart=True
                restart_from=int(L[1])
            elif L[0]=='conver':
                for i in range(len(L)-1):
                    conver[i]=float(L[1+i])
            elif L[0]=='linear':
                if int(L[1])>=1:
                    insert=int(L[1])
                else:
                    print('At least one insert point. Set to default: 1.')
                try:
                    if L[1]=='int' or L[1]=='internal':
                        internal=True
                except:
                    pass
            elif L[0]=='add_int' or L[0]=='add_internal':
                if len(L)==4 and L[1]=='B':
                    addB.append([int(2),int(3)])
                if len(L)==5 and L[1]=='A':
                    addA.append([int(2),int(3),int(4)])
                if len(L)==6 and L[1]=='D':
                    addD.append([int(2),int(3),int(4),int(5)])
                add_ints=[addB,addA,addD]
                insert_list=[insert,internal,add_ints]
    readinput.close()
    opt_wordlist=[method,maxcycle,maxstep,restart,restart_from,conver,insert_list]
    return opt_wordlist

def section_field(inputfile):
    field=False

    ### >>> Parameters expected in external electrostatic field section >>> ###
    esf=False                       # Turn on external electrostatic field.
    plane=None                      # Atoms used to define external electrostatic field.
    intensity=0.0                   # Intensity of external electrostatic field.
    charge='prmtop'                 # Charge file. File with suffix of '.prmtop' of '.chg'.
    direction=None                  # When len(plane)>=3, used to define direction. (optional when len(plane)>3)
    field_wordlist=[plane,intensity,charge,direction,esf]
    ### <<< Parameters expected in external electrostatic field section <<< ###

    def parse_id(ids):
        ids=ids.replace(',',' ')
        from geom import is_number,all_numbers
        idlist=[]
        for i in ids.split():
            if is_number(i):
                idlist.append(int(i))
            elif '-' in i:
                l=i.split('-')
                if all_numbers(l) and len(l)==2:
                    for j in range(int(l[0]),1+int(l[1])):
                        idlist.append(j)
        idlist.sort()
        return idlist

    readinput=open(inputfile,'r')
    for line in readinput.readlines():
        L=line.split()
        if len(L)!=0:
            if L[0]=='#section' and L[1]=='field':
                field=True
                continue
            elif L[0]=='#end':
                field=False
        if field and len(L)!=0:
            if L[0]=='plane':
                plane=parse_id(line)
            elif L[0]=='intensity':
                intensity=float(L[1])
            elif L[0]=='charge':
                charge=L[1]
            elif L[0]=='direction':
                direction=int(L[1])
            elif L[0]=='esf' and L[1]=='true':
                esf=True
    readinput.close()
    field_wordlist=[plane,intensity,charge,direction,esf]
    return field_wordlist
