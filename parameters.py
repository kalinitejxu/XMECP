import os,sys,numpy,math
inp=sys.argv[1]

with open('input.tmp','w') as gen:
    with open(inp,'r') as rawinp:
        for line in rawinp.readlines():
            gen.write(line.lower())
os.system(f'mv input.tmp {inp}')

sasf=False

### >>> Constants >>> ###
bohr=0.52917721067
### <<< Constants <<< ###

import inp_proc
### >>> Input Parameters >>> ###
job=inp_proc.jobtype(inp)
algorithm       =                job[0]
layer           =                job[1]

control_wordlist=inp_proc.section_control(inp)
job_name        =   control_wordlist[0]
system_name     =   control_wordlist[1]
nproc           =   control_wordlist[2]
memory          =   control_wordlist[3]

qm_wordlist=inp_proc.section_qm(inp)
qm_software     =        qm_wordlist[0]
qm_path         =        qm_wordlist[1]
qm_charge       =        qm_wordlist[2]
qm_mult         =        qm_wordlist[3]

mecp_wordlist=inp_proc.section_mecp(inp)
state1_mult     =      mecp_wordlist[0]
state2_mult     =      mecp_wordlist[1]
state1_root     =      mecp_wordlist[2]
state2_root     =      mecp_wordlist[3]
if_td1          =      mecp_wordlist[4]
if_td2          =      mecp_wordlist[5]
uks             =      mecp_wordlist[6]
sf              =      mecp_wordlist[7]
dynamic         =      mecp_wordlist[8]
pf_alpha        =      mecp_wordlist[9]
pf_sigma        =      mecp_wordlist[10]
pf_thresh       =      mecp_wordlist[11]
sf_thresh       =      mecp_wordlist[12]

pf_thresh_step  =          pf_thresh[0]
pf_thresh_grad  =          pf_thresh[1]

templates_wordlist=inp_proc.section_templates(inp)
qm_template     = templates_wordlist[0]
chemsh_template = templates_wordlist[1]

opt_wordlist=inp_proc.section_opt(inp)
opt_method      =       opt_wordlist[0]
maxcycle        =       opt_wordlist[1]
maxstep         =       opt_wordlist[2]
restart         =       opt_wordlist[3]
restart_from    =       opt_wordlist[4]
conver          =       opt_wordlist[5]
insert          =       opt_wordlist[6][0]
internal        =       opt_wordlist[6][1]
add_ints        =       opt_wordlist[6][2]

thresh_de       =             conver[0]
thresh_rms      =             conver[1]
thresh_max_dis  =             conver[2]
thresh_max_g    =             conver[3]
thresh_rms_g    =             conver[4]

field_wordlist=inp_proc.section_field(inp)
plane           =     field_wordlist[0]
intensity       =     field_wordlist[1]
charge          =     field_wordlist[2]
direction       =     field_wordlist[3]
esf             =     field_wordlist[4]
### <<< Input Parameters <<< ###

### >>> Running Parameters >>> ###
scratch_path='scratch-'+job_name
xyz_backup,ene_backup='xyz_backup.xyz','ene_backup.ene'
xyz_backup_qm='xyz_backup_qm.xyz'
# Linear insert.
dEs=[]

# Initial geometry.
import geom
atm,geom0=[],[]
# Read restart geometry if required.
if restart:
    print('Restarted from old geometry backup file: "xyz_backup.xyz"')
    print('Current initial geometry is the No.'+str(restart_from)+' geometry in backup file.')
    atm,xyzs=geom.readxyzs('xyz_backup.xyz')
    geom0=xyzs[restart_from-1]
else:
    mol=sys.argv[2]
    atm,xyzs=geom.readxyzs(mol)
    geom0=xyzs[-1]

natms=len(atm)

e1s,g1s,e2s,g2s=[],[],[],[]

# Updating branching plane.
xs=[]
ys=[]
gdiffs=[]
gmeans=[]

# Penalty function.
a=pf_alpha
s=pf_sigma
pfF=[]

# Convergence.
de=0.0              # Energy Gap
rms_g=0.0           # RMS Gradient
max_g=0.0           # Maximium Gradient
rms=0.0             # RMS Displacement
max_dis=0.0         # Maximium Displacement

# SF-TD-DFT.
sf_ref_mult=qm_mult
if sf:
    if state1_root==0 or state2_root==0:
        print('The target root 0 stands for spin-flip reference state.')
        print('This can sometimes be meaningless. Please pay attention!')
        for i in range(2):
            if mecp_wordlist[2+i]==0:
                sf_ref_mult=mecp_wordlist[i]
                if state1_mult+state2_mult!=2*sf_ref_mult-2:
                    print('The state other than the reference state must corresponds to a flip-down state! Exit program!')
                    exit()
                break
    elif state1_mult==state2_mult:
        sf_ref_mult=state1_mult+2
    elif abs(state1_mult-state2_mult)==2:
        sf_ref_mult=max(state1_mult,state2_mult)
    else:
        print('Required states are not available for spin-flip framework at the same time! Exit program!')
        exit()

# External electrostatic field.
chg,fielddef=None,None
if esf:
    from electrostatic_field import readchg,readchg_from_prmtop
    # Get charges.
    chgfile=f'{system_name}.{charge}'
    if charge=='prmtop':
        chg=readchg_from_prmtop(chgfile)
    elif charge=='chg':
        chg=readchg(chgfile)
    fielddef=[intensity,plane]
    import geom
    try:
        if geom.is_number(direction):
            fielddef.append(direction)
    except:
        pass
### >>> Running Parameters >>> ###
