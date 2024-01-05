import optimizer
from parameters import *

print('Entering XMECP program.')
print('Report error: xujiawei@fjirsm.ac.cn\n')

if layer=='qm':
    if qm_software=='g16':
        from rung16 import *
    elif qm_software=='orca':
        from runorca import *
    elif qm_software=='qchem':
        from runqchem import *
    elif qm_software=='pyscf':
        from runpyscf import *
        keywords(qm_template)
    else:
        print('Unsupported QM software! Exit program!')
    optimizer.opt(geom0)
elif layer=='qmmm':
    if qm_software=='orca':
        from runchemsh_orca import *
    if qm_software=='pyscf':
        from runchemsh_pyscf import *
    optimizer.opt(geom0)
else:
    print('Unsupported layer definition! Exit program!')
