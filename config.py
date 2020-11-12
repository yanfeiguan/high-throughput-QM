import os

# SQL paramters
SQL_DBNAME = 'nmrtest'
SQL_USER = 'nmruser'
SQL_PASSWORD = '123456'
SQL_HOST = 'rmg'
SQL_PORT = '5432'

# conformer searching parameters
MIN_CONFORMERS = 10
MAX_CONFORMERS = 20
CONFORMER_E_THRESHOLD = 2.5                 # kcal/mol
CONFORMER_RMSD_THRESHOLD = 0.4

# running parameters
SCRATCHDIR = os.environ['WORKDIR']          # by default the WORKDIR environment varialbe is set in the submitting bash script, but you can provide your own scarch folder
PROJECTDIR = '/home/yanfeig/nmr/test'
WALLTIME = 48                               # time limit for job in hours
NPROCS = 16                                 # number of processes
MEM = '35000mb'                             # memory with unit (used in qm_worker)

# optimization parameters
XTB_CMD = 'xtb'                             # if the xtb is installed in the conda environment, just provide the command here (xtb)
XTB_CONFORMER_E_THRESHOLD = 2.5             # kcal/mol

# Gaussian parameters
QM_CMD = 'g16'
LEVEL_OF_THEORY = "mPW1PW91/6-31G"
SOLVENT = 'dimethylsulfoxide'
ADDITIONAL_PARAMS = 'scf=(xqc,maxconventionalcycles=400) nosymm'
