# Imports
###############################################################################
import os
import argparse
import progressbar
import Bio.PDB as bpp
import textwrap as tw

# License
###############################################################################
'''

ProtCHOIR: A tool for generation of homo oligomers from pdb structures

Authors: Torres, P.H.M.; Blundell, T.L.

[The University of Cambridge]

Contact info:
Department Of Biochemistry
University of Cambridge
80 Tennis Court Road
Cambridge CB2 1GA
E-mail address: monteirotorres@gmail.com

This project is licensed under Creative Commons license (CC-BY-4.0)

'''

# Description
###############################################################################
description = tw.dedent("""
    \033[1;93m+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    +-+-+-+-+-+-+-+-+-+-+-+-+-+  \033[1;95mProtCHOIR  \033[1;93m+-+-+-+-+-+-+-+-+-+-+-+-+-+
    +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    \033[1;95m
                  Because proteins only sing...
    \033[1;93m
                                 ...when they are together!
    \033[1;0m
      Copyright (C) 2018  Torres, P.H.M.; Blundell, T.L.

    \033[1;93m
                      [The University of Cambridge]

    \033[1;0m
          This program comes with ABSOLUTELY NO WARRANTY

      This program was devised to create homo-oligomeric structures
              based on selected subsets of the PDB databank.

    The software relies on well documented computational biology tools
        such as Gesamt(1), PISA(2), Parasail(3), Biopython(4) and
           Modeller(5) to generate and assess the final models.


     Please cite:
     ProtCHOIR: Protein Complexes and Homo-Oligomeric Interfaces
     Resolver. Torres,P.H.M.; Blundell, T.L.\033[1;0m

     +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
          """)

epilogue = tw.dedent("""
    Protchoir generates a summary file, whose columns are ordered as follows:

    1 Input           7 H30Score      13 BestModel      19 ProtCHOIR
    2 Seq.Mode        8 Template      14 Molprobity     20 Runtime
    3 Templated          9 Chains        15 RMSD           21 Exit
    4 Length         10 Identity      16 Quality
    5 TMSpans        11 Coverage      17 Surface
    6 LikelyState    12 Av.QScore     18 Interfaces
    """)


# functions
###############################################################################


def create_choir_conf():
    conf_file = "CHOIR.cfg"
    with open(conf_file, 'w') as cf:
        cf.write(tw.dedent("""
        # PyMol Executable
        pymol_exe = /usr/local/bin/pymol

        # CCP4 installation directory
        ccp4_base = /opt/ccp4/ccp4-7.0

        # PSI-Blast Executable
        psiblast_exe = /usr/bin/psiblast

        # blastdbcmd Executable
        blastdbcmd_exe = /usr/bin/blastdbcmd

        # makeblastdb_exe
        makeblastdb_exe = /usr/bin/makeblastdb

        # MAFFT Executable
        mafft_exe = /usr/bin/mafft

        # TMHMM-2.0 executable
        tmhmm2_exe = /opt/tmhmm-2.0c/bin/tmhmm

        # Root directory for the ProtCHOIR Database
        choirdb = /data/choirdb"""))

    print(clrs['g']+'Configuration file created!'+clrs['n']+' Please'+clrs['y']+' EDIT ITS CONTENTS '+clrs['n']+'to match your environment and run ProtCHOIR again.')


# Define Global Variables
###############################################################################
global p
global io
global args
global maxasa
global aa_bgf
global choirdb
global widgets
global workdir
global uniref50
global homoblast
global oligo_dict
global pdb_archive
global pdb1_archive
global oligo_dict_inv
global ges_homo_archive
global pdb_homo_archive

# Parse command line arguments
###############################################################################
def argument_parsing():
    parser = argparse.ArgumentParser(prog='ProtCHOIR',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=description,
                                     epilog=epilogue)

    parser.add_argument('--version', action='version',
                    version='%(prog)s 1.2.14')

    parser.add_argument('-f', '--file',
                        dest='input_file',
                        type=str,
                        metavar='',
                        help='File containing the paths of pdbs to be oligomerized')

    parser.add_argument('-c, --max-candidates',
                        dest='max_candidates',
                        type=int, default=5,
                        metavar='',
                        help='Maximum number of candidate PDB structures to consider')

    parser.add_argument('-m', '--models',
                        dest='models',
                        type=int, default=1,
                        metavar='',
                        help='Number of Modeller-generated models')

    parser.add_argument('-q', '--qscore-cutoff',
                        dest='qscore_cutoff',
                        type=float, default=0.4,
                        metavar='',
                        help='Determines the cut-off for average GESAMT Q-score When comparing the protomer with each chain of template oligomers.  Ignored when using the ignore-protomer option.')

    parser.add_argument('-t', '--tolerance',
                        dest='tolerance',
                        type=float, default=15,
                        metavar='',
                        help='Determines the maximum allowed percent difference between PSI-BLAST scores in order to attempt modelling')

    parser.add_argument('-s', '--similarity-cutoff',
                        dest='similarity_cutoff',
                        type=float, default=25,
                        metavar='',
                        help='Determines the cut-off for similarity over protein stretches to still attempt modelling')


    parser.add_argument('--bad-stretches',
                        dest='bad_streches',
                        type=int, default=2,
                        metavar='',
                        help='Determines the number of low similarity stretches per chain that still allows modelling')

    parser.add_argument('-r', '--refine-level',
                        dest='refine_level',
                        type=int, default=0,
                        metavar='',
                        help='[0] Refine using very-fast VTFM.[1] Refine using fast VTFM.[2] Refine using normal VTFM.[3] Refine using slow VTFM.[4] Refine using slow VTFM and slow MD protocols.')


    parser.add_argument('-a', '--assessment',
                        dest='assessment',
                        metavar='',
                        type=str, default='I',
                        help='Define oligomer assessment methods. Input a string containing any combination of the following capital letters: [G]Gesamt, [M]Molprobity, [I]Interfaces (not compatible with sequence mode) or [N]None')

    parser.add_argument('--sequence-mode',
                        dest='sequence_mode',
                        action='store_true',
                        default=False,
                        help='Instructs ProtCHOIR NOT to consider the input protomer (sequence mode)')

    parser.add_argument('--ignore-templated',
                        dest='ignoretemplated',
                        action='store_true',
                        default=False,
                        help='Instructs ProtCHOIR NOT to consider Templated-selected templates')

    parser.add_argument('--plot-topologies',
                        dest='plot_topologies',
                        action='store_true',
                        default=False,
                        help='Plot oligmerization topologies of all candidate oligomeric templates')

    parser.add_argument('--psiblast-threads',
                        dest='psiblast_threads',
                        type=int,
                        metavar='',
                        help='Number of threads to use for PSI-BLAST, defaults to the number of processors')

    parser.add_argument('--modeller-threads',
                        dest='modeller_threads',
                        type=int,
                        metavar='',
                        help='Number of threads to use for MODELLER, defaults to the smallest vaule between the number of processors or number of models')

    parser.add_argument('--multiprocess',
                        dest='multiprocess',
                        action='store_true',
                        default=False,
                        help='Defines whether python multiprocessing should be enabled for compatible lenghty tasks')

    parser.add_argument('--single-core',
                        dest='force_single_core',
                        action='store_true',
                        default=False,
                        help='Forces all individual tasks and programs to run in a single core. Disables multiprocessing.')

    parser.add_argument('--skip-conservation',
                        dest='skip_conservation',
                        action='store_true',
                        default=False,
                        help='Skip PSI-Blast against UniRef50 and entropy calculations')

    parser.add_argument('--allow-monomers',
                        dest='allow_monomers',
                        action='store_true',
                        default=False,
                        help='Instructs ProtCHOIR to build monomeric structures if template oligomers are not top hits. Only available in sequence mode.')

    parser.add_argument('--symmetry',
                        dest='symmetry',
                        action='store_true',
                        default=False,
                        help='Run Modeller with symmetry constraints for all chains.WARNING: Enforcing symmetry might take a very long time!')

    parser.add_argument('--repeat-optimization',
                        dest='repeat_opt',
                        type=int, default=0,
                        metavar='',
                        help='Defines how many times the Modeller optimization protocol should be executed')

    parser.add_argument('--generate-report',
                        dest='generate_report',
                        action='store_true',
                        default=False,
                        help='Creates a final HTML report for each generated model (forces -a MIG and --plot-topologies)')

    parser.add_argument('-z', '--zip-output',
                        dest='zip_output',
                        type=int, default=0,
                        metavar='',
                        help='Defines the compression level. [0] No compression, [1] partial compression, [2] full compression')

    parser.add_argument('-u', '--update-databases',
                        dest='update',
                        action='store_true',
                        default=False,
                        help='Updates databases')

    parser.add_argument('-v', '--verbose',
                        dest='verbosity',
                        action='count',
                        default=0,
                        help='Controls verbosity')

    parser.add_argument('--conf',
                        dest='config_file',
                        type=str,
                        metavar='',
                        help='Configuration file containing external executable paths')

    initial_args = parser.parse_args()

    return initial_args

initial_args = argument_parsing()


# Define Global Variables
###############################################################################
# Dictionary relating number of chains and oligomeric name up to 60
oligo_dict = {1: 'MONOMERIC', 2: 'DIMERIC', 3: 'TRIMERIC', 4: 'TETRAMERIC',
              5: 'PENTAMERIC', 6: 'HEXAMERIC', 7: 'HEPTAMERIC', 8: 'OCTAMERIC',
              9: 'NONAMERIC', 10: 'DECAMERIC', 11: 'UNDECAMERIC',
              12: 'DODECAMERIC', 13: 'TRIDECAMERIC', 14: 'TETRADECAMERIC',
              15: 'PENTADECAMERIC', 16: 'HEXADECAMERIC', 17: 'HEPTADECAMERIC',
              18: 'OCTADECAMERIC', 19: 'NONADECAMERIC', 20: 'EICOSAMERIC',
              21: '21-MERIC', 22: '22-MERIC', 23: '23-MERIC', 24: '24-MERIC',
              25: '25-MERIC', 26: '26-MERIC', 27: '27-MERIC', 28: '28-MERIC',
              29: '29-MERIC', 30: '30-MERIC', 31: '31-MERIC', 32: '32-MERIC',
              33: '33-MERIC', 34: '34-MERIC', 35: '35-MERIC', 36: '36-MERIC',
              37: '37-MERIC', 38: '38-MERIC', 39: '39-MERIC', 40: '40-MERIC',
              41: '41-MERIC', 42: '42-MERIC', 43: '43-MERIC', 44: '44-MERIC',
              45: '45-MERIC', 46: '46-MERIC', 47: '47-MERIC', 48: '48-MERIC',
              49: '49-MERIC', 50: '50-MERIC', 51: '51-MERIC', 52: '52-MERIC',
              53: '53-MERIC', 54: '54-MERIC', 55: '55-MERIC', 56: '56-MERIC',
              57: '57-MERIC', 58: '58-MERIC', 59: '59-MERIC', 60: '60-MERIC',}

# Inverted dictionary relating oligomeric name and number of chains
oligo_dict_inv = {v: k for k, v in oligo_dict.items()}

# Dictionary relating letters to numbers
alphanum = {'A':'1','B':'2','C':'3','D':'4','E':'5','F':'6','G':'7','H':'8',
            'I':'9','J':'10','K':'11','L':'12','M':'13','N':'14','O':'15',
            'P':'16','Q':'17','R':'18','S':'19','T':'20','U':'21','V':'22',
            'W':'23','X':'24','Y':'25','Z':'26','a':'27','b':'28','c':'29',
            'd':'30','e':'31','f':'32','g':'33','h':'34','i':'35','j':'36',
            'k':'37','l':'38','m':'39','n':'40','o':'41','p':'42','q':'43',
            'r':'44','s':'45','t':'46','u':'47','v':'48','w':'49','x':'50',
            'y':'51','z':'52','1':'53','2':'54','3':'55','4':'56','5':'57',
            '6':'58','7':'59','8':'60'}

# Inverted dictionary relating numbers to letters
numalpha = {v: k for k, v in alphanum.items()}

# Pdb format
pdb_format = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n"

# Van der Waals radius of most common atoms
VdWRadius = {'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'S': 1.8, 'CL': 1.75 }

# Pdb parser class from Biopython
p = bpp.PDBParser(PERMISSIVE=0, QUIET=True)

# Pdb writer class from Biopython
io = bpp.PDBIO()

# Dictionary for the output colors
clrs = {'r': "\033[1;91m",
        'g': "\033[1;92m",
        'y': "\033[1;93m",
        'b': "\033[1;94m",
        'p': "\033[1;95m",
        'c': "\033[1;96m",
        'n': "\033[1;0m"}

# Format for the progressbar shown on stdout
widgets = [' [', progressbar.SimpleProgress(), '] ',
           progressbar.Bar(),
           progressbar.Percentage(),
           ' (', progressbar.AdaptiveETA(), ') ']

aa3to1 = {'CYS': 'C', 'ASP': 'D', 'GLN': 'Q', 'ILE': 'I',
          'ALA': 'A', 'TYR': 'Y', 'TRP': 'W', 'HIS': 'H',
          'LEU': 'L', 'ARG': 'R', 'VAL': 'V', 'GLU': 'E',
          'PHE': 'F', 'GLY': 'G', 'MET': 'M', 'ASN': 'N',
          'PRO': 'P', 'SER': 'S', 'LYS': 'K', 'THR': 'T',
          'MSE': 'M', 'CSE': 'U', 'GLH': 'E', 'HID': 'H',
          'HIE': 'H', 'HIP': 'H', 'HYP': 'P', 'ASX': 'B',
          'GLX': 'Z', 'MME': 'M', 'LYZ': 'K'}


maxasa = {'ALA':129.0,'ARG':274.0,'ASN':195.0,'ASP':193.0,
          'CYS':167.0,'GLU':223.0,'GLN':225.0,'GLY':104.0,
          'HIS':224.0,'ILE':197.0,'LEU':201.0,'LYS':236.0,
          'MET':224.0,'PHE':240.0,'PRO':159.0,'SER':155.0,
          'THR':172.0,'TRP':285.0,'TYR':263.0,'VAL':174.0,
          'MME':224.0,'MSE':224.0, 'UNK':129.0}
'''
Theoretical maximum solvent accessible area values retrieved from:
TIEN, M. Z., MEYER, A. G., SYDYKOVA, D. K., SPIELMAN, S. J., & WILKE, C. O.
Maximum allowed solvent accessibilites of residues in proteins. PloS one, 2013.
* B(ASX) takes the same value as D(ASP)
* Z(GLX) takes the same value as E(GLU)
'''

#------------------------------------------------------------------------------
hydro = {'ALA': -0.17, 'ARG': -0.81, 'ASN': -0.42, 'ASP': -1.23,
         'CYS': 0.24, 'GLN': -2.02, 'GLU': -0.58, 'GLY': -0.01,
         'HIS': -0.96, 'ILE': 0.31, 'LEU': 0.56, 'LYS': -0.99,
         'MET': 0.23, 'PHE': 1.13, 'PRO': -0.45, 'SER': -0.13,
         'THR': -0.14, 'TRP': 1.85, 'TYR': 0.94, 'VAL': -0.07,
         'MME': 0.23, 'MSE': 0.23, 'UNK': -0.17}
'''
Hydropobicity values retrieved from:
WIMLEY, W. C.; WHITE, S. H. Experimentally determined hydrophobicity scale for
proteins at membrane interfaces. Nature Structural Biology, 1996.
* B(ASX) takes the same value as D(ASP)
* Z(GLX) takes the same value as E(GLU)
'''

#------------------------------------------------------------------------------
aa_bgf = {'A': 0.078, 'C': 0.024, 'D': 0.052, 'E': 0.059, 'F': 0.044,
          'G': 0.083, 'H': 0.025, 'I': 0.062, 'K': 0.056, 'L': 0.092,
          'M': 0.024, 'N': 0.041, 'P': 0.043, 'Q': 0.034, 'R': 0.051,
          'S': 0.059, 'T': 0.055, 'V': 0.072, 'W': 0.014, 'Y': 0.034,
          'B': 0.052, 'Z': 0.059}

'''
Amino acids background frequencies retrieved from:
CAPRA, J. A.; SINGH, M. Predicting functionally important residues from sequence
conservation. Bioinformatics, 2007.
* B(ASX) takes the same value as D(ASP)
* Z(GLX) takes the same value as E(GLU)
'''

# Initialise
###############################################################################
print(description)

# Retrieve the paths from provided configuration file
if (not initial_args.config_file or not os.path.isfile(initial_args.config_file)) and not os.path.isfile('CHOIR.cfg'):
    print('ProtCHOIR configuration file not found in the provided path')
    create_config = input('Do you wish to create it? (y/n)')
    if create_config == 'y' or create_config == 'Y' or create_config == 'YES' or create_config == 'yes' or create_config == 'Yes':
        create_choir_conf()
        quit()
    else:
        print('\n\nNo positive confirmation, please provide a valid configuration file.\n')
        quit()

elif not initial_args.config_file and os.path.isfile('CHOIR.cfg'):
    config_file = 'CHOIR.cfg'

elif initial_args.config_file:
    assert os.path.isfile(initial_args.config_file), clrs['r']+'\n\n Not able to find configuration file.\n\n Does "'+initial_args.config_file+'" exist?'+clrs['n']
    config_file = initial_args.config_file

for line in open(config_file, 'r'):
    if line.startswith('pymol_exe'):
        pymol_exe = line.split('=')[1].strip()
    if line.startswith('ccp4_base'):
        ccp4_base = line.split('=')[1].strip()
        pisa_exe = os.path.join(ccp4_base, 'bin', 'pisa')
        gesamt_exe = os.path.join(ccp4_base, 'bin', 'gesamt')
        molprobity_exe = os.path.join(ccp4_base, 'bin', 'molprobity.molprobity')
    if line.startswith('psiblast_exe'):
        psiblast_exe = line.split('=')[1].strip()
    if line.startswith('blastdbcmd_exe'):
        blastdbcmd_exe = line.split('=')[1].strip()
    if line.startswith('makeblastdb_exe'):
        makeblastdb_exe = line.split('=')[1].strip()
    if line.startswith('mafft_exe'):
        mafft_exe = line.split('=')[1].strip()
    if line.startswith('tmhmm2_exe'):
        tmhmm2_exe = line.split('=')[1].strip()
    if line.startswith('choirdb'):
        choirdb = line.split('=')[1].strip()

# Root directory for CHOIR module
choir_path = os.path.dirname(os.path.abspath( __file__ ))

# Directory containing the pdb mirror in "divided" scheme
pdb_archive = os.path.join(choirdb, 'pdb')

# Directory containing the pdb biounit mirror in "divided" scheme
pdb1_archive = os.path.join(choirdb, 'pdb1')

# Directory containing the ProtCHOIR database
ges_homo_archive = os.path.join(choirdb, 'gesamt_homo_archive')

# Directory containing the GESAMT archive
pdb_homo_archive = os.path.join(choirdb, 'pdb_homo_archive')

# Blast database on homo-oligomers
homoblast = os.path.join(pdb_homo_archive, 'sequences/homodb')

# Blast database on monomers
monoblast = os.path.join(pdb_homo_archive, 'sequences/monodb')

# Blast database on hetero-oligomers
heteroblast = os.path.join(pdb_homo_archive, 'sequences/heterodb')

# Blast database on UniRef50
uniref50 = os.path.join(choirdb, 'uniref50/uniref50')

# Blast database on seqres
seqres = os.path.join(choirdb, 'seqres/seqres')

# Blast database on pdbseq
pdbseq = os.path.join(choirdb, 'pdbseq/pdbseq')
