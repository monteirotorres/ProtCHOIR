# Imports
###############################################################################
import os
import sys
import gzip
import itertools
import collections
import Bio.PDB as bpp
import ProtCHOIR.Toolbox as pctools
from ProtCHOIR.Initialise import *

# pdb_file = '/data/choirdb/pdb1/c8/2c8i.pdb1.gz'
# pdb = gzip.open(pdb, 'rt')
# structure = p.get_structure('2c8i', pdb)
# structure.child_dict
# len(bpp.Selection.unfold_entities(structure, 'C'))
# structure2, test = pctools.split_states(structure)
#
# pctools.extract_seqs(structure2, 0)
# LICENSE
###############################################################################
'''

ProtCHOIR: A tool for generation of homo oligomers from pdb structures

Authors: Torres, P.H.M.; Malhotra, S.; Blundell, T.L.

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

# Dictionaries
###############################################################################

# Global Variables
###############################################################################

# Classes
###############################################################################

# Functions
###############################################################################
def make_local_template(best_oligo_template):
    middle_letters_best = best_oligo_template[1:3]
    best_template_file = os.path.join(pdb_homo_archive, middle_letters_best, best_oligo_template+".pdb.gz")
    clean_template_file = os.path.join(workdir, best_oligo_template+"_CHOIR_CleanTemplate.pdb")
    with open(clean_template_file, 'w') as f:
        for line in gzip.open(best_template_file, 'rt').readlines():
            if line.startswith('ATOM'):
                f.write(line)
    return clean_template_file


def extract_relevant_chains(pdb_file, relevant_chains):
    template_name = os.path.basename(pdb_file).split('_CHOIR_')[0]
    pname, structure, nchains = pctools.parse_any_structure(pdb_file)
    relevant_chains_file = os.path.join(workdir, template_name+"_CHOIR_RelevantChains.pdb")
    chains = bpp.Selection.unfold_entities(structure, 'C')
    io.set_structure(structure)
    io.save(relevant_chains_file, pctools.SelectChains(relevant_chains))

    return relevant_chains_file


def rename_relevant_chains(pdb_file):
    template_name = os.path.basename(pdb_file).split('_CHOIR_')[0]
    pname, structure, nchains = pctools.parse_any_structure(pdb_file)
    renamed_chains_file = os.path.join(workdir, template_name+"_CHOIR_RenamedChainsTemplate.pdb")
    chains = bpp.Selection.unfold_entities(structure, 'C')
    chains_dict = {}
    n = 1
    for chain in chains:
        original = chain.id
        new = numalpha[str(n)]
        chain.id = 'X'+new
        n += 1
        chains_dict[original] = new
    for chain in chains:
        chain.id = chain.id[1]
    io.set_structure(structure)
    io.save(renamed_chains_file)

    return renamed_chains_file, chains_dict


def restore_chain_identifiers(pdb_file, chains_dict, full_residue_mapping):
    pname, structure, nchains = pctools.parse_any_structure(pdb_file)
    restored_chains_file = os.path.join(workdir, pname+"_CHOIR_CorrectedChains.pdb")
    chains = bpp.Selection.unfold_entities(structure, 'C')
    str_id = structure.id
    new_structure = bpp.Structure.Structure(str_id)
    new_model = bpp.Model.Model(0)
    for original, current in chains_dict.items():
        for chain in chains:
            if chain.id == current:
                new_chain = bpp.Chain.Chain(current)
                new_chain.id = original
                for residue in chain:
                    new_residue = bpp.Residue.Residue(residue.id, residue.get_resname(), residue.get_segid())
                    if type(full_residue_mapping[current]) is collections.OrderedDict:
                        for atom in residue:
                            new_residue.add(atom)
                            new_residue.id = (' ', full_residue_mapping[current][residue.id[1]], ' ')
                    if type(full_residue_mapping[current]) is int:
                        for atom in residue:
                            new_residue.add(atom)
                            new_residue.id = (' ', full_residue_mapping[current]+residue.id[1], ' ')
                    new_chain.add(new_residue)
                new_model.add(new_chain)
    new_structure.add(new_model)
    io.set_structure(new_structure)
    io.save(restored_chains_file)
    return restored_chains_file


def alignment_from_sequence(best_oligo_template, renamed_chains_file, input_fasta, current_chain):
    input_basename = os.path.basename(input_fasta).split('_CHOIR_MonomerSequence.fasta')[0]
    print('Running '+clrs['b']+'MODELLER'+clrs['n']+' to align '+clrs['y']+input_basename+clrs['n']+' to '+clrs['y']+best_oligo_template+clrs['n']+' - Chain '+clrs['y']+current_chain+clrs['n'])
    output_basename = input_basename+'_'+best_oligo_template+current_chain
    genali_file = os.path.join(workdir, output_basename+'_CHOIR_Genali.py')
    temp_out = os.path.join(workdir, output_basename+'_CHOIR_Align2Dtemp.fasta')
    fasta_out = os.path.join(workdir, output_basename+'_CHOIR_Align2D.fasta')
    with open(genali_file, 'w') as f:
        f.write('from modeller import *\nenv = environ()\naln = alignment(env)\n')
        f.write("mdl = model(env, file='"+renamed_chains_file+"', model_segment=('FIRST:"+current_chain+"','LAST:"+current_chain+"'))\n")
        f.write("aln.append_model(mdl, align_codes='"+os.path.basename(renamed_chains_file)+"("+current_chain+")', atom_files='"+renamed_chains_file+"')\n")
        f.write("aln.append(file='"+input_fasta+"', align_codes='"+input_basename+"',alignment_format='FASTA', remove_gaps=False)\n")
        f.write('aln.align2d()\n')
        f.write("aln.write(file='"+temp_out+"', alignment_format='FASTA')")
    alignment_log = genali_file.replace('.py', '.log')
    script_name = os.path.basename(genali_file).split('.py')[0]
    temp = sys.stdout
    sys.stdout = open(alignment_log, 'w')
    __import__(script_name)
    sys.stdout.close()
    sys.stdout = temp
    with open(temp_out, 'r') as infile, open(fasta_out, 'w') as outfile:
        n = 1
        for line in infile.readlines():
            if line.startswith('>'):
                if n == 1:
                    outfile.write(line)
                if n > 1:
                    outfile.write('\n'+line)
            else:
                outfile.write(line.replace('\n', ''))
            n += 1
    os.remove(temp_out)
    return fasta_out


def generate_ali(alignments, best_oligo_template, residue_index_mapping, args):
    best_oligo_template_file = best_oligo_template+"_CHOIR_RenamedChainsTemplate"
    final_alignment = os.path.join(workdir, input_name+'_'+best_oligo_template+'_CHOIR_Alignment.ali')
    getseq = False
    alignment_dict = {}
    full_residue_mapping = {}
    # Parse individual GESAMT alignments and organize in a per-chain dictionary
    for fasta_alignment in alignments:
        getseq = False
        template = False
        chain = None
        entryseq_dict = {}
        for line in open(fasta_alignment, 'r').readlines():
            # Only record sequence if line above starts with >
            if getseq is True:
                getseq = False
                seq = line.replace('\n', '')
                # If this is the template, count leading and trailing gaps
                if template is True:
                    template = False
                    leading_gaps = 0
                    for r in seq:
                        if r == '-':
                            leading_gaps += 1
                        else:
                            break
                    trailing_gaps = 0
                    for r in seq[::-1]:
                        if r == '-':
                            trailing_gaps += 1
                        else:
                            break
                assert seq is not None, 'Sequence is None'
                assert seq != '', 'Sequence is empty'
                entryseq_dict[entry] = seq.upper()
                del seq
            # If it is an entry line, get details and expect sequence
            if line.startswith('>'):
                entry = line.split('>')[1].split('(')[0].split('.pdb')[0].replace('\n', '')
                # If entry is template, use chain as reference
                if entry == best_oligo_template_file:
                    chain = line.split('(')[1].split(')')[0]
                    template = True
                getseq = True

        # Remove leading and trailing gaps from the alignment for both template and query
        if trailing_gaps == 0:
            for entry, seq in entryseq_dict.items():
                entryseq_dict[entry] = leading_gaps*'-'+seq[leading_gaps:]
        else:
            for entry, seq in entryseq_dict.items():
                entryseq_dict[entry] = leading_gaps*'-'+seq[leading_gaps:-trailing_gaps]+trailing_gaps*'-'
        if residue_index_mapping is not None:
            full_residue_mapping[chain] = collections.OrderedDict()
            for res, i in residue_index_mapping.items():
                full_residue_mapping[chain][res] = i+leading_gaps
        else:
            full_residue_mapping[chain] = leading_gaps

        alignment_dict[chain] = entryseq_dict
        pctools.printv('Removed '+clrs['c']+str(leading_gaps)+clrs['n']+' leading gaps and '+clrs['c']+str(trailing_gaps)+clrs['n']+' trailing gaps from chain '+clrs['c']+chain+clrs['n']+' alignment.\n', verbosity)

    # If symmetry is desired, reduce all chains to match the size of the smallest
    if args.symmetry:
        max_leading_gaps = 0
        max_trailing_gaps = 0
        for chain, seqs in alignment_dict.items():
            for entry, seq in seqs.items():
                if entry == best_oligo_template_file:
                    leading_gaps = 0
                    for r in seq:
                        if r == '-':
                            leading_gaps += 1
                        else:
                            break
                    if leading_gaps > max_leading_gaps:
                        max_leading_gaps = leading_gaps
                    trailing_gaps = 0
                    for r in seq[::-1]:
                        if r == '-':
                            trailing_gaps += 1
                        else:
                            break
                    if trailing_gaps > max_trailing_gaps:
                        max_trailing_gaps = trailing_gaps
        pctools.printv('To cope with symmetry restraints, the modelled sequence will contain '+clrs['c']+str(max_leading_gaps)+clrs['n']+' leading gaps and '+clrs['c']+str(max_trailing_gaps)+clrs['n']+' trailing gaps'+clrs['n']+'.\n', verbosity)
        print(max_trailing_gaps)
        for chain, seqs in alignment_dict.items():
            if max_trailing_gaps == 0:
                seqs[entry] = max_leading_gaps*'-'+seqs[entry][max_leading_gaps:]
            else:
                seqs[entry] = max_leading_gaps*'-'+seqs[entry][max_leading_gaps:-max_trailing_gaps]+max_trailing_gaps*'-'


    # Find out first and last chains
    first_chain = sorted(alignment_dict)[0]
    last_chain = sorted(alignment_dict)[-1]

    # Create strings to write in alignment file
    alignment_string_dict = {}
    for entry in [input_name, best_oligo_template_file]:
        if entry == input_name:
            alignment_string_dict[entry] = ">P1;"+input_name+"\nsequence:"+input_name+":FIRST:"+first_chain+":LAST:"+last_chain+"::::\n"
        elif entry == best_oligo_template_file:
            alignment_string_dict[entry] = ">P1;"+best_oligo_template_file+".pdb\nstructureX:"+best_oligo_template_file+".pdb:FIRST:"+first_chain+":LAST:"+last_chain+"::::\n"
        for chain, entryseq in sorted(alignment_dict.items()):
            if chain == last_chain:
                alignment_string_dict[entry] += entryseq[entry]+'*\n'
            else:
                alignment_string_dict[entry] += entryseq[entry]+'/\n'

    # Write alignment file
    with open(final_alignment, 'w') as f:
        for entry, entrystring in alignment_string_dict.items():
            pctools.printv(entrystring, verbosity)
            f.write(entrystring)

    print('Modeller Alignment written to '+clrs['g']+os.path.basename(final_alignment)+clrs['n']+'\n')
    return final_alignment, full_residue_mapping


def create_genmodel(final_alignment, best_oligo_template, chains, nmodels, refine_level, symmetry):
    genmodel_file = os.path.join(workdir, input_name+'_'+best_oligo_template+'_CHOIR_Genmodel.py')
    for line in open(final_alignment, 'r').readlines():
        if line.startswith('sequence:'):
            sequence_name = line.split(':')[1]
        if line.startswith('structureX:'):
            template_structure = line.split(':')[1]
    with open(genmodel_file, 'w') as f:
        f.write('from modeller import *\nfrom modeller.automodel import *\n\nlog.verbose()\n\n\n')
        if symmetry is True:
            f.write('class symmodel(automodel):\n\tdef special_restraints(self,aln):\n')
            for chain_i, chain_j in itertools.combinations(chains, 2):
                f.write("\t\tself.restraints.symmetry.append(symmetry(selection(self.chains['"+chain_i+"']), selection(self.chains['"+chain_j+"']), 1.0))\n")
            f.write('\tdef user_after_single_model(self):\n\t\tself.restraints.symmetry.report(1.0)\n\n')
            f.write('\tdef special_patches(self,aln):\n\t\tself.rename_segments(segment_ids='+str(sorted(chains))+', renumber_residues='+str([1]*len(chains))+')\n\n')
            f.write('env = environ()\n\n')
            f.write("a = symmodel(env, alnfile='"+final_alignment+"', knowns=('"+template_structure+"'), sequence='"+sequence_name+"', assess_methods=(assess.DOPE, assess.GA341))\n")
        else:
            f.write('class renumbermodel(automodel):\n')
            f.write('\tdef special_patches(self,aln):\n\t\tself.rename_segments(segment_ids='+str(sorted(chains))+', renumber_residues='+str([1]*len(chains))+')\n\n')
            f.write('env = environ()\n\n')
            f.write("a = renumbermodel(env, alnfile='"+final_alignment+"', knowns=('"+template_structure+"'), sequence='"+sequence_name+"', assess_methods=(assess.DOPE, assess.GA341))\n")
        f.write("a.starting_model = 1\na.ending_model = "+str(nmodels)+"\n")
        if refine_level == 0:
            f.write("a.library_schedule = autosched.very_fast\na.max_var_iterations = 1000\n")
        elif refine_level == 1:
            f.write("a.library_schedule = autosched.fast\na.max_var_iterations = 1000\n")
        elif refine_level == 2:
            f.write("a.library_schedule = autosched.normal\na.max_var_iterations = 1000\n")
        elif refine_level == 3:
            f.write("a.library_schedule = autosched.slow\na.max_var_iterations = 1000\n")
        elif refine_level == 4:
            f.write("a.library_schedule = autosched.slow\na.max_var_iterations = 1000\na.md_level = refine.slow\n")
        f.write("a.make()")
    print('Modeller script written to '+clrs['g']+os.path.basename(genmodel_file)+clrs['n']+'\n')
    expected_models = [input_name+'.B9999000'+str(n)+'.pdb' for n in range(1, nmodels+1)]
    return genmodel_file, expected_models


def run_modeller(genmodel_file):
    print('Running '+clrs['b']+'MODELLER'+clrs['n']+' for '+clrs['y']+os.path.basename(genmodel_file)+clrs['n']+'\n')
    script_name = os.path.basename(genmodel_file).split('.py')[0]
    genmodel_log = os.path.join(workdir,script_name+'.log')
    temp = sys.stdout
    sys.stdout = open(genmodel_log, 'w')
    __import__(script_name)
    sys.stdout.close()
    sys.stdout = temp
    print('Done running '+clrs['b']+'MODELLER'+clrs['n']+'\n')


# Main Function
###############################################################################
def make_oligomer(input_file, largest_oligo_complexes, report, args, residue_index_mapping=None):
    # input_file = '/home/torres/work/protchoir/ml0002/ml0002_CHOIR_MonomerSequence.fasta'
    # with open('/home/torres/work/protchoir/ml0002/CHOIR_OligoComplexes.pickle', 'rb') as p:
    #     largest_oligo_complexes = pickle.load(p)
    # with open('/home/torres/work/protchoir/ml0002/CHOIR_args.pickle', 'rb') as p:
    #     args = pickle.load(p)
    # args.sequence_mode = True
    global workdir
    global input_name
    global verbosity
    verbosity = args.verbosity
    workdir = os.getcwd()
    symmetry = args.symmetry
    if args.sequence_mode is False:
        input_name = os.path.basename(input_file).split(".pdb")[0].replace('.', '_')
        candidate_qscores = {}
        # Select structurally best oligomeric template using GESAMT
        pctools.print_section(2, 'OLIGOMER ASSEMBLING')
        pctools.print_subsection('2[a]', 'Structural template selection')
        for hitchain, chains in largest_oligo_complexes.items():
            template, hit_chain = hitchain.split(':')
            middle_letters = template[1:3]
            template_file = os.path.join(pdb_homo_archive, middle_letters, template+".pdb.gz")
            sum_qscore = 0
            chain_n = 0
            for chain in chains:
                chain_n += 1
                qscore, rmsd, fasta_out = pctools.run_gesamt(template, template_file, input_name, input_file, chain, args, delfasta=True)
                sum_qscore += float(qscore)

            average_qscore = sum_qscore/chain_n
            report['hits'][hitchain]['qscore'] = round(average_qscore, 3)
            print('--\n\nAverage Q-Score for all candidate chains is '+clrs['c']+str(average_qscore)+clrs['n']+'\n')
            print('-------------------------------------------------------------------\n')
            candidate_qscores[hitchain] = average_qscore

        best_oligo_template = max(candidate_qscores.keys(), key=(lambda x: candidate_qscores[x]))
        if candidate_qscores[best_oligo_template] >= args.qscore_cutoff:
            print('Structurally, the best template is: '+clrs['y']+best_oligo_template+clrs['n']+'. Using that!\n')
            report['best_template'] = best_oligo_template.split(':')[0]
        else:
            print('No template had an average Q-score above cut-off of'+clrs['c']+str(args.qscore_cutoff)+clrs['n']+'\nTry lowering the cutoff or running in sequence mode.\n')
            pctools.print_sorry()
            quit()


    elif args.sequence_mode is True:
        if input_file.endswith('.pdb'):
            input_name = os.path.basename(input_file).split(".pdb")[0].replace('.', '_')
            input_file = os.path.join(workdir, input_name+'_CHOIR_MonomerSequence.fasta')

        elif input_file.endswith('_CHOIR_MonomerSequence.fasta'):
            input_name = os.path.basename(input_file).split("_CHOIR_MonomerSequence.fasta")[0]
        pctools.print_section(2, 'OLIGOMER ASSEMBLING - SEQUENCE MODE')
        print(clrs['y']+"Skipping section 2[a] - Structural template selection"+clrs['n']+"\n")
        best_oligo_template = list(largest_oligo_complexes)[0]
        report['best_template'] = best_oligo_template.split(':')[0]
    report['topology_figure'] = './'+best_oligo_template.replace(':', '_')+'_CHOIR_Topology.png'
    template_chains = largest_oligo_complexes[best_oligo_template]
    best_oligo_template_code = best_oligo_template.split(':')[0]
    clean_template_file = make_local_template(best_oligo_template_code)
    relevant_chains_file = extract_relevant_chains(clean_template_file, template_chains)
    report['template_figure'] = pctools.pymol_screenshot(relevant_chains_file, args)
    renamed_chains_file, chains_dict = rename_relevant_chains(relevant_chains_file)
    relevant_chains = [chains_dict[template_chain] for template_chain in template_chains]

    pctools.print_subsection('2[b]', 'Generating alignment')
    # Generate per chain alignment files
    alignment_files = []
    if args.sequence_mode is False:
        for original_chain, current_chain in chains_dict.items():
            qscore, rmsd, fasta_out = pctools.run_gesamt(best_oligo_template_code, renamed_chains_file, input_name, input_file, current_chain, args, delfasta=False)
            alignment_files.append(fasta_out)
    elif args.sequence_mode is True:
        for original_chain, current_chain in chains_dict.items():
            fasta_out = alignment_from_sequence(best_oligo_template_code, renamed_chains_file, input_file, current_chain)
            alignment_files.append(fasta_out)
    print('Alignment files:\n'+clrs['g']+('\n').join([os.path.basename(i) for i in alignment_files])+clrs['n'])

    # Generate final alignment which will be the input for Modeller
    final_alignment, full_residue_mapping = generate_ali(alignment_files, best_oligo_template_code, residue_index_mapping, args)
    pctools.print_subsection('2[c]', 'Generating models')
    genmodel_file, expected_models = create_genmodel(final_alignment, best_oligo_template_code, relevant_chains, args.models, args.refine_level, symmetry)
    run_modeller(genmodel_file)

    # Record list of oligomers built
    nmodels = 0
    built_oligomers = []
    for model in expected_models:
        built_oligomers.append(restore_chain_identifiers(model, chains_dict, full_residue_mapping))
        nmodels += 1
    print(clrs['b']+'ProtCHOIR'+clrs['n']+' built '+clrs['c']+str(nmodels)+clrs['n']+' model oligomers:')
    for model in built_oligomers:
        print(clrs['g']+model+clrs['n'])

    return best_oligo_template, built_oligomers, report