# Imports
###############################################################################
import os
import sys
import time
import gzip
import itertools
import subprocess
import collections
import Bio.PDB as bpp
import importlib.util
import ProtCHOIR.Toolbox as pctools
from Bio import SeqIO
from multiprocessing import Pool
from ProtCHOIR.Initialise import *
from Bio.SubsMat import MatrixInfo

# LICENSE
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




# Functions
###############################################################################
def make_local_template(best_oligo_template):
    middle_letters_best = best_oligo_template[1:3]
    if g_args.allow_monomers:
        best_template_file = os.path.join(pdb_archive, middle_letters_best, 'pdb'+best_oligo_template+".ent.gz")
        pdb_name, contents = pctools.parse_pdb_contents(best_template_file)
        is_nmr = pctools.is_nmr(contents)
        if is_nmr:
            print(clrs['r']+'\n\n Selected template '+best_oligo_template+' is an NMR structure \n Will try a a different candidate.\n\n'+clrs['n'])
            raise


    else:
        best_template_file = os.path.join(pdb_homo_archive, middle_letters_best, best_oligo_template+".pdb.gz")
    clean_template_file = os.path.join(workdir, best_oligo_template+"_CHOIR_CleanTemplate.pdb")
    pdb_name, structure, nchains = pctools.parse_any_structure(best_template_file)
    io.set_structure(structure)
    io.save(clean_template_file, pctools.SelectIfCA())
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


def alignment_from_sequence(current_chain):
    output = []
    input_basename = os.path.basename(g_input_file).split('_CHOIR_MonomerSequence.fasta')[0]
    output.append('Running '+clrs['b']+'MODELLER'+clrs['n']+' to align '+clrs['y']+input_basename+clrs['n']+' to '+clrs['y']+best_oligo_template_code+clrs['n']+' - Chain '+clrs['y']+current_chain+clrs['n']+'...')
    output_basename = input_basename+'_'+best_oligo_template_code+current_chain
    genali_file = os.path.join(workdir, output_basename+'_CHOIR_Genali.py')
    temp_out = os.path.join(workdir, output_basename+'_CHOIR_Align2Dtemp.fasta')
    fasta_out = os.path.join(workdir, output_basename+'_CHOIR_Align2D.fasta')
    with open(genali_file, 'w') as f:
        f.write('from modeller import *\nenv = environ()\naln = alignment(env)\n')
        f.write('\nlog.verbose()\n\n\n')
        f.write("mdl = model(env, file='"+renamed_chains_file+"', model_segment=('FIRST:"+current_chain+"','LAST:"+current_chain+"'))\n")
        f.write("aln.append_model(mdl, align_codes='"+os.path.basename(renamed_chains_file)+"("+current_chain+")', atom_files='"+renamed_chains_file+"')\n")
        f.write("aln.append(file='"+g_input_file+"', align_codes='"+input_basename+"',alignment_format='FASTA', remove_gaps=False)\n")
        f.write('aln.align2d()\n')
        f.write("aln.write(file='"+temp_out+"', alignment_format='FASTA')")
    alignment_log = genali_file.replace('.py', '.log')
    spec = importlib.util.spec_from_file_location("genmodel", genali_file)
    genali_module = importlib.util.module_from_spec(spec)
    temp = sys.stdout
    sys.stdout = open(alignment_log, 'w')
    spec.loader.exec_module(genali_module)
    #sys.stdout.close()
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
    output.append('DONE\n')
    return fasta_out, '\n'.join(output)


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


def create_genmodel(final_alignment, best_oligo_template, chains, args):
    genmodel_file = os.path.join(workdir, input_name+'_'+best_oligo_template+'_CHOIR_Genmodel.py')
    modelclass_file = 'CHOIR_ModelClass'
    for line in open(final_alignment, 'r').readlines():
        if line.startswith('sequence:'):
            sequence_name = line.split(':')[1]
        if line.startswith('structureX:'):
            template_structure = line.split(':')[1]
    with open(genmodel_file, 'w') as f:
        f.write('import os\nimport sys\nsys.path.append(os.path.dirname(os.path.abspath( __file__ )))\nfrom modeller import *\nfrom modeller.automodel import *\nfrom '+os.path.basename(modelclass_file)+' import CHOIRModel\n')
        if args.modeller_threads > 1:
            f.write('from modeller.parallel import *\n')
        f.write('\nlog.verbose()\n\n\n')
        if args.modeller_threads > 1:
            f.write('j = job()\nfor i in range('+str(args.modeller_threads)+'):\n\tj.append(local_slave())\n\n')
        f.write('env = environ()\n\n')
        f.write("a = CHOIRModel(env, alnfile='"+final_alignment+"', knowns=('"+template_structure+"'), sequence='"+sequence_name+"', assess_methods=(assess.DOPE, assess.GA341))\n")
        f.write("a.starting_model = 1\na.ending_model = "+str(args.models)+"\n")
        if args.refine_level == 0:
            f.write("a.library_schedule = autosched.very_fast\na.max_var_iterations = 1000\n")
        elif args.refine_level == 1:
            f.write("a.library_schedule = autosched.fast\na.max_var_iterations = 1000\n")
        elif args.refine_level == 2:
            f.write("a.library_schedule = autosched.normal\na.max_var_iterations = 1000\n")
        elif args.refine_level == 3:
            f.write("a.library_schedule = autosched.slow\na.max_var_iterations = 1000\n")
        elif args.refine_level == 4:
            f.write("a.library_schedule = autosched.slow\na.max_var_iterations = 1000\na.md_level = refine.slow\n")
        if args.repeat_opt > 0:
            f.write("a.repeat_optimization = "+str(args.repeat_opt)+"\n")
        if args.modeller_threads > 1:
            f.write("a.use_parallel_job(j)\n")
        f.write("a.make()")

    with open(modelclass_file+'.py', 'w') as f:
        f.write('from modeller import *\nfrom modeller.automodel import *\n')
        if args.symmetry is True:
            f.write('class CHOIRModel(automodel):\n\tdef special_restraints(self,aln):\n')
            for chain_i, chain_j in itertools.combinations(chains, 2):
                f.write("\t\tself.restraints.symmetry.append(symmetry(selection(self.chains['"+chain_i+"']), selection(self.chains['"+chain_j+"']), 1.0))\n")
            f.write('\tdef user_after_single_model(self):\n\t\tself.restraints.symmetry.report(1.0)\n\n')
            f.write('\tdef special_patches(self,aln):\n\t\tself.rename_segments(segment_ids='+str(sorted(chains))+', renumber_residues='+str([1]*len(chains))+')\n\n')
        else:
            f.write('class CHOIRModel(automodel):\n')
            f.write('\tdef special_patches(self,aln):\n\t\tself.rename_segments(segment_ids='+str(sorted(chains))+', renumber_residues='+str([1]*len(chains))+')\n\n')



    print('Modeller script written to '+clrs['g']+os.path.basename(genmodel_file)+clrs['n']+'\n')
    print('Modeller class written to '+clrs['g']+os.path.basename(modelclass_file)+'.py'+clrs['n']+'\n')
    expected_models = [input_name+'.B9999'+str('{0:04d}'.format(n))+'.pdb' for n in range(1, args.models+1)]
    return genmodel_file, expected_models


def run_modeller(genmodel_file):
    print('Running '+clrs['b']+'MODELLER'+clrs['n']+' for '+clrs['y']+os.path.basename(genmodel_file)+clrs['n']+'\n')
    script_name = os.path.basename(genmodel_file).split('.py')[0]
    genmodel_log = os.path.join(workdir,script_name+'.log')
    spec = importlib.util.spec_from_file_location("genmodel", genmodel_file)
    genmodel_module = importlib.util.module_from_spec(spec)
    temp = sys.stdout
    sys.stdout = open(genmodel_log, 'w')
    spec.loader.exec_module(genmodel_module)
    sys.stdout.close()
    sys.stdout = temp
    print('Done running '+clrs['b']+'MODELLER'+clrs['n']+'\n')


def analyse_largest_complexes(item):
    output = []
    hitchain, chains = item
    template, hit_chain = hitchain.split(':')
    middle_letters = template[1:3]
    template_file = os.path.join(pdb_homo_archive, middle_letters, template+".pdb.gz")
    sum_qscore = 0
    chain_n = 0
    for chain in chains:
        chain_n += 1
        qscore, rmsd, fasta_out, gesamt_output = pctools.run_gesamt(template, template_file, input_name, g_input_file, chain, g_args)
        sum_qscore += float(qscore)
        output.append(gesamt_output)

    average_qscore = sum_qscore/chain_n
    output.append('--\n\nAverage Q-Score for all candidate chains is '+clrs['c']+str(average_qscore)+clrs['n']+'\n')
    output.append('-------------------------------------------------------------------\n')

    return hitchain, average_qscore, '\n'.join(output)


def run_gesamt_parallel(chain):
    output = []
    output.append('Running '+clrs['b']+'GESAMT'+clrs['n']+' to align '+clrs['y']+input_name+clrs['n']+' to '+clrs['y']+best_oligo_template_code+clrs['n']+' - Chain '+clrs['y']+chain+clrs['n'])
    fasta_out = input_name+"_"+best_oligo_template_code+chain+'_CHOIR_Gesamt.fasta'
    gesamtcmd = [gesamt_exe, renamed_chains_file, '-s', chain, g_input_file, '-a', fasta_out]
    if g_args.verbosity == 1:
        output.append(clrs['b']+'GESAMT'+clrs['n']+' command line: '+' '.join(gesamtcmd))
    gesout = subprocess.check_output(gesamtcmd).decode('utf-8').split('\n')
    time.sleep(1)
    for line in gesout:
        if line.startswith(' Q-score          :'):
            qscore = line.split(':')[1].strip()
            output.append(line.strip())
        if line.startswith(' RMSD             :'):
            rmsd = line.split(':')[1].strip()
            output.append(line.strip())
        if line.startswith(' Aligned residues :'):
            al_res = line.split(':')[1].strip()
            output.append(line.strip())
        if line.startswith(' Sequence Id:     :'):
            id = line.split(':')[1].strip()
            output.append(line.strip())

    if not os.path.isfile(fasta_out):
        output.append(clrs['r']+'GESAMT FAILED'+clrs['n'])
        qscore = 0
        rmsd = 0
        fasta_out = 'NONE'
    output.append('Done running '+clrs['b']+'GESAMT'+clrs['n']+'. Alignment written to '+clrs['g']+os.path.basename(fasta_out)+clrs['n']+'\n')

    return qscore, rmsd, fasta_out, '\n'.join(output)

def score_match(pair, matrix):
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]+4
    else:
        return matrix[pair]+4

def score_pairwise(seq1, seq2, matrix, gap_s, gap_e):
    score = 0
    gap = False
    ipos = 0
    fpos = 30
    nwindows = -(-len(seq1)//30)
    pctools.printv('Number of 30-residue segments: '+str(nwindows), g_args.verbosity)
    wscores = []
    for window in range(nwindows):
        wscore = 0
        if fpos > len(seq1):
            fpos = len(seq1)
        pctools.printv(str(ipos+1)+' '+seq1[ipos:fpos]+' '+str(fpos), g_args.verbosity)
        pctools.printv(str(ipos+1)+' '+seq2[ipos:fpos]+' '+str(fpos), g_args.verbosity)
        for i in range(len(seq1))[ipos:fpos]:
            pair = (seq1[i], seq2[i])
            if not gap:
                if pair == ('-', '-'):
                    score += 4
                    wscore += 4
                elif '-' in pair:
                    gap = True
                    score += gap_s
                    wscore += gap_s
                else:
                    score += score_match(pair, matrix)
                    wscore += score_match(pair, matrix)
            else:
                if '-' not in pair:
                    gap = False
                    score += score_match(pair, matrix)
                    wscore += score_match(pair, matrix)
                else:
                    score += gap_e
                    wscore += gap_e

        ipos += 30
        fpos += 30
        pctools.printv('Segment score: '+str(wscore), g_args.verbosity)
        wscores.append(wscore)

    return score, wscores


def score_alignment(alignment_file):
    print(clrs['b']+'SCORING ALIGNMENT'+clrs['n']+' in '+clrs['y']+os.path.basename(alignment_file)+clrs['n']+'\n')
    sequences = list(SeqIO.parse(alignment_file, "pir"))
    query_chains = str(sequences[0].seq).split('/')
    template_chains = str(sequences[1].seq).split('/')
    trimmed_query_chains = []
    trimmed_template_chains = []
    for query_chain, template_chain in zip(query_chains, template_chains):

        leading_gaps = 0
        for r in query_chain:
            if r == '-':
                leading_gaps += 1
            else:
                break
        trailing_gaps = 0
        for r in query_chain[::-1]:
            if r == '-':
                trailing_gaps += 1
            else:
                break

        if trailing_gaps == 0:
            trimmed_query_chains.append(query_chain[leading_gaps:])
            trimmed_template_chains.append(template_chain[leading_gaps:])
        else:
            trimmed_query_chains.append(query_chain[leading_gaps:-trailing_gaps])
            trimmed_template_chains.append(template_chain[leading_gaps:-trailing_gaps])

    relative_wscores = []
    relative_scores = []
    for q_chain, t_chain in zip(trimmed_query_chains, trimmed_template_chains):
        pctools.printv('\nCalculating '+clrs['y']+'maximum scores'+clrs['n']+' for chain segments:', g_args.verbosity)
        max_score, max_wscores = score_pairwise(t_chain, t_chain, MatrixInfo.blosum62, 0, 0)
        pctools.printv('\nCalculating '+clrs['y']+'actual scores'+clrs['n']+' for chain segments:', g_args.verbosity)
        score, wscores = score_pairwise(q_chain, t_chain, MatrixInfo.blosum62, 0, 0)
        relative_scores.append(round(score*100/max_score, 2))

        for max_wscore, wscore in zip(max_wscores, wscores):
            if max_wscore != 0:
                relative_wscore = round(wscore*100/max_wscore, 2)
            else:
                relative_wscore = 100
            relative_wscores.append(relative_wscore)

    relative_score = sum(relative_scores) / len(relative_scores)
    string = ''
    for relative_wscore in relative_wscores:
        if relative_wscore > g_args.similarity_cutoff:
            color = 'g'
        else:
            color = 'r'
        if string == '':
            string += (clrs[color]+str(relative_wscore)+clrs['n'])
        else:
            string += (' ~ '+clrs[color]+str(relative_wscore)+clrs['n'])
    print('\nRelative score per 30-res segment: '+string+clrs['n'])
    return relative_score, relative_wscores, len(query_chains)



# Main Function
###############################################################################
def make_oligomer(input_file, largest_oligo_complexes, report, args, residue_index_mapping=None):
    global workdir
    global input_name
    global verbosity
    global g_input_file
    global g_args
    global best_oligo_template_code
    global renamed_chains_file
    g_input_file = input_file
    g_args = args
    verbosity = args.verbosity
    workdir = os.getcwd()
    symmetry = args.symmetry

    # Subsection 2[a] #######################################################################
    if args.sequence_mode is False:
        input_name = os.path.basename(input_file).split(".pdb")[0].replace('.', '_')
        candidate_qscores = {}
        # Select structurally best oligomeric template using GESAMT
        pctools.print_section(2, 'OLIGOMER ASSEMBLING')
        pctools.print_subsection('2[a]', 'Structural template selection')
        if args.multiprocess is True:
            p = Pool()
            for hitchain, average_qscore, output in p.map_async(analyse_largest_complexes, largest_oligo_complexes.items()).get():
                candidate_qscores[hitchain] = average_qscore
                report['hits'][hitchain]['qscore'] = round(average_qscore, 3)
                print(output)
            p.close()
            p.join()
        else:
            for item in largest_oligo_complexes.items():
                hitchain, average_qscore, output = analyse_largest_complexes(item)
                candidate_qscores[hitchain] = average_qscore
                report['hits'][hitchain]['qscore'] = round(average_qscore, 3)
                print(output)

        best_oligo_template = max(candidate_qscores.keys(), key=(lambda x: candidate_qscores[x]))
        if candidate_qscores[best_oligo_template] >= args.qscore_cutoff:
            print('Structurally, the best template is: '+clrs['y']+best_oligo_template+clrs['n']+'. Using that!\n')
            report['best_template'] = best_oligo_template.split(':')[0]
            report['best_id'] = report['hits'][best_oligo_template]['id']
            report['best_cov'] = report['hits'][best_oligo_template]['coverage']
            report['best_qscore'] = report['hits'][best_oligo_template]['qscore']
            report['best_nchains'] = report['hits'][best_oligo_template]['final_homo_chains']
        else:
            print('No template had an average Q-score above cut-off of '+clrs['c']+str(args.qscore_cutoff)+clrs['n']+'\nTry lowering the cutoff or running in sequence mode.\n')
            report['exit'] = '4'
            return None, None, report
        report['topology_figure'] = './'+best_oligo_template.replace(':', '_')+'_CHOIR_Topology.png'
        template_chains = largest_oligo_complexes[best_oligo_template]
        best_oligo_template_code = best_oligo_template.split(':')[0]
        clean_template_file = make_local_template(best_oligo_template_code)


    elif args.sequence_mode is True:
        if input_file.endswith('.pdb'):
            input_name = os.path.basename(input_file).split(".pdb")[0].replace('.', '_')
            input_file = os.path.join(workdir, input_name+'_CHOIR_MonomerSequence.fasta')
            g_input_file = input_file

        elif input_file.endswith('_CHOIR_MonomerSequence.fasta'):
            input_name = os.path.basename(input_file).split("_CHOIR_MonomerSequence.fasta")[0]

        pctools.print_section(2, 'OLIGOMER ASSEMBLING - SEQUENCE MODE')
        print(clrs['y']+"Skipping section 2[a] - Structural template selection"+clrs['n']+"\n")
        attempt = 0
        while attempt < len(largest_oligo_complexes):
            try:
                best_oligo_template = list(largest_oligo_complexes)[attempt]
                report['best_template'] = best_oligo_template.split(':')[0]
                report['best_id'] = report['hits'][best_oligo_template]['id']
                report['best_cov'] = report['hits'][best_oligo_template]['coverage']
                report['best_qscore'] = 'NA'
                report['best_nchains'] = report['hits'][best_oligo_template]['final_homo_chains']
                report['topology_figure'] = './'+best_oligo_template.replace(':', '_')+'_CHOIR_Topology.png'
                template_chains = largest_oligo_complexes[best_oligo_template]
                best_oligo_template_code = best_oligo_template.split(':')[0]
                clean_template_file = make_local_template(best_oligo_template_code)
                break
            except:
                attempt += 1
                if attempt < len(largest_oligo_complexes):
                    print('Attempt '+str(attempt)+' failed, trying a differente template candidate.')
                if attempt == len(largest_oligo_complexes):
                    print('Failed to find templates in local databases.')
                    report['exit'] = '5'
                    return None, None, report


    relevant_chains_file = extract_relevant_chains(clean_template_file, template_chains)
    if args.generate_report is True:
        report['template_figure'], pymol_output = pctools.pymol_screenshot(relevant_chains_file, args)
        print(pymol_output)
    renamed_chains_file, chains_dict = rename_relevant_chains(relevant_chains_file)
    relevant_chains = [chains_dict[template_chain] for template_chain in template_chains]

    # Subsection 2[b] #######################################################################
    pctools.print_subsection('2[b]', 'Generating alignment')
    # Generate per chain alignment files
    alignment_files = []
    if args.sequence_mode is False:
        if args.multiprocess is True:
            p = Pool()
            for qscore, rmsd, fasta_out, gesamt_output in p.map_async(run_gesamt_parallel, chains_dict.values()).get():
                alignment_files.append(fasta_out)
                print(gesamt_output)
            p.close()
            p.join()
        else:
            for chain in chains_dict.values():
                qscore, rmsd, fasta_out, gesamt_output = run_gesamt_parallel(chain)
                alignment_files.append(fasta_out)
                print(gesamt_output)

    elif args.sequence_mode is True:
        if args.multiprocess is True:
            p = Pool()
            for fasta_out, output in p.map_async(alignment_from_sequence, chains_dict.values()).get():
                alignment_files.append(fasta_out)
                print(output)
        else:
            for current_chain in chains_dict.values():
                fasta_out, output = alignment_from_sequence(current_chain)
                alignment_files.append(fasta_out)
                print(output)
    print('Alignment files:\n'+clrs['g']+('\n').join([os.path.basename(i) for i in alignment_files])+clrs['n'])


    # Generate final alignment which will be the input for Modeller
    final_alignment, full_residue_mapping = generate_ali(alignment_files, best_oligo_template_code, residue_index_mapping, args)
    # Score said alignment and enforce treshold
    report['relative_alignment_score'], relative_wscores, nchains = score_alignment(final_alignment)
    print('\nFinal average relative score for alignment: '+str(round(report['relative_alignment_score'], 2))+'%')
    bad_streches = 0
    for wscore in relative_wscores:
        if wscore < args.similarity_cutoff:
            bad_streches += 1
    if bad_streches >= args.bad_streches*nchains:
        if args.sequence_mode is True:
            print('\nThe alignment score was unacceptable for '+clrs['r']+str(bad_streches)+clrs['n']+' 30-res segments of the protein complex.\nTry running the default (structure) mode.\n')
        else:
            print('\nThe alignment score was unacceptable for '+clrs['r']+str(bad_streches)+clrs['n']+' 30-res segments of the protein complex.\nTry increasing the number of candidate templates or tweaking the similarity cut-offs.\n')
        report['exit'] = '6'
        return None, None, report

    # Subsection 2[c] #######################################################################
    pctools.print_subsection('2[c]', 'Generating models')
    genmodel_file, expected_models = create_genmodel(final_alignment, best_oligo_template_code, relevant_chains, args)
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
