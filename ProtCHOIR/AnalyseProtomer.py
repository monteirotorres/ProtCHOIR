# Imports
###############################################################################
import os
import re
import math
import pickle
import subprocess
import matplotlib
import collections
import numpy as np
import pandas as pd
matplotlib.use('Agg')
import textwrap as tw
import networkx as nx
from Bio import SeqIO
import matplotlib.pyplot as plt
from multiprocessing import Pool
from ProtCHOIR.Initialise import *
import xml.etree.ElementTree as et
from matplotlib.lines import Line2D
import ProtCHOIR.Toolbox as pctools
from progressbar import progressbar as pg
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

# Dictionaries
###############################################################################

# Global Variables
###############################################################################

# Classes
###############################################################################

# Functions
###############################################################################
def parse_templated_model(sequence, modelfile):
    print(clrs['g']+"Templated model! New homologue search not necessary!\n"+clrs['n'])
    pattern = 'REMARK   6 TEMPLATE:'
    templates_ids = {}
    with open(modelfile, 'r') as f:
        for line in f:
            if re.search(pattern, line):
                template = line.split()[3]
                templated_pid = line.split()[-1].replace('%', '')
                if len(template) == 5:
                    template_pdbcode = template[:4].lower()
                    templated_chain = template[4]
                elif len(template) == 7:
                    if template.startswith('d'):
                        template_pdbcode = template[1:5].lower()
                        templated_chain = template[5]
                    else:
                        template_pdbcode = template[:4].lower()
                        templated_chain = template[4]
                template_chain, pid = pctools.find_most_similar_chain(sequence, template_pdbcode)
                if template_chain is None:
                    print(clrs['r']+'Hit '+str(template_pdbcode)+' not found in oligomeric database. Disregarding it.'+clrs['n']+'\n')
                    continue
                else:
                    templates_ids[template_pdbcode+":"+template_chain] = float(pid)
                    print('Templated-determined hit: '+clrs['p']+template_pdbcode+clrs['n']+' Chain: '+clrs['y']+templated_chain+clrs['n']+' id: '+clrs['c']+str(templated_pid)+clrs['n'])
                    print('Most identical chain found in CHOIRdb oligomer: '+clrs['y']+template_chain+clrs['n']+' id: '+clrs['c']+str(pid)+clrs['n'])
    return templates_ids

def write_fasta(sequence):
    # Write protomer sequence to file
    fasta_file = os.path.join(workdir, pdb_name+'_CHOIR_MonomerSequence.fasta')
    with open(fasta_file, 'w') as f:
        f.write('>'+pdb_name+'\n')
        wrapped_seq = "\n".join(tw.wrap(sequence))
        f.write(wrapped_seq+'\n')
    return fasta_file

def blast_protomer(fasta_file, database, nhits, nit, nthreads, params, verbosity):
    print('\nRunning '+clrs['b']+'PSI-BLAST'+clrs['n']+' ( '+clrs['c']+os.path.basename(fasta_file)+clrs['n']+' x '+clrs['c']+os.path.basename(database)+clrs['n']+' )')
    matrix, gapopen, gapextend = params.split('-')
    blast_cmd = [psiblast_exe,
                 '-query', fasta_file,
                 '-db', database,
                 '-gapopen', gapopen,
                 '-gapextend', gapextend,
                 '-matrix', matrix,
                 '-word_size', '3',
                 '-num_threads', str(nthreads),
                 '-inclusion_ethresh', '0.005',
                 '-num_iterations', str(nit),
                 '-outfmt', '5',
                 '-num_alignments', '500',
                 '-comp_based_stats', '1',
                 '-max_hsps', '1']
    if database == uniref50:
        blast_cmd += ['-qcov_hsp_perc', '75']
    pctools.printv(clrs['b']+'BLAST'+clrs['n']+' command line: '+' '.join(blast_cmd), verbosity)
    blast_res = subprocess.Popen(blast_cmd, stdout=subprocess.PIPE)
    blast_xml = blast_res.stdout.read().decode()
    with open(os.path.join(workdir, os.path.basename(database)+'-tree.xml'), 'w') as f:
        f.write(blast_xml)
    blast_tree = et.fromstring(blast_xml)
    blast_list = {}
    for iteration in blast_tree:
        if iteration.tag == 'BlastOutput_query-len':
            query_length = iteration.text
        if iteration.tag == 'BlastOutput_iterations':
            for hits in iteration[-1]:
                if hits.tag == 'Iteration_iter-num':
                    last_iter = hits.text
                if hits.tag == 'Iteration_hits':
                    for hit in hits:
                        for hsps in hit:
                            if database == homoblast or database == monoblast or database == heteroblast:
                                if hsps.tag == 'Hit_def':
                                    code = ':'.join(hsps.text.split('|')[-2:])
                            elif database == uniref50:
                                if hsps.tag == 'Hit_id':
                                    code = hsps.text.split(' ')[0]
                            elif database == seqres:
                                if hsps.tag == 'Hit_id':
                                    code = ':'.join(hsps.text.split('|')[1:3])
                            if hsps.tag == 'Hit_hsps':
                                for hsp in hsps:
                                    for data in hsp:
                                        if data.tag == 'Hsp_score':
                                            score = data.text
                                        if data.tag == 'Hsp_identity':
                                            nid = data.text
                                        if data.tag == 'Hsp_align-len':
                                            ali_length = data.text
                                        if data.tag == 'Hsp_query-from':
                                            from_res = data.text
                                        if data.tag == 'Hsp_query-to':
                                            to_res = data.text
                        cov_length = int(to_res)-int(from_res)
                        cov = round((int(cov_length)/int(query_length))*100, 2)
                        pid = round((int(nid)/int(ali_length))*100, 2)
                        blast_list[code] = [int(score), float(pid), float(cov)]
    print('Last '+clrs['b']+'PSI-BLAST'+clrs['n']+' iteration was: '+clrs['c']+last_iter+clrs['n'])
    topn = collections.OrderedDict()
    if blast_list:
        for i in range(1,nhits+1):
            if blast_list:
                top = max(blast_list.keys(), key=(lambda x: blast_list[x][0]))
                topn[top]=blast_list[top]
                del blast_list[top]
        return topn
    else:
        return False

def generate_msa_input(topn, fasta_file, verbosity):
    multi_fasta = os.path.join(workdir, pdb_name+'_CHOIR_UniRef50Hits.fasta')
    if os.path.isfile(multi_fasta):
        os.remove(multi_fasta)
    for hit in topn:
        pctools.printv(clrs['p']+hit+clrs['y']+' score:'+clrs['c']+str(topn[hit][0])+clrs['y']+' id:'+clrs['c']+str(topn[hit][1])+clrs['y']+' cov:'+clrs['c']+str(topn[hit][2])+clrs['n'], verbosity)
        blastdbcmd = [blastdbcmd_exe, '-entry', hit, '-db', uniref50]
        pctools.printv(clrs['b']+'Blastdbcmd'+clrs['n']+' command line: '+' '.join(blastdbcmd), verbosity)
        blastdb_out = subprocess.Popen(blastdbcmd, stdout=subprocess.PIPE)
        blastdb_res = blastdb_out.stdout.read().decode()
        with open(multi_fasta, 'a') as f:
            f.write(blastdb_res+'\n')
    with open(fasta_file, 'r') as fin:
        with open(multi_fasta, 'a') as fout:
            for line in fin:
                fout.write(line)
    return multi_fasta


def run_mafft(multi_fasta, args):
    if args.force_single_core is True:
        mafft_threads = '1'
    else:
        mafft_threads = str(args.available_cores)
    print('\nRunning '+clrs['b']+'MAFFT'+clrs['n']+'...')
    msa_file = os.path.join(workdir, pdb_name+'_CHOIR_UniRef50Hits.msa')
    mafftcmd = [mafft_exe, '--localpair', '--maxiterate', '1000', '--quiet', '--thread', mafft_threads, '--anysymbol', multi_fasta]
    pctools.printv(clrs['b']+'MAFFT'+clrs['n']+' command line: '+' '.join(mafftcmd), args.verbosity)
    with open(msa_file, 'w') as f:
        subprocess.run(mafftcmd, stdout=f, check=True)
    print('Done running '+clrs['b']+'MAFFT'+clrs['n']+'. MSA file written to '+clrs['g']+os.path.basename(msa_file)+clrs['n'])
    return msa_file


def parse_msa(msa_file):
    msa_dict = {}
    entry = None
    with open(msa_file, 'r') as fin:
        for line in fin:
            if line.startswith('>'):
                if entry:
                    msa_dict[entry]=seq
                entry = line.split('>')[1].replace('\n','')
                seq = ''
            else:
                seqline = line.replace('\n','')
                seq += seqline
        msa_dict[entry]=seq
    return msa_dict

def trim_msa(msa_file):
    trimmed_msa = os.path.join(workdir, pdb_name+'_CHOIR_UniRef50HitsTrim.msa')
    if os.path.isfile(trimmed_msa):
        os.remove(trimmed_msa)
    # Create initial dictionary
    msa_dict = parse_msa(msa_file)
    # Get start and end positions of pdb
    for entry, seq in msa_dict.items():
        if entry == pdb_name:
            start = 0
            end = len(seq)
            for i in seq:
                if i != '-':
                    break
                start += 1
            for i in seq[::-1]:
                if i != '-':
                    break
                end -= 1
    # Create trimmed dictionary
    msa_dict_trimmed = {}
    for entry, seq in msa_dict.items():
        msa_dict_trimmed[entry] = seq[start:end]

    # Write trimmed MSA file
    with open(trimmed_msa, 'a') as fout:
        for entry, seq in msa_dict_trimmed.items():
            wrapped_seq = "\n".join(tw.wrap(seq,break_on_hyphens=False))
            fasta_entry = '>'+entry+'\n'+wrapped_seq+'\n\n'
            fout.write(fasta_entry)
    print('Trimmed MSA file written to '+clrs['g']+os.path.basename(trimmed_msa)+clrs['n'])
    return trimmed_msa

def shannon_entropy(msa_dict_trim,  surface_residues):
    # Turn sequences into arrays
    for entry, seq in msa_dict_trim.items():
        msa_dict_trim[entry] = np.array(list(seq))
    # Make dataframe
    df = pd.DataFrame.from_dict(msa_dict_trim, orient='index')
    # Remove gaps from reference sequence (the protomer sequence)
    df = df.loc[: , df.loc[pdb_name].ne('-')]
    # Fix column numbering
    x = []
    if surface_residues:
        for res in surface_residues.keys():
            x.append(int(res[3:]))

    else:
        for res in range(1,len(df.columns)+1):
            x.append(res)
    df.columns = x
    # Calculate entropies
    entropies = collections.OrderedDict()
    for column_number in df.columns:
        column = np.array(df[column_number])
        entropy = 0
        frequencies = {}
        for aa in set(column):
            count = 0
            for residue in column:
                if aa == residue:
                    count+=1
            entropy += -(count / len(column))*math.log2(count / len(column))
        entropies[column_number+1] = entropy
    return entropies

def relative_entropy(msa_dict_trim, surface_residues):
    # Turn sequences into arrays
    for entry, seq in msa_dict_trim.items():
        msa_dict_trim[entry] = np.array(list(seq))
    # Make dataframe
    df = pd.DataFrame.from_dict(msa_dict_trim, orient='index')
    # Remove gaps from reference sequence (the protomer sequence)
    df = df.loc[: , df.loc[pdb_name].ne('-')]
    # Fix column numbering
    x = []
    if surface_residues:
        for res in surface_residues.keys():
            x.append(int(res[3:]))

    else:
        for res in range(1,len(df.columns)+1):
            x.append(res)
    df.columns = x

    # Calculate entropies
    entropies = collections.OrderedDict()
    for column_number in df.columns:
        column = np.array(df[column_number])
        entropy = 0
        frequencies = {}
        for aa in set(column):
            if aa != '-' and aa != 'X':
                count = 0
                for residue in column:
                    if aa == residue:
                        count+=1
                # Use tryptophan background distribution if residue is non-standard
                if aa in aa_bgf:
                    entropy += (count / len(column))*math.log2((count / len(column))/aa_bgf[aa])
                else:
                    entropy += (count / len(column))*math.log2((count / len(column))/aa_bgf['W'])
        entropies[column_number] = entropy
    return entropies


def calc_z_scores(entropies):
    #Entropy Z Scores
    average = sum(entropies.values()) / len(entropies)
    sum_dev2 = 0
    for col, entropy in entropies.items():
        dev2 = (entropy - average)**2
        sum_dev2 += dev2
    stdev = math.sqrt(sum_dev2/len(entropies))
    z_entropies = collections.OrderedDict()
    for column_number, entropy in entropies.items():
        z = (entropy - average)/stdev
        z_entropies[column_number] = z
    return z_entropies


def map_conservation(structure, z_entropies):
    for atom in structure.get_atoms():
        for res, entropy in z_entropies.items():
            if int(atom.get_parent().id[1]) == int(res):
                atom.bfactor = entropy
    io.set_structure(structure)
    outfile = os.path.join(workdir, pdb_name+'_CHOIR_Conservation.pdb')
    io.save(outfile)
    return outfile


def confirm_homo_state(hit_code, candidate_chains, interfaces_list):
    output = []
    # Create cluster dictionary and flag chains with cluster number
    clustered_nodes = {}
    for node in candidate_chains:
        clustered_nodes[node] = 0

    # Define list of edges from interfaces calculated by PISA
    edges = []
    for interface in interfaces_list:
        edge = tuple(interface['chains'])
        edges.append(edge)

    # Use edges to determine clusters
    n = 1
    while any(cluster == 0 for chain, cluster in clustered_nodes.items()):
        for chain, cluster in clustered_nodes.items():
            if cluster == 0:
                clustered_nodes[chain] = n
                break
        marked = True
        while marked is True:
            marked = False
            for edge in edges:
                for chain in edge:
                    if clustered_nodes[chain] == n:
                        for chain in edge:
                            if clustered_nodes[chain] != n:
                                clustered_nodes[chain] = n
                                marked = True

        n += 1

    # Conclusion
    clusters = set([cluster for chain, cluster in clustered_nodes.items()])
    nclusters = len(clusters)
    cluster_dict = collections.OrderedDict()
    # List clusters and choose largest
    current_largest = 0
    for cluster_n in clusters:
        cluster_chains = []
        for chain, cluster in clustered_nodes.items():
            if cluster == cluster_n:
                cluster_chains.append(chain)
        if len(cluster_chains) > current_largest:
            current_largest = len(cluster_chains)
            largest_cluster = cluster_n
        cluster_dict[cluster_n] = cluster_chains

    if nclusters == 1:
        chain_string = (', ').join(cluster_dict[cluster_n])
        output.append('All candidate chains ('+clrs['y']+chain_string+clrs['n']+') are in close contact.')
        output.append('If used as template, would yield a '+clrs['y']+' HOMO-'+oligo_dict[len(candidate_chains)]+clrs['n']+' structure.')
    elif nclusters > 1:
        output.append('Not all candidate chains in structure '+clrs['c']+hit_code+clrs['n']+' are in close contact...')
        output.append('Instead, they are divided in '+clrs['r']+str(nclusters)+' clusters:'+clrs['n'])

        for cluster_n, chains in cluster_dict.items():
            if len(chains) > 1:
                chain_string = (', ').join(chains)
                output.append('Cluster '+clrs['y']+str(cluster_n)+clrs['n']+' contains chains: '+clrs['y']+chain_string+clrs['n']+'.')
                monomeric_cluster = False
            elif len(chains) == 1:
                chain_string = (', ').join(chains)
                output.append('Cluster '+clrs['y']+str(cluster_n)+clrs['n']+' contains a single chain: '+clrs['y']+chain_string+clrs['n']+'.')
                monomeric_cluster = True
        if monomeric_cluster is False:
            output.append('Largest cluster (Cluster No. '+clrs['y']+str(largest_cluster)+clrs['n']+'), would yield a '+clrs['y']+'HOMO-'+oligo_dict[len(cluster_dict[largest_cluster])]+clrs['n']+' structure.')
        elif monomeric_cluster is True:
            output.append('Largest cluster (Cluster No. '+clrs['y']+str(largest_cluster)+clrs['n']+'), would yield a '+clrs['y']+oligo_dict[len(cluster_dict[largest_cluster])]+clrs['n']+' structure.')

    return cluster_dict, largest_cluster, '\n'.join(output)


def plot_topology(complex_name, interfaces_list, cluster_dict):
    '''

    '''
    complex_name, chain_name = complex_name.split(':')
    G = nx.Graph()
    labels = {}
    for cluster, chains in cluster_dict.items():
        for chain in sorted(chains):
            G.add_nodes_from(chain)
            labels[chain] = chain

    # Fetch maximum and minimum areas of interaction from interfaces list
    if interfaces_list:
        maxw = max([abs(interface['interface area']) for interface in interfaces_list])
        minw = min([abs(interface['interface area']) for interface in interfaces_list])
    else:
        minw = 0
        maxw = 0

    # Initalize plots
    p, ax = plt.subplots(figsize=(8, 8))
    if interfaces_list:
        ax.title.set_text('Edge Weight: '+str(round(minw, 2))+' Å² - '+str(round(maxw, 2))+' Å²')
        ax.title.set_fontweight('bold')
        ax.title.set_fontsize(20)
    plt.axis('off')
    plt.suptitle('Topology of '+complex_name, fontsize=30, fontweight='bold')
    pos = nx.drawing.circular_layout(G)
    nodes = nx.draw_networkx_nodes(G, pos, node_size=3000, alpha=0.9, linewidths=4, node_color='skyblue' )
    nodes.set_edgecolor('k')
    nx.draw_networkx_labels(G,pos,labels,font_size=25, font_weight="bold",)

    # # Create edges and set weights
    if interfaces_list:
        for interface in interfaces_list:
            edge = tuple(interface['chains'])
            w = abs(interface['interface area'])**(1/2)/abs(maxw)**(1/2)*15
            G.add_edges_from([edge])
            nx.draw_networkx_edges(G,pos,edgelist=[edge],alpha=0.9,width=w)

    # Save figure
    outpath = os.path.join(workdir, complex_name+'_'+chain_name+"_CHOIR_Topology.png")
    plt.savefig(outpath, dpi=600)
    # Close figure
    plt.close()

    # Be verbose about it
    print('Topology plot for '+complex_name+' generated : '+clrs['g']+os.path.basename(outpath)+clrs['n']+'\n')


def analyse_hits(hit):
    output = []
    hit_code, hit_chain = hit.split(':')
    if templatedmodel is False:
        output.append('\nHit '+clrs['p']+hit_code.lower()+clrs['n']+', Chain: '+clrs['y']+hit_chain+clrs['n']+', Score: '+clrs['c']+str(hits[hit][0])+clrs['n']+', %id: '+clrs['c']+str(hits[hit][1])+clrs['n']+', Coverage: '+clrs['c']+str(hits[hit][2])+clrs['n'])
    else:
        output.append('\nTemplated Hit '+clrs['p']+hit_code.lower()+clrs['n']+', Chain: '+clrs['y']+hit_chain+clrs['n']+', %id: '+clrs['c']+str(hits[hit])+clrs['n'])
    middle_letters = hit_code[1:3]
    hit_pdb = os.path.join(pdb_homo_archive, middle_letters, hit_code+'.pdb.gz')
    if not os.path.isfile(hit_pdb):
        output.append(clrs['r']+'Hit '+str(hit_code)+' not found in oligomeric database. Disregarding it.'+clrs['n']+'\n')
        return hit_code, None, None, None, '\n'.join(output)
    hit_pdb_name, hit_structure, hit_nchains = pctools.parse_pdb_structure(hit_pdb)
    if not hit_structure:
        output.append(clrs['r']+'Hit '+clrs['y']+hit_pdb+clrs['r']+' could not be parsed.\n'+clrs['n'])
        return hit, None, None, None, None, None, None, '\n'.join(output)
    hit_nchains, hit_seqs, hit_chain_ids = pctools.extract_seqs(hit_structure, 0)
    hit_pids, protein_bool = pctools.get_pairwise_ids(hit_seqs, hit_nchains)
    n = 1
    candidate_chains = set()
    for id in hit_pids:
        if (id[1] == hit_chain or id[2] == hit_chain) and id[0] >= 90:
            if verbosity > 0:
                output.append('Identity between chains '+clrs['y']+str(id[2])+clrs['n']+' and '+clrs['y']+str(id[1])+clrs['n']+' is '+clrs['g']+str(round(id[0], 2))+"%"+clrs['n']+".")
            n += 1
            candidate_chains.add(id[1])
            candidate_chains.add(id[2])
    if all(id[0] > 90 for id in hit_pids) and protein_bool is True:
        hetero_complex = False
    else:
        hetero_complex = True

    initial_homo_chains = n
    if n == 1:
        output.append('No similar chains. Would yield '+clrs['y']+oligo_dict[n]+clrs['n']+' model.\n')
        return hit, None, None, None, hit_nchains, initial_homo_chains, hetero_complex, '\n'.join(output)
    else:
        chain_string = ','.join(candidate_chains)
        output.append('Similar chains '+clrs['y']+chain_string+clrs['n']+' would yield '+clrs['y']+'HOMO-'+oligo_dict[n]+clrs['n']+' model.')
        output.append('\nRunning '+clrs['b']+'PISA'+clrs['n']+' for '+clrs['p']+hit_code+clrs['n']+'...')
        pisa_output, pisa_error = pctools.run_pisa(hit_pdb, hit_chain, verbosity, gen_monomer_data=False, gen_oligomer_data=True)
        output.append(pisa_output)
        if pisa_error is True:
            output.append(clrs['r']+'Disregarding hit '+str(hit_code)+clrs['n']+'\n')
            return hit, None, None, None, hit_nchains, initial_homo_chains, hetero_complex, '\n'.join(output)
        xml_out = os.path.join(workdir, hit_code+hit_chain+'_CHOIR_PisaInterfaces.xml')
        try:
            interfaces_list, interfaces_output = pctools.parse_interfaces(xml_out, candidate_chains, verbosity)
        except et.ParseError:
            output.append(clrs['r']+'Failed to parse xml file for PISA interfaces. Disregarding hit.'+clrs['n'])
            return hit, None, None, None, hit_nchains, initial_homo_chains, hetero_complex, '\n'.join(output)
        if verbosity > 0:
            output.append(interfaces_output)
        cluster_dict, largest_cluster, confirm_homo_output = confirm_homo_state(hit_code, candidate_chains, interfaces_list)
        output.append(confirm_homo_output)
        if candidate_chains == set(cluster_dict[largest_cluster]):
            output.append('Homo-oligomeric state '+clrs['g']+'CONFIRMED'+clrs['n']+' by '+clrs['b']+' PISA'+clrs['n']+'.')
        else:
            output.append('Homo-oligomeric state '+clrs['r']+'NOT CONFIRMED'+clrs['n']+' by '+clrs['b']+' PISA'+clrs['n']+'.')
            if len(cluster_dict[largest_cluster]) == 1:
                output.append('A '+clrs['r']+oligo_dict[len(cluster_dict[largest_cluster])]+clrs['n']+' cluster means no homo-oligomeric interface was detected. Disregarding hit.')
                return hit, None, None, None, hit_nchains, initial_homo_chains, hetero_complex, '\n'.join(output)
        return hit, cluster_dict[largest_cluster], interfaces_list, cluster_dict, hit_nchains, initial_homo_chains, hetero_complex, '\n'.join(output)


def run_tmhmm2(fasta, args):
    tmhmm2_cmd = [tmhmm2_exe, fasta]
    print('\nRunning '+clrs['b']+'TMHMM'+clrs['n']+' to assess likely transmembrane residues...')
    tmhmm2_out = os.path.join(workdir, pdb_name+'_CHOIR_TMHMM2.dat')
    tmhmm2_cmd = [tmhmm2_exe, fasta]
    pctools.printv(clrs['b']+'TMHMM2'+clrs['n']+' command line: '+' '.join(tmhmm2_cmd), args.verbosity)
    with open(tmhmm2_out, 'w') as f:
        subprocess.run(tmhmm2_cmd, stdout=f, check=True)
    print('Done running '+clrs['b']+'TMHMM'+clrs['n']+'. Data file written to '+clrs['g']+os.path.basename(tmhmm2_out)+clrs['n'])
    return tmhmm2_out


def parse_tmhmm_output(tmhmm_out, residue_index_mapping, args):
    tmspans = []
    for line in open(tmhmm_out, 'r'):
        if not line.startswith("#"):
            if line.split()[2] == 'TMhelix':
                tmspans.append((residue_index_mapping[int(line.split()[3])], residue_index_mapping[int(line.split()[4])]))
    print(clrs['b']+'TMHMM'+clrs['n']+' found '+clrs['c']+str(len(tmspans))+clrs['n']+' transmembrane helices')
    tmresidues = []
    n = 0
    for tmspan in tmspans:
        n += 1
        print('TM Helix '+str(n)+': '+str(tmspan[0])+' <-> '+ str(tmspan[1]))
        for i in range(tmspan[0], tmspan[1]):
            tmresidues.append(i)
    return tmspans, tmresidues


def search_homologues(fasta_file, report, args):
    # Search for homologous proteins in all three CHOIR databases
    best_score = {}

    # Search homo-oligomers database
    hits = blast_protomer(fasta_file, homoblast, args.max_candidates, args.psiblast_iterations, args.psiblast_threads, args.psiblast_params, args.verbosity)
    if hits:
        best_score['HOMO-OLIGOMERIC'] = list(hits.items())[0][1][0]
        for hit in hits:
            hit_code, hit_chain = hit.split(':')
            print('Hit '+clrs['p']+hit_code.lower()+clrs['n']+', Chain: '+clrs['y']+hit_chain+clrs['n']+', Score: '+clrs['c']+str(hits[hit][0])+clrs['n']+', %id: '+clrs['c']+str(hits[hit][1])+clrs['n']+', Coverage: '+clrs['c']+str(hits[hit][2])+clrs['n'])
    else:
        print('No hit found in Homo-oligomers database')
        best_score['HOMO-OLIGOMERIC'] = 0

    # Search monomers database
    mono_hits = blast_protomer(fasta_file, monoblast, 3, args.psiblast_iterations, args.psiblast_threads, args.psiblast_params, args.verbosity)
    if mono_hits:
        best_score['MONOMERIC'] = list(mono_hits.items())[0][1][0]
        for hit in mono_hits:
            hit_code, hit_chain = hit.split(':')
            print('Hit '+clrs['p']+hit_code.lower()+clrs['n']+', Chain: '+clrs['y']+hit_chain+clrs['n']+', Score: '+clrs['c']+str(mono_hits[hit][0])+clrs['n']+', %id: '+clrs['c']+str(mono_hits[hit][1])+clrs['n']+', Coverage: '+clrs['c']+str(mono_hits[hit][2])+clrs['n'])
    else:
        print('No hit found in monomers database')
        best_score['MONOMERIC'] = 0

    # Search hetero-oligomers database
    hetero_hits = blast_protomer(fasta_file, heteroblast, 3, args.psiblast_iterations, args.psiblast_threads, args.psiblast_params, args.verbosity)
    if hetero_hits:
        best_score['HETERO-OLIGOMERIC'] = list(hetero_hits.items())[0][1][0]
        for hit in hetero_hits:
            hit_code, hit_chain = hit.split(':')
            print('Hit '+clrs['p']+hit_code.lower()+clrs['n']+', Chain: '+clrs['y']+hit_chain+clrs['n']+', Score: '+clrs['c']+str(hetero_hits[hit][0])+clrs['n']+', %id: '+clrs['c']+str(hetero_hits[hit][1])+clrs['n']+', Coverage: '+clrs['c']+str(hetero_hits[hit][2])+clrs['n'])
    else:
        print('No hit found in Hetero-oligomers database')
        best_score['HETERO-OLIGOMERIC'] = 0

    # Define highest scoring state
    highest_scoring_state = max(best_score, key=best_score.get)
    report['highest_scoring_state'] = highest_scoring_state

    # If hits were found in any of the three databases: Get H3O-Score
    if max(best_score.values()) > 0:
        print('\nThe input protein is likely '+clrs['y']+highest_scoring_state+clrs['n']+'.')
        report['homo_oligomeric_over_other_score'] = round(best_score['HOMO-OLIGOMERIC']/max(best_score.values()), 2)
    else:
        report['homo_oligomeric_over_other_score'] = 'NA'
        print('PSI-BLAST found hits in none of the three local databases')
        if not args.allow_monomers:
            return None, report

    # Report and proceed
    if highest_scoring_state != 'HOMO-OLIGOMERIC' and args.allow_monomers:
        print(clrs['y']+'\nMonomer/Protomer building option selected. Proceeding accordingly.'+clrs['n'])
        seqres_hits = blast_protomer(fasta_file, seqres, args.max_candidates, args.psiblast_iterations, args.psiblast_threads, args.psiblast_params, args.verbosity)
        if seqres_hits:
            for hit in seqres_hits:
                hit_code, hit_chain = hit.split(':')
                print('Hit '+clrs['p']+hit_code.lower()+clrs['n']+', Chain: '+clrs['y']+hit_chain+clrs['n']+', Score: '+clrs['c']+str(seqres_hits[hit][0])+clrs['n']+', %id: '+clrs['c']+str(seqres_hits[hit][1])+clrs['n']+', Coverage: '+clrs['c']+str(seqres_hits[hit][2])+clrs['n'])
            print(clrs['y']+'\nCandidate templates found for monomer/protomer building.'+clrs['n'])
            if args.symmetry is True:
                print(clrs['y']+'Disabling option for symmetry constraints.'+clrs['n'])
                args.symmetry = False
            return seqres_hits, report
        else:
            print('PSI-BLAST found NO hits in seqres database')
            return None, report
    else:
        args.allow_monomers = False

    if not hits:
        print('PSI-BLAST found NO hits in Homo-Oligomeric database')
        print('The H3O-Score is '+clrs['c']+str(report['homo_oligomeric_over_other_score'])+clrs['n']+'.')
    elif hits and (mono_hits or hetero_hits):
        if highest_scoring_state == 'HOMO-OLIGOMERIC':
            print('\nHighest scoring hit forms Homo-Oligomeric interfaces.')
            print('The H3O-Score is '+clrs['c']+str(report['homo_oligomeric_over_other_score'])+clrs['n']+'.\n')
        else:
            print('\nHighest scoring hit does not form Homo-Oligomeric interfaces.'+clrs['n']+'.')
            tolerated_score = (100-args.tolerance) * max(best_score.values())/100
            if best_score['HOMO-OLIGOMERIC'] >= tolerated_score:
                print('Nevertheless, ProtCHOIR is still going to attempt building models, since the best Homo-Oligomeric score is within the tolerance range ('+clrs['c']+str(args.tolerance)+'%'+clrs['n']+' of the best hit\'s score).')
                print('The H3O-Score is '+clrs['c']+str(report['homo_oligomeric_over_other_score'])+clrs['n']+'.\n')
                for hit in hits:
                    hit_code, hit_chain = hit.split(':')
                    print('\nHit '+clrs['p']+hit_code.lower()+clrs['n']+', Chain: '+clrs['y']+hit_chain+clrs['n']+', Score: '+clrs['c']+str(hits[hit][0])+clrs['n']+', %id: '+clrs['c']+str(hits[hit][1])+clrs['n']+', Coverage: '+clrs['c']+str(hits[hit][2])+clrs['n'])

            else:
                print('ProtCHOIR is NOT going to attempt building models, since the best Homo-Oligomeric score is  NOT within the tolerance range ('+clrs['c']+str(args.tolerance)+'%'+clrs['n']+' of the best hit\'s score).')
                print('The H3O-Score is '+clrs['c']+str(report['homo_oligomeric_over_other_score'])+clrs['n']+'.\n')
                return None, report
    else:
        print('Highest scoring hit forms Homo-Oligomeric interfaces.\n')
        print('The H3O-Score is '+clrs['c']+str(report['homo_oligomeric_over_other_score'])+clrs['n']+'.')
        for hit in hits:
            hit_code, hit_chain = hit.split(':')
            print('\nHit '+clrs['p']+hit_code.lower()+clrs['n']+', Chain: '+clrs['y']+hit_chain+clrs['n']+', Score: '+clrs['c']+str(hits[hit][0])+clrs['n']+', %id: '+clrs['c']+str(hits[hit][1])+clrs['n']+', Coverage: '+clrs['c']+str(hits[hit][2])+clrs['n'])

    return hits, report


# Main Function
###############################################################################
def analyze_protomer(input_file, report, args):
    global workdir
    global pdb_name
    global templatedmodel
    global hits
    global verbosity
    verbosity = args.verbosity
    workdir = os.getcwd()
    ignoretemplated = args.ignoretemplated
    max_candidates = args.max_candidates
    largest_oligo_complexes = collections.OrderedDict()
    report['sequence_mode'] = str(args.sequence_mode)

    # Default Mode ############################################################################################################################
    if args.sequence_mode is False:
        pctools.print_section(1, 'PROTOMER ANALYSIS')
        # Subsection 1[a] #######################################################################
        pctools.print_subsection('1[a]', 'Loading input files')
        filename = os.path.basename(input_file)
        print('Will now begin analysis of '+clrs['p']+filename+clrs['n']+'\n')
        print('Loading structure...')
        pattern = 'REMARK   6 TEMPLATE:'
        templates_ids = {}

        # Check if we are dealing with a templated model!
        templatedmodel = False
        if ignoretemplated is False:
            for line in open(input_file, 'r'):
                if re.search(pattern, line):
                    print(clrs['y']+'Dealing with a templated model! Kudos!'+clrs['n'])
                    templatedmodel = True
                    break
        report['templatedmodel'] = str(templatedmodel)
        # Check number of chains in input and clean input file
        pdb_name, structure, nchains = pctools.parse_any_structure(input_file)
        clean_input_file = os.path.join(workdir, pdb_name+'_Clean.pdb')
        io.set_structure(structure)
        io.save(clean_input_file, pctools.SelectIfCA())
        # Reload clean file
        pdb_name, structure, nchains = pctools.parse_any_structure(clean_input_file)
        if nchains == 1:
            print('Structure '+clrs['p']+pdb_name+clrs['n']+' is '+clrs['y']+'MONOMERIC'+clrs['n']+' as expected')
        else:
            print('Structure '+clrs['p']+pdb_name+clrs['n']+' contains '+clrs['y']+str(nchains)+clrs['n']+' chains.')
            print('Will consider only first chain')
            clean_input_file = pctools.split_chains(pdb_name, structure, workdir)
            pdb_name, structure, nchains = pctools.parse_any_structure(clean_input_file)

        # Extract sequence of (first) chain in structure
        nchains, seqs, chain_ids = pctools.extract_seqs(structure, 0)
        sequence = seqs[0][1]
        report['protomer_residues'] = str(len(sequence))
        fasta_file = write_fasta(sequence)

        # Subsection 1[b] #######################################################################
        pctools.print_subsection('1[b]', 'Oligomeric homologues search')
        # If not a Templated model, search for homologous proteins in all three CHOIR databases
        if templatedmodel is False:
            hits, report = search_homologues(fasta_file, report, args)
            if not hits:
                report['exit'] = '2'
                return None, report, args
        else:
            hits = parse_templated_model(sequence, input_file)
            report['highest_scoring_state'] = 'NA'
            report['homo_oligomeric_over_other_score'] = 'NA'

            if not hits:
                print('No Templated-determined hits were found in the homo-oligomeric database. Try using --ignore-templated.\n')
                report['exit'] = '1'
                return None, report, args

        # Subsection 1[c] #######################################################################
        pctools.print_subsection('1[c]', 'Protomer structure check')
        # Run PISA for monomer and get surface residues
        print('\nRunning '+clrs['b']+'PISA'+clrs['n']+' for '+clrs['p']+pdb_name+clrs['n']+'...')
        output, pisa_error, monomer_data = pctools.run_pisa(clean_input_file, '', args.verbosity, gen_monomer_data=True, gen_oligomer_data=False)
        if output:
            print(output)
        protomer_surface_residues = pctools.get_areas(monomer_data)
        residue_index_mapping = pctools.map_residue_index(protomer_surface_residues)
        print('Done running '+clrs['b']+'PISA'+clrs['n']+'.')


        # Get likely transmembrane residues
        tmhmm_out = run_tmhmm2(fasta_file, args)
        tmdata = parse_tmhmm_output(tmhmm_out, residue_index_mapping, args)
        report['tmspans'] = str(len(tmdata[0]))

        # Run Molprobity for monomer
        protomer_molprobity, molprobity_output = pctools.run_molprobity(clean_input_file, args)
        print(molprobity_output)
        print(clrs['y']+'MOLPROBITY ASSESSMENT'+clrs['n'])
        print('Rama. Fav.\t'+str(protomer_molprobity['rama_fav']))
        print('Rama. Out.\t'+str(protomer_molprobity['rama_out']))
        print('Rot. Out.\t'+str(protomer_molprobity['rot_out']))
        print('CBeta Dev.\t'+str(protomer_molprobity['cb_dev']))
        print('Clashscore\t'+str(protomer_molprobity['clashscore']))
        print('Molprob. Score\t'+str(protomer_molprobity['molprobity_score']))
        report['protomer_molprobity'] = protomer_molprobity['molprobity_score']
        report['protomer_clashscore'] = protomer_molprobity['clashscore']

        # Subsection 1[d] #######################################################################
        if not args.skip_conservation:
            pctools.print_subsection('1[d]', 'Sequence conservation analysis')
            # Use PSI-BLAST to search UniRef50 and return hits
            uni50hits = blast_protomer(fasta_file, uniref50, 50, 1, args.psiblast_threads, args.psiblast_params, args.verbosity)
            if not uni50hits:
                print('PSI-BLAST found NO hits in Uniref50 database. Skipping conservation analysis.')
                args.skip_conservation = True
                entropies = None
                z_entropies = None
                minx = None
                maxx = None
                report['protomer_plot'], report['protomer_exposed_area'], report['protomer_hydrophobic_area'], report['protomer_conserved_area'], minx, maxx, analysis_output = pctools.plot_analysis(pdb_name, protomer_surface_residues, None, None, tmdata, args)

            if uni50hits:
                multi_fasta = generate_msa_input(uni50hits, fasta_file, args.verbosity)
                msa_file = run_mafft(multi_fasta, args)
                trimmed_msa = trim_msa(msa_file)
                msa_dict_trim = parse_msa(trimmed_msa)
                entropies = relative_entropy(msa_dict_trim, protomer_surface_residues)
                z_entropies = calc_z_scores(entropies)
                report['protomer_plot'], report['protomer_exposed_area'], report['protomer_hydrophobic_area'], report['protomer_conserved_area'], minx, maxx, analysis_output = pctools.plot_analysis(pdb_name, protomer_surface_residues, entropies, z_entropies, tmdata, args)
                print(analysis_output)
                monomer_conservation = map_conservation(structure, z_entropies)
                report['protomer_figure'] = pctools.pymol_screenshot_mono(monomer_conservation, z_entropies, args)
        else:
            print(clrs['y']+"Skipping section 1[d] - Sequence conservation analysis"+clrs['n']+"\n")
            report['protomer_plot'], report['protomer_exposed_area'], report['protomer_hydrophobic_area'], report['protomer_conserved_area'], minx, maxx, analysis_output = pctools.plot_analysis(pdb_name, protomer_surface_residues, None, None, tmdata, args)


    # Sequence Mode ###########################################################################################################################
    elif args.sequence_mode is True:
        pctools.print_section(1, 'PROTOMER ANALYSIS - SEQUENCE MODE')
        # Subsection 1[a] #######################################################################
        pctools.print_subsection('1[a]', 'Loading input files')
        if input_file.lower().endswith('.fasta'):
            fasta_file = clean_input_file = input_file
            pdb_name = os.path.basename(fasta_file).split("_CHOIR_MonomerSequence.fasta")[0].replace('.', '_')
            records = list(SeqIO.parse(fasta_file, "fasta"))
            sequence = records[0].seq
            report['protomer_residues'] = str(len(sequence))

        elif args.input_file.endswith('.pdb'):
            report['input_filename'] = os.path.basename(input_file)
            # Check number of chains in input and clean input file
            pdb_name, structure, nchains = pctools.parse_any_structure(input_file)
            clean_input_file = os.path.join(workdir, pdb_name+'_Clean.pdb')
            io.set_structure(structure)
            io.save(clean_input_file, pctools.SelectIfCA())
            # Reload clean file
            pdb_name, structure, nchains = pctools.parse_any_structure(clean_input_file)
            print(pdb_name)
            nchains, seqs, chain_ids = pctools.extract_seqs(structure, 0)
            sequence = seqs[0][1]
            report['protomer_residues'] = str(len(sequence))
            fasta_file = write_fasta(sequence)
            if nchains == 1:
                print('Structure '+clrs['p']+pdb_name+clrs['n']+' is '+clrs['y']+'MONOMERIC'+clrs['n']+' as expected')
            else:
                print('Structure '+clrs['p']+pdb_name+clrs['n']+' contains '+clrs['y']+str(nchains)+clrs['n']+' chains.')
                print('Will consider only first chain')
        templatedmodel = False
        report['templatedmodel'] = str(templatedmodel)

        residue_index_mapping = {}
        for i, aa in enumerate(sequence):
            residue_index_mapping[i] = i
        # Get likely transmembrane residues
        tmhmm_out = run_tmhmm2(fasta_file, args)
        tmdata = parse_tmhmm_output(tmhmm_out, residue_index_mapping, args)
        report['tmspans'] = str(len(tmdata[0]))

        # Subsection 1[b] #######################################################################
        pctools.print_subsection('1[b]', 'Oligomeric homologues search')
        # Search for homologous proteins in all three CHOIR databases
        hits, report = search_homologues(fasta_file, report, args)
        if not hits:
            report['exit'] = '2'
            return None, report, args

        # Subsection 1[c] #######################################################################
        # Inform there will be no structural analysis
        print(clrs['y']+"\nSkipping section 1[c] - Protomer structure check"+clrs['n']+"\n")

        # Subsection 1[d] #######################################################################
        if not args.skip_conservation:
            pctools.print_subsection('1[d]', 'Sequence conservation analysis')
            # Use PSI-BLAST to search UniRef50 and return hits
            uni50hits = blast_protomer(fasta_file, uniref50, 50, 1, args.psiblast_threads, args.psiblast_params, args.verbosity)
            if not uni50hits:
                print('PSI-BLAST found NO hits in Uniref50 database. Skipping conservation analysis')
                args.skip_conservation = True
                entropies = None
                z_entropies = None
                minx = None
                maxx = None
                report['protomer_plot'] = pctools.plot_entropy_only(pdb_name, None, None, tmdata, args)
            if uni50hits:
                multi_fasta = generate_msa_input(uni50hits, fasta_file, args.verbosity)
                msa_file = run_mafft(multi_fasta, args)
                trimmed_msa = trim_msa(msa_file)
                msa_dict_trim = parse_msa(trimmed_msa)
                entropies = relative_entropy(msa_dict_trim, None)
                z_entropies = calc_z_scores(entropies)
                report['protomer_plot'] = pctools.plot_entropy_only(pdb_name, entropies, z_entropies, tmdata, args)
        else:
            print(clrs['y']+"Skipping section 1[d] - Sequence conservation analysis"+clrs['n']+"\n")
            report['protomer_plot'] = pctools.plot_entropy_only(pdb_name, None, None, tmdata, args)



    # Analyse hits ############################################################################################################################
    report['hits'] = {}
    if hits and not args.allow_monomers:
        # Subsection 1[e] #######################################################################
        pctools.print_subsection('1[e]', 'Oligomeric homologues analysis')
        if args.multiprocess is True:
            p = Pool()
            interfaces_dict = {}
            for hitchain, largest_cluster, interfaces_list, cluster_dict, total_chains, initial_homo_chains, hetero_complex, output in p.map_async(analyse_hits, hits).get():
                report['hits'][hitchain] = {}
                report['hits'][hitchain]['hetero_complex'] = hetero_complex
                report['hits'][hitchain]['total_chains'] = total_chains
                report['hits'][hitchain]['initial_homo_chains'] = initial_homo_chains
                print(output)

                if largest_cluster is None:
                    report['hits'][hitchain]['qscore'] = 'NA'
                    report['hits'][hitchain]['final_homo_chains'] = '0'
                    print('-------------------------------------------------------------------')
                    continue
                largest_oligo_complexes[hitchain] = largest_cluster
                report['hits'][hitchain]['final_homo_chains'] = str(len(largest_cluster))
                if args.plot_topologies is True:
                    plot_topology(hitchain, interfaces_list, cluster_dict)
                interfaces_dict[hitchain] = interfaces_list
                print('-------------------------------------------------------------------')
            p.close()
            p.join()
        else:
            interfaces_dict = {}
            for hit in hits:
                hitchain, largest_cluster, interfaces_list, cluster_dict, total_chains, initial_homo_chains, hetero_complex, output = analyse_hits(hit)
                report['hits'][hitchain] = {}
                report['hits'][hitchain]['hetero_complex'] = hetero_complex
                report['hits'][hitchain]['total_chains'] = total_chains
                report['hits'][hitchain]['initial_homo_chains'] = initial_homo_chains
                print(output)

                if largest_cluster is None:
                    report['hits'][hitchain]['qscore'] = 'NA'
                    report['hits'][hitchain]['final_homo_chains'] = '0'
                    print('-------------------------------------------------------------------')
                    continue
                largest_oligo_complexes[hitchain] = largest_cluster
                report['hits'][hitchain]['final_homo_chains'] = str(len(largest_cluster))
                if args.plot_topologies is True:
                    plot_topology(hitchain, interfaces_list, cluster_dict)
                interfaces_dict[hitchain] = interfaces_list
                print('-------------------------------------------------------------------')

        if not largest_oligo_complexes:
            print('**ProtCHOIR'+clrs['r']+' failed '+clrs['n']+'to select good oligomeric templates.\n')
            report['exit'] = '3'
            return None, report, args

        pickle.dump(largest_oligo_complexes, open('CHOIR_OligoComplexes.pickle', 'wb'))

    elif args.allow_monomers:
        print(clrs['y']+"Skipping section 1[e] - Oligomeric homologues analysis"+clrs['n']+"\n")
        # Read Chain correspondences
        stats_dir = os.path.join(pdb_homo_archive, 'stats')
        for hitchain in hits:
            hit_code, hit_chain = hitchain.split(':')
            report['hits'][hitchain] = {}
            report['hits'][hitchain]['hetero_complex'] = 'NA'
            report['hits'][hitchain]['total_chains'] = 'NA'
            report['hits'][hitchain]['initial_homo_chains'] = 'NA'
            report['hits'][hitchain]['qscore'] = 'NA'
            report['hits'][hitchain]['final_homo_chains'] = '1'
            interfaces_dict = None
            largest_oligo_complexes[hitchain] = hit_chain

    for hitchain in hits:
        hit_code, chain = hitchain.split(':')
        report['hits'][hitchain]['hit_code'] = hit_code
        report['hits'][hitchain]['chain'] = chain
        if templatedmodel is False:
            report['hits'][hitchain]['score'] = str(hits[hitchain][0])
            report['hits'][hitchain]['id'] = str(hits[hitchain][1])
            report['hits'][hitchain]['coverage'] = str(hits[hitchain][2])
        else:
            report['hits'][hitchain]['score'] = 'NA'
            report['hits'][hitchain]['id'] = str(round(hits[hitchain], 1))
            report['hits'][hitchain]['coverage'] = 'NA'



    if args.sequence_mode:
        if args.skip_conservation:
            return (pdb_name, clean_input_file, largest_oligo_complexes, interfaces_dict, tmdata), report, args
        elif not args.skip_conservation:
            return (pdb_name, clean_input_file, largest_oligo_complexes, interfaces_dict, entropies, z_entropies, tmdata), report, args

    elif not args.sequence_mode:
        if args.skip_conservation:
            return (pdb_name, clean_input_file, largest_oligo_complexes, interfaces_dict, residue_index_mapping, tmdata), report, args
        elif not args.skip_conservation:
            return (pdb_name, clean_input_file, largest_oligo_complexes, interfaces_dict, entropies, z_entropies, residue_index_mapping, minx, maxx, tmdata), report, args
