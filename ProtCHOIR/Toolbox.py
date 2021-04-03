# Imports
###############################################################################
import os
import re
import gzip
import shutil
import jinja2
import string
import secrets
import parasail
import subprocess
import matplotlib
import collections
import numpy as np
import pandas as pd
matplotlib.use('Agg')
import textwrap as tw
import Bio.PDB as bpp
import itertools as it
import matplotlib.pyplot as plt
import xml.etree.ElementTree as et
from datetime import datetime
from matplotlib.lines import Line2D
import Bio.PDB.Polypeptide as bpp_poly
from ProtCHOIR.Initialise import *



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
'''
Sets of classes and functions that are used by the ProtCHOIR Oligomerization
Pipeline. They are imported as:

import ProtCHOIR.Toolbox as pctools
'''


# Classes
###############################################################################
class FileFormatError(Exception):
    """Exception for invalid file format.
    Thrown when an unexpected value type is read from the file stream
    or when the end of file is reached unexpectedly."""
    pass

class SelectAA(bpp.Select):
    '''
    Biopython select class to select only aminoacids
    Called by: clean_pdb()
    '''
    def accept_residue(self, residue):
        if bpp_poly.is_aa(residue.get_resname(), standard=True):
            return 1
        else:
            return 0

class SelectIfCA(bpp.Select):
    '''
    Biopython select class to select only aminoacids that contain the alfa-carbon
    Called by: clean_pdb()
    '''
    def accept_residue(self, residue):
        if residue.has_id('CA') and bpp_poly.is_aa(residue.get_resname(), standard=True) and residue.id[0] == " ":
            return 1
        else:
            return 0


class SelectChain(bpp.Select):
    '''
    Biopython select class to select a single chain id.
    Called by: single_chain()
    '''
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        if chain.id == self.chain_id:
            return 1
        else:
            return 0


class SelectChains(bpp.Select):
    '''
    Biopython select class to select a single chain id.
    Called by: single_chain()
    '''
    def __init__(self, chain_ids):
        self.chain_ids = chain_ids

    def accept_chain(self, chain):
        if chain.id in self.chain_ids:
            return 1
        else:
            return 0


class SelectResidues(bpp.Select):
    '''
    Biopython select class to select all residues in a given list.
    Called by: single_chain()
    '''
    def __init__(self, reslist):
        self.reslist = reslist

    def accept_residue(self, residue):
        if residue in self.reslist:
            return 1
        else:
            return 0


class SelectAtoms(bpp.Select):
    '''
    Biopython select class to select all atoms in a given list.
    Called by: single_chain()
    '''
    def __init__(self, atomlist):
        self.atomlist = atomlist

    def accept_atom(self, atom):
        if atom in self.atomlist:
            return 1
        else:
            return 0


# Functions
###############################################################################
def printv(text, verbosity):
    '''
    Function to print if verbosity is invoked.
    '''
    if verbosity == 1:
        print(text)


def print_section(n, name):
    '''
    Function to print section header.
    '''
    print('\n'+clrs['y']+'+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n' +
          clrs['r']+'|' +
          clrs['y']+'S'+clrs['r']+'|' +
          clrs['y']+'E'+clrs['r']+'|' +
          clrs['y']+'C'+clrs['r']+'|' +
          clrs['y']+'T'+clrs['r']+'|' +
          clrs['y']+'I'+clrs['r']+'|' +
          clrs['y']+'O'+clrs['r']+'|' +
          clrs['y']+'N'+clrs['r']+'|' +
          clrs['c']+' '+str(n)+clrs['r']+' | ' +
          clrs['c']+str(name)+'\n' +
          clrs['y']+'+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n' +
          clrs['n'])
    if name == "Runtime Arguments":
        with open('CHOIR_Progress.out', 'w') as f:
            f.write(datetime.now().strftime("%H:%M:%S")+": Starting new ProtCHOIR run\n")
    else:
        with open('CHOIR_Progress.out', 'a') as f:
            f.write("\n"+datetime.now().strftime("%H:%M:%S")+": "+str(name)+"...\n")


def section(n, name):
    '''
    Function to return section header.
    '''
    section_string = str('\n'+clrs['y']+'+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n' +
                         clrs['r']+'|' +
                         clrs['y']+'S'+clrs['r']+'|' +
                         clrs['y']+'E'+clrs['r']+'|' +
                         clrs['y']+'C'+clrs['r']+'|' +
                         clrs['y']+'T'+clrs['r']+'|' +
                         clrs['y']+'I'+clrs['r']+'|' +
                         clrs['y']+'O'+clrs['r']+'|' +
                         clrs['y']+'N'+clrs['r']+'|' +
                         clrs['c']+' '+str(n)+clrs['r']+' | ' +
                         clrs['c']+str(name)+'\n' +
                         clrs['y']+'+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n' +
                         clrs['n'])

    return section_string


def print_subsection(n, name):
    '''
    Function to print section header.
    '''
    print('\n'+clrs['r']+'|' +
          clrs['y']+'S' +
          clrs['y']+'u' +
          clrs['y']+'b' +
          clrs['y']+'s' +
          clrs['y']+'e' +
          clrs['y']+'c' +
          clrs['y']+'t' +
          clrs['y']+'i' +
          clrs['y']+'o' +
          clrs['y']+'n' +
          clrs['c']+' '+str(n)+clrs['r']+'| ' +
          clrs['c']+str(name)+'\n' +
          clrs['y']+'+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n' +
          clrs['n'])
    with open('CHOIR_Progress.out', 'a') as f:
        f.write(datetime.now().strftime("%H:%M:%S")+": "+str(name)+"...\n")


def subsection(n, name):
    '''
    Function to return section header.
    '''
    subsection_string = str('\n'+clrs['r']+'|' +
                            clrs['y']+'S' +
                            clrs['y']+'u' +
                            clrs['y']+'b' +
                            clrs['y']+'s' +
                            clrs['y']+'e' +
                            clrs['y']+'c' +
                            clrs['y']+'t' +
                            clrs['y']+'i' +
                            clrs['y']+'o' +
                            clrs['y']+'n' +
                            clrs['c']+' '+str(n)+clrs['r']+'| ' +
                            clrs['c']+str(name)+'\n' +
                            clrs['y']+'+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n' +
                            clrs['n'])

    return subsection_string

def print_sorry():
    print('**We are '+clrs['y']+'t'+clrs['r']+'e'+
          clrs['y']+'r'+clrs['r']+'r'+clrs['y']+'i'+
          clrs['r']+'b'+clrs['y']+'l'+clrs['r']+'y'+
          clrs['n']+' sorry... =(\n')


def three_to_one(i):
    '''
    Function that returns the one letter code of a given amino acid, based on
    the aa3to1 dictionary in GlobalVars.py
    Called by: extract_seqs()
    '''
    return aa3to1[i]


def gzip_pdb(file):
    '''
    Runs gunzip to compress file. DELETES ORIGINAL.
    Called by: UpdateDatabases.py:curate_homoDB()
    '''
    with open(file, 'rb') as infile, gzip.open(file+".gz", 'wb') as outfile:
        outfile.writelines(infile)
    os.remove(file)


def parse_pdb_contents(pdb):
    '''
    Returns the pdb code and the full text content.
    Called by: UpdateDatabases.py:curate_homoDB()
    '''
    if pdb.endswith(".ent") or pdb.endswith(".pdb") or pdb.endswith(".ent.gz") or pdb.endswith(".pdb1") or pdb.endswith(".pdb1.gz") or pdb.endswith(".pdb.gz"):
        pdb_name = pdb.split('.')[0].split("/")[-1][-4:]
        if pdb.endswith(".gz"):
            contents = gzip.open(pdb, 'rt').readlines()
        else:
            contents = open(pdb, 'rt').readlines()
    return pdb_name, contents


def parse_pdb_structure(pdb):
    '''
    Returns the pdb code and the biopython structure object and the number
    of chains in a pdb file.
    Called by: UpdateDatabases.py:curate_homoDB()
    '''
    if pdb.endswith(".ent") or pdb.endswith(".pdb") or pdb.endswith(".ent.gz") or pdb.endswith(".pdb1") or pdb.endswith(".pdb1.gz") or pdb.endswith(".pdb.gz"):
        pdb_name = pdb.split('.')[0].split("/")[-1][-4:]
        if pdb.endswith(".gz"):
            pdb_file = gzip.open(pdb, 'rt')
        else:
            pdb_file = open(pdb)
        try:
            structure = p.get_structure(pdb_name, pdb_file)
            nchains = len(bpp.Selection.unfold_entities(structure, 'C'))
        except:
            print("Structure "+pdb_name+" could not be strictly parsed.")
            structure = None
            nchains = None
    return pdb_name, structure, nchains


def parse_any_structure(pdb):
    '''
    Returns the pdb code and the biopython structure object and the number
    of chains in a pdb file.
    Called by: UpdateDatabases.py:curate_homoDB()
    '''
    if pdb.endswith(".ent") or pdb.endswith(".pdb") or pdb.endswith(".ent.gz") or pdb.endswith(".pdb1") or pdb.endswith(".pdb1.gz") or pdb.endswith(".pdb.gz"):
        pdb_name = os.path.basename(pdb).split(".pdb")[0].replace('.', '_')
        if pdb.endswith(".gz"):
            pdb_file = gzip.open(pdb, 'rt')
        else:
            pdb_file = open(pdb)
        try:
            structure = p.get_structure(pdb_name, pdb_file)
            nchains = len(bpp.Selection.unfold_entities(structure, 'C'))
        except:
            print("Structure "+pdb_name+" could not be strictly parsed.")
            structure = None
            nchains = None
    else:
        print('\n'+clrs['r']+'The provided file is not in a supported structure format.'+clrs['n'])
        structure = None
        nchains = None
    return pdb_name, structure, nchains


def split_states(structure):
    chains = bpp.Selection.unfold_entities(structure, 'C')
    str_id = structure.id
    new_structure = bpp.Structure.Structure(str_id)
    new_structure.add(bpp.Model.Model(0))
    chain_correspondence_dict = {}
    n = 1
    for chain in chains:
        original = str(chain.id)
        new = numalpha[str(n)]
        chain.id = 'X'+new
        chain_correspondence_dict[new] = original
        n += 1
        if n > 60:
            break
    n = 1
    for chain in chains:
        chain.id = str(chain.id)[1]
        new_structure[0].add(chain)
        n += 1
        if n > 60:
            break

    return new_structure, chain_correspondence_dict


def split_chains(pdb_name, structure, outdir):
    chains = structure.get_chains()
    single_chain_files = []
    for chain in chains:
        chain_id = chain.id
        single_chain_file = os.path.join(outdir, pdb_name+'_'+chain_id+".pdb")
        io.set_structure(structure)
        io.save(single_chain_file, SelectChain(chain_id))
        single_chain_files.append(single_chain_file)
    return single_chain_files[0]


def count_chains(structure):
    '''
    Uses Biopython to obtain the number of chains.
    Called by: NONE
    '''
    nchains = 0
    for model in structure:
        if model.id == 0:
            for chain in model:
                nres = 0
                nchains += 1
                for residue in chain:
                    nres += 1
    return nchains, nres


def extract_seqs(structure, defmodel):
    '''
    Uses Biopython to count the numer of chains and to extract the
    each chain's sequence as a list of sequences.
    Called by: UpdateDatabases.py:curate_homoDB()
    '''
    nchains = 0
    nprotchains = 0
    for model in structure:
        if model.id == defmodel:
            seqs = []
            chain_ids = []
            for chain in model:
                nchains += 1
                seqlist = []
                for residue in chain:
                    if bpp_poly.is_aa(residue.get_resname(), standard=False):
                        try:
                            seqlist.append(three_to_one(residue.get_resname()))
                        except:
                            seqlist.append('X')
                seq = str("".join(seqlist))
                if seq != '':
                    nprotchains += 1
                seqs.append([chain.id, seq])
                chain_ids.append(chain.id)
    return nprotchains, seqs, chain_ids


def get_pairwise_ids(seqs, nchains):
    '''
    Receives a list of sequences and calculates the identity matrix
    among them, as a list, using Biopython and itertools.
    Called by:  UpdateDatabases.py | curate_homoDB()
                AnalyzeProtomer.py | analyze_protomer()
    '''
    ids = []
    are_proteins_list = []
    for a, b in it.combinations(range(nchains), 2):
        are_proteins = False
        if seqs[b][1] and seqs[a][1]:
            if set(seqs[b][1]) != {'X'} and seqs[a][1] != {'X'}:
                are_proteins = True
                alignment = parasail.sg_stats_striped_16(seqs[b][1], seqs[a][1], 10, 1, parasail.blosum62)
                # Discard small alignments (smaller than 10% the length of longest sequence)
                if alignment.length < 0.1*max([len(seqs[b][1]), len(seqs[a][1])]):
                    percent_id = 0
                else:
                    percent_id = (alignment.matches)/alignment.length*100
            else:
                percent_id = 0
        else:
            percent_id = 0
        ids.append([percent_id, seqs[b][0], seqs[a][0]])
        are_proteins_list.append(are_proteins)
    if all(are_proteins_list) is False:
        proteinpair = False
    else:
        proteinpair = True
    return ids, proteinpair

def find_most_similar_chain(sequence, pdbcode):
    middle_letters = pdbcode[1:3]
    reference_file = os.path.join(pdb_homo_archive, middle_letters, pdbcode+'.pdb.gz')
    if not os.path.isfile(reference_file):
        return None, None
    else:
        structure = parse_any_structure(reference_file)[1]
        seqs = extract_seqs(structure, 0)[1]
        chain_identities = {}
        for ref_seq in seqs:
            if ref_seq[1]:
                if set(ref_seq[1]) != {'X'}:
                    alignment = parasail.sg_stats_striped_16(sequence, ref_seq[1], 10, 1, parasail.blosum62)
                    # Discard small alignments (smaller than 10% the length of longest sequence)
                    if alignment.length < 0.1*max([len(sequence), len(ref_seq[1])]):
                        chain_identities[ref_seq[0]] = 0
                    else:
                        percent_id = (alignment.matches)/alignment.length*100
                        chain_identities[ref_seq[0]] = percent_id
                else:
                    chain_identities[ref_seq[0]] = 0
            else:
                chain_identities[ref_seq[0]] = 0
        hit_chain = max(chain_identities.keys(), key=(lambda x: chain_identities[x]))
        pid = chain_identities[hit_chain]
    return hit_chain, pid


def is_valid_sequence(seq):
    if set(seq) != {'X'} and seq != '' and len(set(seq)) > 4:
        return True
    else:
        return False


def get_homo_chains(hit_sequence, pdb):
    '''
    Receives a sequence and a PDB file  and calculates the identity of
    given sequence and all pdb chains using Biopython Parasail.
    Returns a list of homo oligomeric chains in the PDB file.
    Called by: AnalyzeProtomer.py | parse_interfaces()
    '''
    ids = []
    are_proteins_list = []
    homochains = []
    pdb_name, structure, nchains = parse_any_structure(pdb)
    nchains, seqs, chain_ids = extract_seqs(structure, 0)
    print(hit_sequence)
    for seq in seqs:
        are_proteins = False
        print(are_proteins)
        if seq[1]:
            print(seq[1])
            if set(seq) != {'X'}:
                are_proteins = True
                print(are_proteins)
                alignment = parasail.sg_stats_striped_16(hit_sequence, seq[1], 10, 1, parasail.blosum62)
                if alignment.length == 0:
                    percent_id = 0
                else:
                    percent_id = (alignment.matches)/alignment.length*100
            else:
                percent_id = 0
        else:
            percent_id = 0
        print(percent_id)
        ids.append([percent_id, 'hit', seq[0]])
        if percent_id > 90 and are_proteins is True:
            homochains.append(seq[0])
            print(homochains)
    homochains = sorted(homochains)
    return homochains


def check_interfaces(structure, ch1, ch2):
    for chain in bpp.Selection.unfold_entities(structure, 'C'):
        if chain.id == ch1:
            ch1 = chain
    for chain in bpp.Selection.unfold_entities(structure, 'C'):
        if chain.id == ch2:
            ch2 = chain
    atoms_ch1 = [atom for atom in bpp.Selection.unfold_entities(ch1, 'A') if bpp_poly.is_aa(atom.get_parent(), standard=True)]
    atoms_ch2 = [atom for atom in bpp.Selection.unfold_entities(ch2, 'A') if bpp_poly.is_aa(atom.get_parent(), standard=True)]
    contacts = []
    for atom in atoms_ch1:
        ns = bpp.NeighborSearch(atoms_ch2)
        center = atom.get_coord()
        neighbors = ns.search(center, 3.5)
        if neighbors:
            contact = []
            contact1 = []
            contact2 = []
            contact1.append(atom)
            contact.append(contact1)
            for neighbor in neighbors:
                contact2.append(neighbor)
            contact.append(contact2)
            contacts.append(contact)
    return contacts


def get_annotated_states(contents):
    author_remark = r"REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT:"
    software_remark = r"REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE:"
    author = []
    software = []
    for line in contents:
        if re.match(author_remark, line):
            try:
                author.append(str(oligo_dict_inv[line.split(": ")[-1].split(" ")[0]]))
            except:
                author.append('?')
        if re.match(software_remark, line):
            try:
                software.append(str(oligo_dict_inv[line.split(": ")[-1].split(" ")[0]]))
            except:
                software.append('?')
    if not author:
        author.append('NA')
    if not software:
        software.append('NA')
    return author, software

def is_nmr(contents):
    exp_remark = r"REMARK 210  EXPERIMENT TYPE                : NMR"
    for line in contents:
        if re.match(exp_remark, line):
            return True


def author_agrees(oligo_dict, contents, nchains):
    '''
    Searches in the original PDB file for the oligomeric
    status determined by the author.
    Called by: clean_and_sort()
    '''
    pattern = r"REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: "+oligo_dict[nchains]
    if re.search(pattern, contents):
        return True
    else:
        return False


def ls(dir):
    obj = subprocess.Popen(
        ['ls', '-ltrh', dir],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        close_fds=True,
        shell=False,
        cwd=None,
        env=None)
    result = obj.communicate(None)[0].decode()
    obj.stdin.close()
    print(result)


def run_pisa(pdb, chain, verbosity, gen_monomer_data=False, gen_oligomer_data=False):

    alphabet = string.ascii_letters + string.digits
    session = ''.join(secrets.choice(alphabet) for i in range(10))
    pisa_error = False
    output = []
    if '_CHOIR_CorrectedChains' in os.path.basename(pdb):
        session_name = os.path.basename(pdb).split('_CHOIR_CorrectedChains')[0]
    else:
        session_name = os.path.basename(pdb).split('.')[0]
    pisa_cfg, pisa_tmp_dir = create_pisa_conf(os.getcwd(), session+chain)

    # Erase ocasional previous failed session
    pisa_cmd0 = [pisa_exe, session, '-erase', pisa_cfg]
    try:
        if verbosity:
            output.append(clrs['b']+'PISA'+clrs['n']+' command line: '+' '.join(pisa_cmd0))
        pisa_res = subprocess.Popen(pisa_cmd0, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        pisa_res.wait()
    except subprocess.CalledProcessError:
        pass

    # Analyse PDB file
    pisa_cmd1 = [pisa_exe, session, '-analyse', pdb, pisa_cfg]
    try:
        if verbosity:
            output.append(clrs['b']+'PISA'+clrs['n']+' command line: '+' '.join(pisa_cmd1))
        pisa_res = subprocess.Popen(pisa_cmd1, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        pisa_res.wait()
    except subprocess.CalledProcessError:
        output.append(clrs['r']+'PISA command '+' '.join(pisa_cmd1)+' returned non-zero status\n'+clrs['n'])
        pisa_error = True

    if gen_oligomer_data is True:

        # Generate Interfaces XML
        pisa_cmd2 = [pisa_exe, session, '-xml', 'interfaces', pisa_cfg]
        try:
            if verbosity:
                output.append(clrs['b']+'PISA'+clrs['n']+' command line: '+' '.join(pisa_cmd2))
            pisa_res = subprocess.Popen(pisa_cmd2, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            pisa_xml = pisa_res.stdout.read().decode()
            with open(session_name+chain+'_CHOIR_PisaInterfaces.xml', 'w') as f:
                f.write(pisa_xml)
        except subprocess.CalledProcessError:
            output.append(clrs['r']+'PISA command '+' '.join(pisa_cmd2)+' returned non-zero status\n'+clrs['n'])
            pisa_error = True
        except UnicodeDecodeError:
            output.append(clrs['r']+'PISA command '+' '.join(pisa_cmd2)+' output could not be properly decoded.\n'+clrs['n'])
            pisa_error = True


        # Generate Assemblies XML
        pisa_cmd3 = [pisa_exe, session, '-xml', 'assemblies', pisa_cfg]
        try:
            if verbosity:
                output.append(clrs['b']+'PISA'+clrs['n']+' command line: '+' '.join(pisa_cmd3))
            pisa_res = subprocess.Popen(pisa_cmd3, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            pisa_xml = pisa_res.stdout.read().decode()
            with open(session_name+chain+'_CHOIR_PisaAssemblies.xml', 'w') as f:
                f.write(pisa_xml)
        except subprocess.CalledProcessError:
            output.append(clrs['r']+'PISA command '+' '.join(pisa_cmd3)+' returned non-zero status\n'+clrs['n'])
            pisa_error = True
        except UnicodeDecodeError:
            output.append(clrs['r']+'PISA command '+' '.join(pisa_cmd3)+' output could not be properly decoded.\n'+clrs['n'])


    if gen_monomer_data is True:
        # Generate monomer data
        pisa_cmd4 = [pisa_exe, session, '-detail', 'monomer', '1', pisa_cfg]
        try:
            if verbosity:
                output.append(clrs['b']+'PISA'+clrs['n']+' command line: '+' '.join(pisa_cmd4))
            pisa_res = subprocess.Popen(pisa_cmd4, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            pisa_out = pisa_res.stdout.read().decode()
            monomer_data = session_name+chain+'_CHOIR_RSA.dat'
            with open(monomer_data, 'w') as f:
                f.write(pisa_out)
        except subprocess.CalledProcessError:
            output.append(clrs['r']+'PISA command '+' '.join(pisa_cmd4)+' returned non-zero status\n'+clrs['n'])
            pisa_error = True
        except UnicodeDecodeError:
            output.append(clrs['r']+'PISA command '+' '.join(pisa_cmd4)+' output could not be properly decoded.\n'+clrs['n'])
            pisa_error = True

    # Erase Session
    pisa_cmd5 = [pisa_exe, session, '-erase', pisa_cfg]
    try:
        if verbosity:
            output.append(clrs['b']+'PISA'+clrs['n']+' command line: '+' '.join(pisa_cmd5))
        subprocess.Popen(pisa_cmd5, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        pisa_res.wait()
    except subprocess.CalledProcessError:
        output.append(clrs['r']+'PISA command '+' '.join(pisa_cmd5)+' returned non-zero status\n'+clrs['n'])
        pisa_error = True

    try:
        shutil.rmtree(pisa_tmp_dir)
    except:
        pass

    if gen_monomer_data is True:
        return '\n'.join(output), pisa_error, monomer_data
    else:
        return '\n'.join(output), pisa_error


def parse_interfaces(interfaces_xml, candidate_chains, verbosity):
    output = []
    tree = et.parse(interfaces_xml)
    root = tree.getroot()
    interfaces_list = []
    for interface in root:
        # List os PISA-detected interactions
        interactions = []
        # Dictionary that will contain the details of the interface
        interface_dict = collections.OrderedDict()
        symop = []
        chain_pair = []
        hbondlist = []
        interfacing_residues = {}
        for item in interface:
            hbons_dict = {}
            if item.tag == 'int_area':
                int_area = round(float(item.text), 2)
            if item.tag == 'int_solv_en':
                int_solv_en = round(float(item.text), 2)
            if item.tag == 'h-bonds':
                for subitem in item:
                    if subitem.tag == 'n_bonds':
                        hbonds = int(subitem.text)
                    if subitem.tag == 'bond':
                        for detail in subitem:
                            if detail.tag == 'chain-1':
                                c1 = detail.text
                            if detail.tag == 'res-1':
                                r1 = detail.text
                            if detail.tag == 'seqnum-1':
                                s1 = detail.text
                            if detail.tag == 'atname-1':
                                a1 = detail.text
                            if detail.tag == 'chain-2':
                                c2 = detail.text
                            if detail.tag == 'res-2':
                                r2 = detail.text
                            if detail.tag == 'seqnum-2':
                                s2 = detail.text
                            if detail.tag == 'atname-2':
                                a2 = detail.text
                        interactions.append({'type':'HB','c1':c1,'r1':r1,'s1':s1,'a1':a1,'c2':c2,'r2':r2,'s2':s2,'a2':a2})
            if item.tag == 'salt-bridges':
                for subitem in item:
                    if subitem.tag == 'n_bonds':
                        saltbridges = int(subitem.text)
                    if subitem.tag == 'bond':
                        for detail in subitem:
                            if detail.tag == 'chain-1':
                                c1 = detail.text
                            if detail.tag == 'res-1':
                                r1 = detail.text
                            if detail.tag == 'seqnum-1':
                                s1 = detail.text
                            if detail.tag == 'atname-1':
                                a1 = detail.text
                            if detail.tag == 'chain-2':
                                c2 = detail.text
                            if detail.tag == 'res-2':
                                r2 = detail.text
                            if detail.tag == 'seqnum-2':
                                s2 = detail.text
                            if detail.tag == 'atname-2':
                                a2 = detail.text
                        interactions.append({'type':'SB','c1':c1,'r1':r1,'s1':s1,'a1':a1,'c2':c2,'r2':r2,'s2':s2,'a2':a2})
            if item.tag == 'ss-bonds':
                for subitem in item:
                    if subitem.tag == 'n_bonds':
                        ssbonds = int(subitem.text)
                    if subitem.tag == 'bond':
                        for detail in subitem:
                            if detail.tag == 'chain-1':
                                c1 = detail.text
                            if detail.tag == 'res-1':
                                r1 = detail.text
                            if detail.tag == 'seqnum-1':
                                s1 = detail.text
                            if detail.tag == 'atname-1':
                                a1 = detail.text
                            if detail.tag == 'chain-2':
                                c2 = detail.text
                            if detail.tag == 'res-2':
                                r2 = detail.text
                            if detail.tag == 'seqnum-2':
                                s2 = detail.text
                            if detail.tag == 'atname-2':
                                a2 = detail.text
                        interactions.append({'type':'SS','c1':c1,'r1':r1,'s1':s1,'a1':a1,'c2':c2,'r2':r2,'s2':s2,'a2':a2})
            if item.tag == 'molecule':
                interface_residues = []
                for subitem in item:
                    if subitem.tag == 'symop':
                        symop.append(subitem.text)
                    if subitem.tag == 'chain_id':
                        chain_id = subitem.text
                        chain_pair.append(subitem.text)
                    if subitem.tag == 'residues':
                        for residue in subitem:
                            if residue.tag == 'residue':
                                for detail in residue:
                                    if detail.tag == 'name':
                                        resname = detail.text
                                    if detail.tag == 'ser_no':
                                        resindex = detail.text
                                    if detail.tag == 'bsa':
                                        if float(detail.text) > 0:
                                            interface_residues.append(resname+resindex)
                interfacing_residues[chain_id] = interface_residues

        if symop:
            if symop[0] == 'x,y,z' and symop[1] == 'X,Y,Z':
                if all(chain in candidate_chains for chain in chain_pair):
                    interface_dict['chains'] = chain_pair
                    interface_dict['interface area'] = int_area
                    interface_dict['interface solvation energy'] = int_solv_en
                    interface_dict['hydrogen bonds'] = hbonds
                    interface_dict['salt bridges'] = saltbridges
                    interface_dict['disulphide bridges'] = ssbonds
                    interface_dict['interactions'] = interactions
                    for chain in chain_pair:
                        interface_dict['interface residues '+chain] = interfacing_residues[chain]
                    output.append(clrs['y']+' <> '.join(chain_pair)+clrs['n'])
                    output.append(clrs['y']+'Interface Area: '+clrs['n']+str(int_area)+' A^2')
                    output.append(clrs['y']+'Interface Solvation Energy: '+clrs['n']+str(int_solv_en)+' kcal/mol')
                    output.append(clrs['y']+'Hydrogen Bonds: '+clrs['n']+str(hbonds))
                    output.append(clrs['y']+'Salt Bridges: '+clrs['n']+str(saltbridges))
                    output.append(clrs['y']+'Disulphide Bridges: '+clrs['n']+str(ssbonds)+"\n")
                    interfaces_list.append(interface_dict)
    return interfaces_list, '\n'.join(output)

def run_gesamt(reference_name, reference_pdb, target_name, target_pdb, chain, args):
    output = []
    if chain is None:
        output.append('Running '+clrs['b']+'GESAMT'+clrs['n']+' to align '+clrs['y']+target_name+clrs['n']+' to '+clrs['y']+reference_name+clrs['n'])
        fasta_out = target_name+"_"+reference_name+'_CHOIR_Gesamt.fasta'
        gesamtcmd = [gesamt_exe, reference_pdb, target_pdb, '-a', fasta_out]
    else:
        output.append('Running '+clrs['b']+'GESAMT'+clrs['n']+' to align '+clrs['y']+target_name+clrs['n']+' to '+clrs['y']+reference_name+clrs['n']+' - Chain '+clrs['y']+chain+clrs['n'])
        fasta_out = target_name+"_"+reference_name+chain+'_CHOIR_Gesamt.fasta'
        gesamtcmd = [gesamt_exe, reference_pdb, '-s', chain, target_pdb, '-a', fasta_out]
    if args.force_single_core is True:
        gesamtcmd.append('-nthreads=1')
    else:
        gesamtcmd.append('-nthreads=auto')
    if args.verbosity == 1:
        output.append(clrs['b']+'GESAMT'+clrs['n']+' command line: '+' '.join(gesamtcmd))


    if 'qscore' in locals():
        del qscore
    if 'rmsd' in locals():
        del rmsd
    if 'id' in locals():
        del id
    if 'al_res' in locals():
        del al_res
    gesout = subprocess.check_output(gesamtcmd).decode('utf-8').split('\n')

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
        fasta_out = 'None'
        qscore = 0
        rmsd = 'NA'
    else:
        output.append('Done running '+clrs['b']+'GESAMT'+clrs['n']+'.\n')
        output.append('Alignment written to '+clrs['g']+os.path.basename(fasta_out)+clrs['n'])


    return qscore, rmsd, fasta_out, '\n'.join(output)


def run_molprobity(structure_file, args):
    output = []
    output.append('\nRunning '+clrs['b']+'Molprobity'+clrs['n']+' for '+clrs['y']+structure_file+clrs['n']+'\n')
    molprobity_cmd = [molprobity_exe, structure_file, 'prefix='+structure_file.split('.')[0]+'_molprobity']
    if args.verbosity == 1:
        output.append(clrs['b']+'Molprobity'+clrs['n']+' command line: '+' '.join(molprobity_cmd))
    molprobity_out = subprocess.check_output(molprobity_cmd).decode('utf-8').split('\n')
    molprobity_results = {}
    start_reading = False
    for line in molprobity_out:
        if start_reading is False:
            if line == '=================================== Summary ===================================':
                start_reading = True
        if start_reading is True:
            if 'Ramachandran outliers =' in line:
                molprobity_results['rama_out'] = float(line.strip().split()[3])
            if 'favored =' in line:
                molprobity_results['rama_fav'] = float(line.strip().split()[2])
            if 'Rotamer outliers      =' in line:
                molprobity_results['rot_out'] = float(line.strip().split()[3])
            if 'C-beta deviations     =' in line:
                molprobity_results['cb_dev'] = int(line.strip().split()[3])
            if 'Clashscore            =' in line:
                molprobity_results['clashscore'] = float(line.strip().split()[2])
            if 'MolProbity score      =' in line:
                molprobity_results['molprobity_score'] = float(line.strip().split()[3])
    output.append('Done running '+clrs['b']+'Molprobity'+clrs['n']+'\n')
    return molprobity_results, '\n'.join(output)

def pymol_screenshot_mono(monomer_structure, z_entropies, args):
    outfile = os.path.basename(monomer_structure).replace('.pdb','.png')
    prot =  os.path.basename(monomer_structure).replace('.pdb','')
    minval = min([i[1] for i in z_entropies.items()])
    maxval = max([i[1] for i in z_entropies.items()])
    pymolrc = 'pymolrc'
    pymol_cmd = [pymol_exe, '-c']
    if args.force_single_core is True:
        max_cores = '1'
    else:
        max_cores = str(args.available_cores)
    with open(pymolrc, 'w') as f:
        f.write(tw.dedent("""
                          set max_threads, """+max_cores+"""
                          load """+monomer_structure+"""
                          set ray_opaque_background,0
                          set ray_shadows, 0
                          set ray_trace_mode, 0
                          set ray_trace_color, 1
                          set antialias,', 1
                          set orthoscopic, 1
                          cmd.spectrum("b", "red_yellow_green", '"""+str(prot)+"""', minimum="""+str(minval)+""", maximum="""+str(maxval)+""")
                          cmd.ramp_new("ramp_obj", '"""+str(prot)+"""', range=["""+str(minval)+""", """+str(maxval)+"""], color="[red, yellow, green]")
                          orient
                          zoom complete=1
                          ray 2000, 2000
                          save """+outfile))
    try:
        subprocess.Popen(pymol_cmd, stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT).wait()
        os.remove(pymolrc)
        print('Protomer image generated: '+clrs['g']+os.path.basename(outfile)+clrs['n']+'\n')
        return './'+outfile
    except subprocess.CalledProcessError:
        os.remove(pymolrc)
        print('Failed to run Pymol to generate image! Does pymol executable exist?')
        print_sorry()


def pymol_screenshot(structure_file, args, putty=False):
    outfile1 = os.path.basename(structure_file).replace('.pdb', '.png')
    outfile2 = os.path.basename(structure_file).replace('.pdb', '_putty0.png')
    outfile3 = os.path.basename(structure_file).replace('.pdb', '_putty90.png')
    outfile4 = os.path.basename(structure_file).replace('.pdb', '_putty180.png')
    outfile5 = os.path.basename(structure_file).replace('.pdb', '_putty270.png')
    prot = os.path.basename(structure_file).replace('.pdb', '')
    pymol_script = os.path.basename(structure_file).replace('.pdb', '.pml')
    pymol_cmd = [pymol_exe, '-c', pymol_script]
    if args.multiprocess is True:
        if int(args.available_cores/args.models) >= 1:
            max_cores = str(int(args.available_cores/args.models))
        else:
            max_cores = '1'
    elif args.force_single_core is True:
        max_cores = '1'
    else:
        max_cores = str(int(args.available_cores))
    with open(pymol_script, 'w') as f:
        f.write(tw.dedent("""
                          set max_threads, """+max_cores+"""
                          load """+structure_file+"""
                          set ray_opaque_background,0
                          set ray_shadows, 0
                          set ray_trace_mode, 0
                          set ray_trace_color, 1
                          set antialias,', 1
                          set orthoscopic, 1
                          util.cbc
                          orient
                          zoom complete=1
                          ray 2000, 2000
                          save """+outfile1))
        if putty:
            f.write(tw.dedent("""
                              ray 2000, 2000
                              save """+outfile2+"""
                              rotate y, 90
                              ray 2000, 2000
                              save """+outfile3+"""
                              rotate y, 90
                              ray 2000, 2000
                              save """+outfile4+"""
                              rotate y, 90
                              ray 2000, 2000
                              save """+outfile5))
                              # preset.b_factor_putty(selection='all')
                              # set cartoon_putty_scale_power, 3
                              # util.cbc
    try:
        subprocess.Popen(pymol_cmd, stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT).wait()
        os.remove(pymol_script)

        if putty:
            output = 'Oligomer images generated:\n'+clrs['g']+'\n'.join([outfile1, outfile2, outfile3, outfile4, outfile5])+clrs['n']+'\n'
            return ['./'+outfile1, './'+outfile2, './'+outfile3, './'+outfile4, './'+outfile5], output
        else:
            output = 'Oligomer image generated: '+clrs['g']+outfile1+clrs['n']+'\n'
            return './'+outfile1, output
    except subprocess.CalledProcessError:
        os.remove(pymolrc)
        print('Failed to run Pymol to generate image! Does pymol executable exist?')
        print_sorry()

def create_pisa_conf(confdir, id):
    pisa_cfg_file = os.path.join(confdir, id+"_pisa.cfg")
    with open(pisa_cfg_file, 'w') as cf:
        cf.write(tw.dedent("""
    # ------------------------------------------------------------------- #
    #                                                                     #
    #          This is configuratrion file for PISA software.             #
    #                                                                     #
    #   When used in command-prompt mode, this file must be specified     #
    #       as last argument in thecommand line, or pointed out by        #
    #              PISA_CONF_FILE environmental variable.                 #
    #                                                                     #
    #    This file may be also used to configure the QtPISA graphical     #
    #   application, by either reading it using the "Load CFG" button   #
    #  in PISA Settings Dialog, or by running QtPISA from command prompt  #
    #    with this file as the only command-line argument. QtPISA needs   #
    #           to be configure only once after installation.             #
    #                                                                     #
    # ------------------------------------------------------------------- #


    #  DATA_ROOT must point on the directory to contain session
    #  directories.
    DATA_ROOT
    """+confdir+"""/"""+id+"""_pisa_tmp/


    #  SRS_DIR must point on the directory containing SRS files.
    SRS_DIR
    """+ccp4_base+"""/share/ccp4srs/


    #  MOLREF_DIR must point on the directory containing MolRef files.
    MOLREF_DIR
    """+ccp4_base+"""/share/pisa/


    #  PISTORE_DIR must point on the directory containing files:
    #          agents.dat
    #          asm_params.dat
    #          intfstats.dat
    #          syminfo_pisa.lib and
    #          rcsb_symops.dat
    PISTORE_DIR
    """+ccp4_base+"""/share/pisa/


    #  HELP_DIR must point on the directory containing HTML-formatted
    #  help files.
    HELP_DIR
    """+ccp4_base+"""/html/pisa/


    #  DOWNLOAD_URL is used only in QtPISA. It must give a valid base URL
    #  for downloading PDB's cif-formatted files. The full downlad URL is
    #  formed by appending lower-case PDB code with extension '.cif' to the
    #  base download URL.
    DOWNLOAD_URL
    http://files.rcsb.org/download/


    #  RASMOL_COM must give the rasmol launch command line
    RASMOL_COM
    /dummy/rasmol


    #  JMOL_COM must give path to jmol.jar
    JMOL_COM
    /dummy/jmol


    #  CCP4MG_COM must give the ccp4mg launch command line
    CCP4MG_COM
    """+ccp4_base+"""/bin/ccp4mg


    #  SESSION_PREFIX is prefix for the names of session directories
    #  created in DATA_PATH ("pisrv_" is used by default). Be sure to
    #  have unique prefixes for each configuration file that is invoked
    #  from a different user login or apache service. Session directories
    #  are regularly removed from DATA_PATH, and SESSION_PREFIX allows
    #  one to avoid permission conflicts between different services.
    SESSION_PREFIX
    __tmp_pisa_


    #  =========================================================================
    #     The following configuration parameters are needed only for jspisa,
    #   which is a web-server application. They may be ignored or even removed
    #        from configuration files for command-prompt pisa and qtpisa.
    #  =========================================================================


    #  PDB_DIR is used only by jspisa (web-server). It must
    #  give absolute path to PDB directory.
    PDB_DIR



    #  HELP_URI should point on the directory containing HTML-formatted
    #  help files.
    HELP_URI
    http://www.domain.com/path/to/pisahelp


    #  PDB_DIR_FORMAT is used only by jspisa (web-server). It must
    #  specify whether the PDB directory has plain structures (all files
    #  at root), or split (files found in sub-directories named by letters
    #  2 and 3 of the PDB code, small case), and whether in PDB or mmCIF
    #  format, gzipped or not. Permissible values include:
    #      PDB_PLAIN_PDB
    #      PDB_PLAIN_mmCIF
    #      PDB_PLAIN_PDB_GZ
    #      PDB_PLAIN_mmCIF_GZ
    #      PDB_SPLIT_PDB
    #      PDB_SPLIT_mmCIF
    #      PDB_SPLIT_PDB_GZ
    #      PDB_SPLIT_mmCIF_GZ
    PDB_DIR_FORMAT
    PDB_PLAIN_PDB_GZ


    #  PHP_URI is used only by jspisa (web-server). It must
    #  start with "http://" and point to directory with php scripts.
    PHP_URI
    http://www.pisa.com/php/


    #  JSRVIEW_URI is used only by jspisa (web-server). It must
    #  start with "http://" and point to directory with javascript
    #  support for jsrview api.
    JSRVIEW_URI
    http://www.pisa.com/jsrview/


    #  EXPIRY_TIME is used only by jspisa (web-server). It gives
    #  time (integer hours) for session to expire. Expired sessions are
    #  erased when necessary unless they are locked.
    EXPIRY_TIME
    10800


    #  ERASE_TIME is used only by jspisa (web-server). It gives
    #  time (integer hours) for session to be erased unconditionally,
    #  whether locked or not.
    ERASE_TIME
    172800
                           """))
    pisa_tmp_dir = os.path.join(confdir, id+'_pisa_tmp')
    os.mkdir(pisa_tmp_dir)
    return pisa_cfg_file, pisa_tmp_dir


def html_report(report, args):
    html_out = report['model_oligomer_name']+'_CHOIR_Report.html'
    templateLoader = jinja2.FileSystemLoader(searchpath=os.path.join(choir_path, 'Contents'))
    templateEnv = jinja2.Environment(loader=templateLoader)
    report_template = "report_template.html"
    template = templateEnv.get_template(report_template)
    outputText = template.render(protchoir_logo=os.path.join(choir_path, 'Contents', 'ProtCHOIRReport.svg'), report=report)

    with open(html_out, 'w') as f:
        f.write(outputText)
    print('ProtCHOIR report available in: '+clrs['g']+os.path.basename(html_out)+clrs['n']+'\n')
    return html_out

def get_areas(monomer_data):
    surface_residues = collections.OrderedDict()
    start_reading = False
    first_line = False
    for line in open(monomer_data, 'r'):
        #while not start_reading:
        if line == ' -----+------------+---------------------------------------------\n':
            #print(line)
            start_reading = True
            first_line = True
        elif line == " -----'------------+---------------------------------------------\n":
            start_reading = False
        if start_reading:
            if not first_line:
                line = line.replace('\n', '')
                splitline = line.split(' | ')
                res = splitline[1].split(':')[1].replace(' ', '')
                bsa = float(splitline[2].split()[1])
                asa = float(splitline[2].split()[0]) - bsa
                if asa > maxasa[res[:3]]:
                    asa = maxasa[res[:3]]
                rsa = round(asa/maxasa[res[:3]]*100, 2)
                surface_residues[res] = (asa, bsa, float(rsa))
        first_line = False
    return surface_residues


def map_residue_index(surface_residues):
    residue_index_mapping = collections.OrderedDict()
    for i, (res, rsa) in enumerate(surface_residues.items()):
        residue_index_mapping[i+1] = int(res[3:])
    return residue_index_mapping


def plot_analysis(pdb_name, surface_residues, entropies, z_entropies, tmdata, args, minx=None, maxx=None):
    output = []
    # Reset Matplotlib parameters to default
    matplotlib.rcParams.update(matplotlib.rcParamsDefault)

    # Plot RSA and calculate total hydrophobic area exposed
    x = []
    for res, areas in surface_residues.items():
        x.append(int(res[3:]))

    if minx is None and maxx is None:
        minx = min(x)
        maxx = max(x)

    for i in range(minx, maxx+1):
        if i not in x:
            surface_residues['XXX'+str(i)] = (-2, -2, -2)
    x = []
    y = []
    z = []
    total_area = 0
    hydrophobic_area = 0
    for res, areas in surface_residues.items():
        if areas[0] > 0:
            total_area += areas[0]
        if areas[2] == -2:
            color = 'grey'
        elif hydro[res[:3]] < hydro['ALA']:
            color = 'blue'
        elif hydro[res[:3]] > hydro['GLY']:
            color = 'red'
            if int(res[3:]) not in tmdata[1]:
                hydrophobic_area += float(areas[0])
            else:
                if args.verbosity == 1:
                    output.append('Ígnoring hydrophobic exposed area from residue: '+res)
        else:
            color = 'black'
        x.append(int(res[3:]))
        y.append(float(areas[2]))
        z.append(color)
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    output.append('\nProtomer surface exposure data:')
    output.append('Total area exposed= '+str(round(total_area, 2))+' A^2')
    output.append('Hydrophobic area exposed = '+str(round(hydrophobic_area, 2))+' A^2')

    # Claculate total conserved area exposed
    conserved_area = 0
    if not args.skip_conservation:
        for (res, areas), (column, zscore) in zip(surface_residues.items(), z_entropies.items()):
            if zscore > 1:
                if int(res[3:]) not in tmdata[1]:
                    conserved_area += float(areas[0])
                else:
                    if args.verbosity == 1:
                        output.append('Ígnoring conserved exposed area from residue: '+res)
        output.append('Conserved area exposed= '+str(round(conserved_area, 2))+' A^2')

    # Plot everything

    p, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, figsize=(18, 9), gridspec_kw={'height_ratios': [1, 1, 0.5, 3.5]})
    plt.suptitle('Analysis of '+pdb_name+' Protomer', fontsize=20, fontweight='bold')


    # Plot Relative Entropy
    ax1.set_ylabel("Relative\nEntropy (bits)"+r"$^*$", fontsize=12)
    ax1.get_yaxis().set_label_coords(-0.02, 0.5)
    if not args.skip_conservation:
        ax1.bar(entropies.keys(), entropies.values(), color='grey', zorder=3)
    ax1.grid(True, linestyle=':', linewidth=0.7, zorder=0, color='k')


    # Plot Entropy Z-Scores
    ax2.set_ylabel("Z Scores", fontsize=12)
    ax2.get_yaxis().set_label_coords(-0.02, 0.5)
    if not args.skip_conservation:
        colors = []
        for column, zscore in z_entropies.items():
            if zscore > 0:
                color = 'darkgreen'
            else:
                color = 'firebrick'
            colors.append(color)

        ax2.bar(z_entropies.keys(), z_entropies.values(), color=colors, zorder=3)
    ax2.grid(True, linestyle=':', linewidth=0.7, zorder=0, color='k')

    # Plot membrane residues
    ax3.set_ylabel("TM Helices", fontsize=12)
    ax3.get_yaxis().set_label_coords(-0.02, 0.5)
    ax3.get_yaxis().set_ticks([])
    ax3.get_yaxis().set_ticklabels([])
    for seg in tmdata[0]:
        ax3.plot( seg, [0, 0],  marker=None, linewidth=12, color='goldenrod')
    ax3.grid(True, linestyle=':', linewidth=0.7, zorder=0, color='k')

    # Plot RSA and Hydrophobicity
    df = pd.DataFrame({'x': x, 'y': y})
    ax4.bar('x', 'y', data=df, color=z, zorder=3)
    ax4.xaxis.set_ticks(np.arange(min(x), max(x), round((max(x)-min(x))*0.025)))
    ax4.set_xlabel('Residue Index', fontsize=16)
    ax4.set_ylabel('Relative Exposed Surface Area (%)', fontsize=12)
    ax4.get_yaxis().set_label_coords(-0.02, 0.5)
    ax4.set_xlim(min(x)-2, max(x)+2)
    ax4.axhline(y=20, linestyle='-', linewidth=0.7, color='k')
    ax4.grid(True, linestyle=':', linewidth=0.7, zorder=0, color='k')
    # Draw legend
    legend_elements = [Line2D([0], [0], color='b', lw=10, label='Hydrophylic'),
                       Line2D([0], [0], color='k', lw=10, label='Mid-range'),
                       Line2D([0], [0], color='r', lw=10, label='Hydrophobic'),
                       Line2D([0], [0], color='grey', lw=10, label='Gap')]
    ax4.legend(handles=legend_elements, loc=9, ncol=4, bbox_to_anchor=(0.5, -0.10), frameon=False)

    # Organize
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    p.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in p.axes[:-1]], visible=False)
    p.text(0.04, 0.04, r'$^*$'+'Entropy loss compared to background distribution', fontsize=10, ha='left')

    # Save plot
    outfile = pdb_name+'_CHOIR_ProtomerPlot.png'
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    # Close figure
    plt.close()
    output.append('\nCombining conservation and surface exposure data...')
    output.append('Analysis plots for '+pdb_name+' generated : '+clrs['g']+os.path.basename(outfile)+clrs['n']+'\n')
    return './'+os.path.basename(outfile), round(total_area, 2), round(hydrophobic_area, 2), round(conserved_area, 2), minx, maxx, '\n'.join(output)


def plot_entropy_only(pdb_name, entropies, z_entropies, tmdata, args):

    # Reset Matplotlib parameters to default
    matplotlib.rcParams.update(matplotlib.rcParamsDefault)
    print('\nPlotting entropy scores per residue...')
    # Plot RSA and calculate total hydrophobic area exposed
    # Plot everything

    p, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize=(18, 9), gridspec_kw={'height_ratios': [1, 2, 1]})
    plt.suptitle('Analysis of '+pdb_name+' Protomer', fontsize=24, fontweight='bold')


    # Plot Relative Entropy
    ax1.set_ylabel("Relative\nEntropy (bits)"+r"$^*$", fontsize=12)
    ax1.get_yaxis().set_label_coords(-0.02, 0.5)
    if not args.skip_conservation:
        ax1.bar(entropies.keys(), entropies.values(), color='grey', zorder=3)
    ax1.grid(True, linestyle=':', linewidth=0.7, zorder=0, color='k')

    # Plot Entropy Z-Scores
    ax2.set_ylabel("Z Scores", fontsize=12)
    ax2.get_yaxis().set_label_coords(-0.02, 0.5)
    if not args.skip_conservation:
        colors = []
        for column, zscore in z_entropies.items():
            if zscore > 0:
                color = 'darkgreen'
            else:
                color = 'firebrick'
            colors.append(color)

        ax2.bar(z_entropies.keys(), z_entropies.values(), color=colors, zorder=3)
    ax2.grid(True, linestyle=':', linewidth=0.7, zorder=0, color='k')

    # Plot membrane residues
    ax3.set_ylabel("TM\nHelices", fontsize=12)
    ax3.get_yaxis().set_label_coords(-0.02, 0.5)
    ax3.get_yaxis().set_ticks([])
    ax3.get_yaxis().set_ticklabels([])
    for seg in tmdata[0]:
        ax3.plot( seg, [0, 0],  marker=None, linewidth=12, color='goldenrod')
    ax3.grid(True, linestyle=':', linewidth=0.7, zorder=0, color='k')

    # Organize
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    p.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in p.axes[:-1]], visible=False)
    p.text(0.04, 0.04, r'$^*$'+'Entropy loss compared to background distribution', fontsize=10, ha='left')

    # Save plot
    outfile = pdb_name+'_CHOIR_ProtomerPlot.png'
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    # Close figure
    plt.close()
    print('Conservation plots for '+pdb_name+' generated : '+clrs['g']+os.path.basename(outfile)+clrs['n']+'\n')
    return './'+os.path.basename(outfile)


# Execute
###############################################################################
