# Imports
###############################################################################
import os
import re
import sys
import time
import gzip
import pickle
import shutil
import parasail
import traceback
import subprocess
import textwrap as tw
from ProtCHOIR.Initialise import *
import ProtCHOIR.Toolbox as pctools
from progressbar import progressbar as pg

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
This file conatins functions used to update and curate the databases that are
used by ProtCHOIR.
'''
# Global Variables
###############################################################################
pdb_server = 'rsync.wwpdb.org::ftp'
pdb_subdir = '/data/structures/divided/pdb/'
pdb_subdir1 = '/data/biounit/PDB/divided/'
seqres_ftp = 'ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz'
port = '33444'
uniref50_ftp = 'ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz'
# Functions
###############################################################################
def create_directories():
    if not os.path.isdir(choirdb):
        os.mkdir(choirdb)
    if not os.path.isdir(pdb_archive):
        os.mkdir(pdb_archive)
    if not os.path.isdir(pdb1_archive):
        os.mkdir(pdb1_archive)
    if not os.path.isdir(ges_homo_archive):
        os.mkdir(ges_homo_archive)
    if not os.path.isdir(pdb_homo_archive):
        os.mkdir(pdb_homo_archive)

def update_pdb(verbosity):
    '''
    Runs rsync to update the local pdb database.
    Called by: update_databases()
    '''
    if verbosity == 1:
        vflag = '-vv'
    else:
        vflag = ''
    subprocess.run(['rsync',
                    '-rlptPW',
                    vflag,
                    '-z',
                    '--delete',
                    '--port='+port,
                    pdb_server+pdb_subdir,
                    pdb_archive])

def update_pdb1(verbosity):
    '''
    Runs rsync to update the local pdb database.
    Called by: update_databases()
    '''
    if verbosity == 1:
        vflag = '-vv'
    else:
        vflag = ''
    subprocess.run(['rsync',
                    '-rlptPW',
                    vflag,
                    '-z',
                    '--delete',
                    '--port='+port,
                    pdb_server+pdb_subdir1,
                    pdb1_archive])


def update_uniref(verbosity):
    '''
    Runs wget to update the local uniref50 database, decompresses it and runs
    makeblastdb.
    Called by: update_databases()
    '''
    uniref50_fasta = os.path.join(choirdb,'uniref50/uniref50.fasta')
    pctools.printv('Fetching uniref50.fasta...', verbosity)
    attempt = 0
    while attempt < 3:
        try:
            wgetout = subprocess.check_output(['wget', '-m',
                                               '-r', '-nH',
                                               '--cut-dirs=4',
                                               '--user=anonymous',
                                               uniref50_ftp,
                                               '-P', choirdb],
                                               stderr=subprocess.STDOUT)
            break
        except:
            attempt += 1
            if attempt < 3:
                print('Attempt '+str(attempt)+' failed, trying again.')
            if attempt == 3:
                print('Failed to download UniRef50 in 3 attempts. Try again later.')

    no_wget = 'uniref50.fasta.gz’ -- not retrieving'

    if no_wget not in wgetout.decode('UTF-8') or not os.path.isfile(uniref50_fasta):
        pctools.printv('Decompressing uniref50.fasta...', verbosity)

        with gzip.open(uniref50_fasta+'.gz', 'rb') as fin, open(uniref50_fasta, 'wb') as fout:
            shutil.copyfileobj(fin,fout)
    if no_wget not in wgetout.decode('UTF-8') or not os.path.isfile(uniref50_fasta+'.pal'):
        subprocess.run([makeblastdb_exe,'-in',
                       uniref50_fasta,
                       '-parse_seqids', '-dbtype', 'prot',
                       '-out', uniref50])

#def split_uniref(verbosity):
def update_seqres(verbosity):
    '''
    Runs wget to update the local seqres database, decompresses it and runs
    makeblastdb.
    Called by: update_databases()
    '''
    seqres_dir = os.path.join(choirdb,'seqres')
    if not os.path.isdir(seqres_dir):
        os.mkdir(seqres_dir)
    seqres_txt = os.path.join(seqres_dir,'pdb_seqres.txt')
    seqres_fasta = os.path.join(seqres_dir,'seqres.fasta')
    pctools.printv('Fetching pdb_seqres.txt...', verbosity)
    attempt = 0
    while attempt < 3:
        try:
            wgetout = subprocess.check_output(['wget', '-m',
                                               '-r', '-nH',
                                               '--cut-dirs=3',
                                               '--user=anonymous',
                                               seqres_ftp,
                                               '-P', seqres_dir],
                                               stderr=subprocess.STDOUT)
            break
        except:
            attempt += 1
            if attempt < 3:
                print('Attempt '+str(attempt)+' failed, trying again.')
            if attempt == 3:
                print('Failed to download seqres in 3 attempts. Try again later.')

    no_wget = 'seqres.txt.gz’ -- not retrieving'

    if no_wget not in wgetout.decode('UTF-8') or not os.path.isfile(seqres_txt):
        pctools.printv('Decompressing pdb_seqres.txt...', verbosity)

        with gzip.open(seqres_txt+'.gz', 'rb') as fin, open(seqres_fasta, 'wb') as fout:
            shutil.copyfileobj(fin,fout)
    if no_wget not in wgetout.decode('UTF-8') or not os.path.isfile(seqres_fasta+'.pal'):
        subprocess.run([makeblastdb_exe,'-in',
                       seqres_fasta,
                       '-parse_seqids', '-dbtype', 'prot', '-blastdb_version', '5',
                       '-out', seqres])



def read_latest_assession(stats_dir):
    '''
    Looks in the stats directory for the latest assession dat file (...-choirdb.dat)
    reads the file and returns a dictionary containing the pdb files as keys and
    assession times as values
    Called by: curate_homoDB()
    '''
    previous_assessions = {}
    for file in os.listdir(stats_dir):
        file = os.path.join(stats_dir,file)
        if file.endswith('choirdb.dat'):
            date = os.path.getctime(file)
            previous_assessions[file] = date
    if previous_assessions:
        latest_assession = max(previous_assessions, key=previous_assessions.get)
        asseession_log = {}
        with open(latest_assession, 'r') as f:
            for line in f:
                items = line.split()
                code, date = items[0], items[1]
                asseession_log[code] = date
    else:
        asseession_log = {}
    return asseession_log

def list_new_files(pdb1_archive, assession_log, verbosity):
    '''
    Taps into the pdb1 local repository and checks if there are new files there
    which should be assessed by curate_homoDB function. It creates a list with
    files from the pdb1 database that are newer than the ones last assessed
    (registered in the dat file); which means that it takes the previous
    assession log as an input in the form of a dictionary where keys are files
    and the values correspond to the last assession time.
    Called by: curate_homoDB()
    '''
    new_files = []
    pctools.printv('Assessing files in PDB1 archive...', verbosity)
    assert os.path.isdir(pdb1_archive), clrs['r']+'\n\n Not able to find PDB archive.\n\n Does "'+pdb1_archive+'" exist?'+clrs['n']
    pdbfiles = [os.path.join(dp, f) for dp, dn, filenames in os.walk(pdb1_archive) for f in filenames if f.endswith(".pdb1.gz")]
    for f in pdbfiles:
        filename = f.split('/')[-1]
        mod_date = os.path.getctime(f)
        if filename not in assession_log or mod_date > float(assession_log[filename]):
            pctools.printv(clrs['y']+f+' should be assessed'+clrs['n']+'...\n', verbosity)
            new_files.append(f)
    return new_files


def record_fasta(pdb_code, seqs, chain_ids, subfolder, type=None):
    if not os.path.isdir(os.path.join(pdb_homo_archive, subfolder)):
        os.mkdir(os.path.join(pdb_homo_archive, subfolder))
    type_folder = os.path.join(pdb_homo_archive, subfolder, type+'_sequences')
    if not os.path.isdir(type_folder):
        os.mkdir(type_folder)
    fasta_file = os.path.join(type_folder, pdb_code+".fasta")
    with open(fasta_file, 'w+') as f:
        for seq, chain_id in zip(seqs, chain_ids):
            if pctools.is_valid_sequence(seq[1]):
                wrapped_seq = "\n".join(tw.wrap(seq[1]))
                fasta_entry = '>'+pdb_code+':'+str(chain_id)+'\n'+wrapped_seq+'\n\n'
                f.write(fasta_entry)

def curate_homoDB(verbosity):
    '''
    Creates homo-oligomeric database from a local pdb repsitory.
    The divided scheme adopted by RCSB, in which the subdirectories
    are the two middle characters in the PDB code, is assumed.
    Each database contains three key files: dat, log and fasta.
    * homodb.dat contains only the pdb codes contained in the database.
    * homodb.log contains summarized relevant information about each entry.
    * homodb.fasta contains the sequences of every chain in the database.
    Called by: update_databases()
    '''
    # Create stats folder if does not exist
    stats_dir = os.path.join(pdb_homo_archive, 'stats')
    if not os.path.isdir(stats_dir):
        os.mkdir(stats_dir)
    # Compare latest assession with new files
    assession_log = read_latest_assession(stats_dir)
    new_files = list_new_files(pdb1_archive, assession_log, verbosity)
    print(clrs['g']+str(len(new_files))+clrs['n']+' new structure files were found and will be processed')
    now = str(time.strftime("%d-%m-%Y@%H.%M.%S"))
    dat_file = os.path.join(stats_dir, now+'-choirdb.dat')
    log_file = os.path.join(stats_dir, now+'-choirdb.log')
    err_file = os.path.join(stats_dir, now+'-choirdb.err')
    if not os.path.isfile(dat_file):
        with open(dat_file, 'w+'):
            pass
    # Write files not to be updated to new dat file
    with open(dat_file, 'a') as f:
        for i in assession_log:
            if i not in new_files:
                f.write(i+" "+assession_log[i]+"\n")
    # Create log file
    if not os.path.isfile(log_file):
        with open(log_file, 'w+') as f:
            f.write('Code, Chains, Author, Software, Date\n')

    # Read Chain correspondences
    chain_correspondences_file = os.path.join(stats_dir, 'chain_correspondences.pickle')
    if os.path.isfile(chain_correspondences_file):
        with open(chain_correspondences_file, 'rb') as p:
            chain_correspondences = pickle.load(p)
    else:
        chain_correspondences = {}

    # Main loop that will populate the ProtCHOIR database
    for pdb in pg(new_files, widgets=widgets):
        filename = pdb.split('/')[-1]
        subfolder = pdb.split('/')[-2]
        # Record assessment in dat file
        with open(dat_file, 'a') as f:
            f.write(filename+" "+str(time.time())+'\n')
        # Start assession
        pctools.printv('\nAssessing '+pdb+'...', verbosity)
        # Reject files larger than 10Mb
        file_size = os.stat(pdb).st_size / 1048576
        pctools.printv('File size: '+clrs['c']+'{0:.1g}'.format(file_size)+' Mb'+clrs['n'], verbosity)
        if file_size > 2:
            pctools.printv(clrs['r']+"File size too large!"+clrs['n'], verbosity)
            pctools.printv(clrs['y']+"Will try to fetch sequences from asymmetric unit."+clrs['n'], verbosity)
            try:
                alternative_pdb = os.path.join(pdb_archive, subfolder, 'pdb'+filename.split('.')[0]+'.ent.gz')
                pdb_code, structure, nchains = pctools.parse_pdb_structure(alternative_pdb)
                structure, chain_correspondences[pdb_code] = pctools.split_states(structure)
                nchainspostsplit, seqs, chain_ids = pctools.extract_seqs(structure, 0)
                # Write in fasta file
                pctools.printv(clrs['y']+"Recording large-pdb sequence"+clrs['n'], verbosity)
                record_fasta(pdb_code, seqs, chain_ids, subfolder, type='largepdb')
            except:
                pctools.printv(clrs['r']+"Failed to fetch sequence!"+clrs['n'], verbosity)
            continue

        try:
            pdb_code, structure, nchains = pctools.parse_pdb_structure(pdb)
            pctools.printv('Number of chains in structure '+clrs['y']+pdb_code+clrs['n']+': '+str(nchains), verbosity)
            # Reject structures with more than 60 chains
            if int(nchains) > 60:
                pctools.printv("Number of chains ("+clrs['y']+str(nchains)+clrs['n']+") larger than 60! "+clrs['r']+"Too many chains!"+clrs['n'], verbosity)
                pctools.printv(clrs['y']+"Will try to fetch sequences anyway."+clrs['n'], verbosity)
                try:
                    pdb_code, structure, nchains = pctools.parse_pdb_structure(pdb)
                    structure, chain_correspondences[pdb_code] = pctools.split_states(structure)
                    nchainspostsplit, seqs, chain_ids = pctools.extract_seqs(structure, 0)
                    pctools.printv(clrs['y']+"Recording large-pdb sequence"+clrs['n'], verbosity)
                    # Write in fasta file
                    record_fasta(pdb_code, seqs, chain_ids, subfolder, type='largepdb')
                except:
                    pctools.printv(clrs['r']+"Failed to fetch sequence!"+clrs['n'], verbosity)
                continue

            structure, chain_correspondences[pdb_code] = pctools.split_states(structure)
            nchainspostsplit, seqs, chain_ids = pctools.extract_seqs(structure, 0)
            pctools.printv('Number of chains ('+clrs['c']+str(nchains)+clrs['n']+') and file size ('+clrs['c']+str(file_size)+clrs['n']+') OK.'+clrs['g']+' Proceeding.'+clrs['n']+'\n', verbosity)
            # Try to get info from the canonic pdb header (homonimous to pdb1)
            canonpdb = "pdb"+pdb_code+".ent.gz"
            try:
                contents = pctools.parse_pdb_contents(os.path.join(pdb_archive, subfolder, canonpdb))[1]
            except:
                pctools.printv(clrs['r']+'\n\n Mismatch between pdb and biounit entries...'+clrs['n'], verbosity)
            author, software = pctools.get_annotated_states(contents)
            pctools.printv('Author determined biological unit = '+str(author), verbosity)
            pctools.printv('Software determined quaternary structure= '+str(software), verbosity)
            # Start assessing sequences and structures (from 2 up to 26 chains)
            if 1 < int(nchains) < 61:
                ids, proteinpair = pctools.get_pairwise_ids(seqs, nchains)
                for id in ids:
                    if id[0] >= 90:
                        color = clrs['g']
                    else:
                        color = clrs['r']
                    pctools.printv('Identity between chains '+clrs['y']+str(id[1])+clrs['n']+' and '+clrs['y']+str(id[2])+clrs['n']+' is '+color+str(id[0])+"%"+clrs['n']+".", verbosity)
                # Save records for pure homo-oligomers
                if all(id[0] > 90 for id in ids) and proteinpair is True:
                    pctools.printv("All identities over 90%. Likely "+clrs['b']+"homo-oligomeric"+clrs['n']+".", verbosity)
                    pctools.printv(clrs['y']+"FETCHING"+clrs['n']+".\n", verbosity)
                    # Write file to database
                    newfile = os.path.join(pdb_homo_archive, subfolder, pdb_code+".pdb")
                    if not os.path.isdir(os.path.join(pdb_homo_archive, subfolder)):
                        os.mkdir(os.path.join(pdb_homo_archive, subfolder))
                    io.set_structure(structure)
                    io.save(newfile)
                    pctools.gzip_pdb(newfile)
                    # Write to log file
                    with open(log_file, 'a') as f:
                        f.write(str(pdb_code)+","+str(nchains)+","+'/'.join(author)+","+'/'.join(software)+","+str(os.path.getctime(newfile+'.gz'))+'\n')
                    # Write in fasta file
                    pctools.printv(clrs['y']+"Recording homo-oligomer sequence."+clrs['n'], verbosity)
                    record_fasta(pdb_code, seqs, chain_ids, subfolder, type='homo')

                # Investigate partial homo-oligomers
                elif any(id[0] > 90 for id in ids) and proteinpair is True:
                    at_least_one_interface = False
                    for id in ids:
                        if id[0] > 90:
                            # Check if similar chains share interfaces
                            if pctools.check_interfaces(structure, id[1], id[2]):
                                at_least_one_interface = True
                                pctools.printv('Contacts found between chains '+clrs['g']+str(id[1])+clrs['n']+' and '+clrs['g']+str(id[2])+clrs['n']+' sharing '+clrs['g']+str(id[0])+clrs['n']+" % identity.", verbosity)
                                pctools.printv("At least one putative "+clrs['b']+"homo-oligomeric "+clrs['n']+"interface found.", verbosity)
                                pctools.printv(clrs['y']+"FETCHING"+clrs['n']+".\n", verbosity)
                                # Write file to database
                                newfile = os.path.join(pdb_homo_archive, subfolder, pdb_code+".pdb")
                                if not os.path.isdir(os.path.join(pdb_homo_archive, subfolder)):
                                    os.mkdir(os.path.join(pdb_homo_archive, subfolder))
                                io.set_structure(structure)
                                io.save(newfile)
                                pctools.gzip_pdb(newfile)
                                # Write to log file
                                with open(log_file, 'a') as f:
                                    f.write(str(pdb_code)+","+str(nchains)+","+'/'.join(author)+","+'/'.join(software)+","+str(os.path.getctime(newfile+'.gz'))+'\n')
                                # Write in fasta file
                                pctools.printv(clrs['y']+"Recording homo-oligomer sequence."+clrs['n'], verbosity)
                                record_fasta(pdb_code, seqs, chain_ids, subfolder, type='homo')

                                break
                    if at_least_one_interface is False:
                        pctools.printv("No homo-oligomeric interface found. Likely "+clrs['r']+"hetero-oligomeric"+clrs['n']+".", verbosity)
                        pctools.printv(clrs['y']+"Recording hetero-oligomer sequence"+clrs['n'], verbosity)
                        # Write in fasta file
                        record_fasta(pdb_code, seqs, chain_ids, subfolder, type='hetero')

                elif proteinpair is False:
                    pctools.printv(clrs['r']+"No proteic chain pairs found"+clrs['n']+".", verbosity)
                    if any([set(seq[1]) != {'X'} for seq in seqs]):
                        pctools.printv(clrs['y']+"Protein sequences found though"+clrs['n'], verbosity)
                        pctools.printv(clrs['y']+"Recording hetero-oligomer sequence"+clrs['n'], verbosity)
                        # Write in fasta file
                        record_fasta(pdb_code, seqs, chain_ids, subfolder, type='hetero')
                    else:
                        pctools.printv(clrs['r']+"Not even a single protein chain. Disregarding."+clrs['n'], verbosity)

                else:
                    pctools.printv("No similar chains found. Likely "+clrs['r']+"hetero-oligomeric"+clrs['n']+".", verbosity)
                    pctools.printv(clrs['y']+"Recording hetero-oligomer sequence"+clrs['n'], verbosity)
                    record_fasta(pdb_code, seqs, chain_ids, subfolder, type='hetero')

            elif int(nchains) == 1:
                pctools.printv("Only one chain found. Likely "+clrs['r']+"monomeric"+clrs['n']+".", verbosity)
                pctools.printv(clrs['y']+"Recording monomer sequence."+clrs['n'], verbosity)
                structure, chain_correspondences[pdb_code] = pctools.split_states(structure)
                nchains, seqs, chain_ids = pctools.extract_seqs(structure, 0)
                record_fasta(pdb_code, seqs, chain_ids, subfolder, type='mono')


        except:
            errtype, errvalue, errtraceback = sys.exc_info()
            errtypeshort = str(errtype).split('\'')[1]
            pctools.printv(clrs['r']+'*'+str(errtypeshort)+': '+str(errvalue)+ ' l.'+str(errtraceback.tb_lineno)+'*'+clrs['n'], verbosity)
            traceback.print_exception(*sys.exc_info())
            if errtypeshort == 'KeyboardInterrupt':
                quit()
            #pctools.printv(clrs['r']+"UNKNOWN FAULT"+clrs['n']+".", verbosity)
            if not os.path.isfile(err_file):
                with open(err_file,'w+') as f:
                    pass
            with open(err_file,'a') as f:
                f.write(filename+'\n')
            continue

    with open(chain_correspondences_file, 'wb') as p:
        pickle.dump(chain_correspondences, p, protocol=pickle.HIGHEST_PROTOCOL)

    if not os.path.isfile(err_file):
        with open(err_file,'w+') as f:
            f.write('\nNo errors. Assessment terminated succesfully.\n')


def collect_fasta(verbosity):
    '''
    Fetches fasta files in the pdb_homo_archive and creates a single fasta file
    within a "sequences" folder. For that, it checks the identity among the
    chains in the original fasta and only keeps track of the unique chains, i.e.
    less than 99% identity to the other chains. This file is later use to make
    the blast database.
    Called by: update_databases()
    '''
    fastafiles = [os.path.join(dp, f) for dp, dn, filenames in os.walk(pdb_homo_archive) for f in filenames if f.endswith(".fasta")]
    seqdir = os.path.join(pdb_homo_archive,'sequences')
    if not os.path.isdir(seqdir):
        os.mkdir(seqdir)

    largepdb_collected_fasta = os.path.join(seqdir,'largepdb_collected.fastas')
    with open(largepdb_collected_fasta, 'w+'):
        pass

    homo_collected_fasta = os.path.join(seqdir,'homo_collected.fastas')
    with open(homo_collected_fasta, 'w+'):
        pass

    mono_collected_fasta = os.path.join(seqdir,'mono_collected.fastas')
    with open(mono_collected_fasta, 'w+'):
        pass

    hetero_collected_fasta = os.path.join(seqdir,'hetero_collected.fastas')
    with open(hetero_collected_fasta, 'w+'):
        pass

    for fasta in pg(fastafiles, widgets=widgets):
        pctools.printv('Assessing '+clrs['y']+fasta+clrs['n']+'...', verbosity)
        contents = open(fasta, 'r').read()
        contentlines = contents.split('>')
        nchains = str(len(re.findall('>',contents)))
        pctools.printv('With '+clrs['y']+nchains+clrs['n']+' chains to be assessed\n', verbosity)
        uniques = []
        for entry in contentlines:
            if entry:
                splitentry = entry.split('\n', 1)
                pdbch = splitentry[0]
                seq = splitentry[1].replace('\n', '')
                if uniques:
                    percent_ids = []
                    for unique in uniques:
                        alignment = parasail.sg_stats_striped_16(seq, unique[1], 10, 1, parasail.blosum62)
                        if alignment.length == 0:
                            percent_ids.append(0)
                        else:
                            percent_ids.append((alignment.matches)/alignment.length*100)
                    if all(percent_id <= 99 for percent_id in percent_ids):
                        uniques.append([pdbch,seq])
                else:
                    uniques.append([pdbch,seq])

        if '/largepdb_sequences/' in fasta:
            with open(largepdb_collected_fasta, 'a') as f:
                for unique in uniques:
                    wrapped_seq = "\n".join(tw.wrap(unique[1]))
                    fasta_entry = '>'+unique[0]+'\n'+wrapped_seq+'\n\n'
                    f.write(fasta_entry)

        elif '/mono_sequences/' in fasta:
            with open(mono_collected_fasta, 'a') as f:
                for unique in uniques:
                    wrapped_seq = "\n".join(tw.wrap(unique[1]))
                    fasta_entry = '>'+unique[0]+'\n'+wrapped_seq+'\n\n'
                    f.write(fasta_entry)

        elif '/hetero_sequences/' in fasta:
            with open(hetero_collected_fasta, 'a') as f:
                for unique in uniques:
                    wrapped_seq = "\n".join(tw.wrap(unique[1]))
                    fasta_entry = '>'+unique[0]+'\n'+wrapped_seq+'\n\n'
                    f.write(fasta_entry)

        elif '/homo_sequences/' in fasta:
            with open(homo_collected_fasta, 'a') as f:
                for unique in uniques:
                    wrapped_seq = "\n".join(tw.wrap(unique[1]))
                    fasta_entry = '>'+unique[0]+'\n'+wrapped_seq+'\n\n'
                    f.write(fasta_entry)

    subprocess.run([makeblastdb_exe, '-in',
                   largepdb_collected_fasta,
                   '-dbtype', 'prot',
                   '-out', os.path.join(seqdir, 'largedb')])

    subprocess.run([makeblastdb_exe, '-in',
                   mono_collected_fasta,
                   '-dbtype', 'prot',
                   '-out', os.path.join(seqdir, 'monodb')])

    subprocess.run([makeblastdb_exe, '-in',
                   hetero_collected_fasta,
                   '-dbtype', 'prot',
                   '-out', os.path.join(seqdir, 'heterodb')])

    subprocess.run([makeblastdb_exe, '-in',
                   homo_collected_fasta,
                   '-dbtype', 'prot',
                   '-out', os.path.join(seqdir, 'homodb')])


def update_gesamt(verbosity):
    '''
    Runs gesamt to update its archive.
    Called by: update_databases()
    '''
    if verbosity == 1:
        vflag = '-v2'
    else:
        vflag = ''
    if not os.listdir(ges_homo_archive):
        subprocess.run([gesamt_exe,
                        '--make-archive',
                        ges_homo_archive,
                        '-pdb',
                        pdb_homo_archive,
                        vflag])
    else:
        subprocess.run([gesamt_exe,
                        '--update-archive',
                        ges_homo_archive,
                        '-pdb',
                        pdb_homo_archive,
                        vflag])


def update_databases(verbosity):
    '''
    Calls the four update functions sequentially (pdb, pdb1, ProtCHOIR, GESAMT)
    Called by: RunProtCHOIR.py:main()
    '''
    print('\n\nWill now update ALL databases.\n')
    create_directories()
    print('Updating PDB database...')
    update_pdb(verbosity)
    print('\n\nDone updating pdb!\n')
    print('Updating PDB1 database...')
    update_pdb1(verbosity)
    print('\n\nDone updating pdb1!\n')
    print('Updating UniRef50 database...')
    update_uniref(verbosity)
    print('\n\nDone updating UniRef50!\n')
    print('Updating seqres database...')
    update_seqres(verbosity)
    print('\n\nDone updating seqres!\n')
    print('Curating ProtCHOIR database...')
    curate_homoDB(verbosity)
    print('\n\nDone Curating ProtCHOIR database!\n')
    print('Collecting fasta files...')
    collect_fasta(verbosity)
    print('\n\nDone collecting fasta files!\n')
    print('\n\nUpdating GESAMT archive...\n')
    update_gesamt(verbosity)
    print('\n\nDone updating GESAMT archive!\n')
