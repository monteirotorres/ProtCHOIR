#!/usr/bin/env python
#
# Imports
###############################################################################
from datetime import datetime
# Mark Initiation time
start_timestamp = datetime.timestamp(datetime.now())-1
start_time = datetime.now()
import os
import sys
import time
import pickle
import shutil
import zipfile
import operator
import argparse
import textwrap as tw
import ProtCHOIR.Toolbox as pctools
from ProtCHOIR.Initialise import *
from multiprocessing import cpu_count
from ProtCHOIR.UpdateDatabases import update_databases
from ProtCHOIR.AnalyseProtomer import analyze_protomer
from ProtCHOIR.MakeOligomer import make_oligomer
from ProtCHOIR.AnalyseOligomer import analyse_oligomers

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

# Dictionaries
###############################################################################

# Global Variables
###############################################################################

# Classes
###############################################################################


# Functions
###############################################################################
def finalize(reports, input_basename, start_time, start_timestamp, args):
    report_data = ['input_filename', 'sequence_mode', 'templatedmodel', 'protomer_residues', 'tmspans', 'highest_scoring_state', 'homo_oligomeric_over_other_score', 'best_template', 'best_nchains', 'best_id', 'best_cov', 'best_qscore', 'model_oligomer_name', 'model_molprobity', 'gesamt_rmsd', 'quality_score', 'surface_score', 'interfaces_score', 'protchoir_score', 'total_runtime', 'exit']
    if type(reports) is list:
        if args.zip_output == 2:
            # Don't prevent compression of anything
            nozip = []
            for report in reports:
                if args.generate_report is True:
                    report['html_report'] = pctools.html_report(report, args)
        else:
            # Prevent compression of files needed for the report and the models
            nozip = [os.path.basename(report['model_filename']) for report in reports]
            for report in reports:
                if args.generate_report is True:
                    report['html_report'] = pctools.html_report(report, args)
                    for key, value in report.items():
                        if key in ['html_report', 'molprobity_radar', 'comparison_plots', 'protomer_figure', 'protomer_plot', 'template_figure', 'topology_figure', 'assemblied_protomer_plot', 'input_filename']:
                            nozip.append(os.path.basename(value))
                        if key == 'model_figures':
                            for figure in value:
                                nozip.append(os.path.basename(figure))


        best_report = sorted(reports, key=operator.itemgetter('protchoir_score'))[-1]

    elif type(reports) is dict:
        nozip = []
        best_report = reports
        for data in report_data:
            if data not in best_report:
                best_report[data] = 'NA'

    # Generate summary tsv file for the best report
    end_time = datetime.now()
    runtime = end_time - start_time
    best_report['total_runtime'] = str(runtime.seconds)
    summary_file = input_basename+'_CHOIR_Summary.tsv'
    nozip.append(summary_file)
    if 'exit' not in best_report:
        best_report['exit'] = '0'
        with open('CHOIR_Progress.out', 'a') as f:
            f.write(datetime.now().strftime("%H:%M:%S")+": Finished running ProtCHOIR!")
    elif best_report['exit'] == '1':
        with open('CHOIR_Progress.out', 'a') as f:
            f.write(datetime.now().strftime("%H:%M:%S")+": ERROR! Indicated template not found in oligomers database...")
    elif best_report['exit'] == '2':
        with open('CHOIR_Progress.out', 'a') as f:
            f.write(datetime.now().strftime("%H:%M:%S")+": ERROR! Failed to find suitable homologues...")
    elif best_report['exit'] == '3':
        with open('CHOIR_Progress.out', 'a') as f:
            f.write(datetime.now().strftime("%H:%M:%S")+": ERROR! Failed to find suitable homo-oligomeri interfaces...")
    elif best_report['exit'] == '4':
        with open('CHOIR_Progress.out', 'a') as f:
            f.write(datetime.now().strftime("%H:%M:%S")+": ERROR! No template had an average Q-score above cut-off...")
    elif best_report['exit'] == '5':
        with open('CHOIR_Progress.out', 'a') as f:
            f.write(datetime.now().strftime("%H:%M:%S")+": ERROR! Failed to find templates in local databases...")
    elif best_report['exit'] == '6':
        with open('CHOIR_Progress.out', 'a') as f:
            f.write(datetime.now().strftime("%H:%M:%S")+": ERROR! Sub-optimal alignment between template and target sequences...")

    with open(summary_file, 'w') as f:
        f.write('Input\tSeq.Mode\tTemplated\tLength\tTMSpans\tLikelyState\tH3OScore\tTemplate\tChains\tIdentity\tCoverage\tAv.QScore\tBestModel\tMolprobity\tRMSD\tQuality\tSurface\tInterfaces\tProtCHOIR\tRuntime\tExit\n')
        f.write('\t'.join([str(best_report[data]) for data in report_data])+'\n')
    # Finalise
    final_end_time = datetime.timestamp(datetime.now())
    time.sleep(1)

    # Compress output
    if args.zip_output > 0:
        try:
            import zlib
            compression = zipfile.ZIP_DEFLATED
        except (ImportError, AttributeError):
            compression = zipfile.ZIP_STORED

        with zipfile.ZipFile(input_basename+'_ProtCHOIR_OUT.zip', 'w', compression=compression) as zipf:
            for f in os.listdir(os.getcwd()):
                if f != input_basename+'_ProtCHOIR_OUT.zip' and os.path.getctime(f) > start_timestamp and os.path.getctime(f) < final_end_time:
                    print('Compressing... '+f)
                    zipf.write(f)
                    if f not in nozip:
                        if os.path.isdir(f):
                            shutil.rmtree(f)
                        elif os.path.isfile(f):
                            os.remove(f)

    print('FINISHED AT: '+datetime.now().strftime("%d-%m-%Y %H:%M"))
    print('TOTAL RUNTIME: '+str(runtime.seconds)+' s')

# Main Function
###############################################################################
def main():

    args = initial_args

    # Define multiprocessing options
    args.available_cores = cpu_count()

    if args.force_single_core is True:
        args.multiprocess = False
        args.psiblast_threads = 1
        args.modeller_threads = 1
    else:
        if args.psiblast_threads is None:
            args.psiblast_threads = args.available_cores
        if args.modeller_threads is None:
            args.modeller_threads = min([args.available_cores, args.models])

    if args.update is True:
        print(tw.dedent("""
                                         !WARNING!

                      You have chosen to updtate the local databases.

              ** The root directory for the database files is: """+clrs['y']+choirdb+clrs['n']+"""

              ** The path to local pdb mirror is: """+clrs['y']+pdb_archive+clrs['n']+"""

              ** The path to local pdb biounit mirror is: """+clrs['y']+pdb1_archive+clrs['n']+"""

              ** The path to local gesamt archive is: """+clrs['y']+ges_homo_archive+clrs['n']+"""

              ** The path to local UniRef50 blast database is: """+clrs['y']+uniref50+clrs['n']+"""


              This could take a long time.

              <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

              """))
        option = input('Do you confirm the information above? (y/n)')
        if option == 'y' or option == 'Y' or option == 'YES' or option == 'yes' or option == 'Yes':
            update_databases(args.verbosity)
            print('\n\nDone updating all databases. Exiting.\n')
        else:
            print('\n\nNo positive confirmation, will not update databases.\n')
            exit()
    # Actually run oligomerization protocol
    else:
        outdir = os.getcwd()
        input_file = args.input_file
        assert os.path.isdir(pdb_archive), clrs['r']+'\n\n Not able to find PDB directory.\n\n Does "'+pdb_archive+'" exist?'+clrs['n']
        assert os.path.isdir(pdb1_archive), clrs['r']+'\n\n Not able to find PDB1 assemblies directory.\n\n Does "'+pdb1_archive+'" exist?'+clrs['n']
        assert os.path.isdir(pdb_homo_archive), clrs['r']+'\n\n Not able to find ProtCHOIR database directory.\n\n Does "'+pdb_homo_archive+'" exist?'+clrs['n']
        assert os.path.isdir(ges_homo_archive), clrs['r']+'\n\n Not able to find GESAMT archive directory.\n\n Does "'+ges_homo_archive+'" exist?'+clrs['n']
        assert args.refine_level in [0, 1, 2, 3, 4], clrs['r']+'\n\n Refinement level must be an integer number from 0 to 4.\n Run ProtCHOIR -h for more information\n\n'+clrs['n']
        assert args.psiblast_params in psiblast_params, clrs['r']+'\n\n PSI-BLAST parameters invalid.\n Run ProtCHOIR -h for more information\n\n'+clrs['n']
        assert input_file is not None, clrs['r']+'\n\n Please inform the input file name.\n Run ProtCHOIR -h for more information.\n\n'+clrs['n']
        assert os.path.isfile(input_file), clrs['r']+'\n\n Not able to find input file.\n\n Does "'+input_file+'" exist?\n'+clrs['n']
        assert args.zip_output in [0, 1, 2], clrs['r']+'\n\n Compression level must be an integer number between 0 and 2.\n Run ProtCHOIR -h for more information\n\n'+clrs['n']
        assert all([i in set('MIG') for i in set(args.assessment)]) or args.assessment == 'N', clrs['r']+'\n\n Oligomer assessment type do not comply.\n Choose any combination of [G]Gesamt, [M]Molprobity, [I]Interfaces or choose [N] for None\n\n'+clrs['n']

        # Force generation of topologies and all assessments if final report is requested
        if args.generate_report is True:
            args.assessment = 'MIG'
            args.plot_topologies = True

        # Deal with dots and dashes in the input file and remove dots
        if input_file.lower().endswith('.pdb'):
            input_basename = os.path.basename(input_file).split('.pdb')[0]
            input_basename = input_basename.replace(".", "_")
            input_basename = input_basename.replace("-", "_")
            new_input_file = input_basename+'.pdb'
            if os.path.basename(input_file) == os.path.basename(new_input_file):
                pass
            else:
                shutil.copy(input_file, new_input_file)

        # Also process filename to fasta header if input file is fasta
        elif input_file.lower().endswith('.fasta'):
            input_basename = os.path.basename(input_file).split('.fasta')[0]
            input_basename = input_basename.replace(".", "_")
            input_basename = input_basename.replace("-", "_")
            new_input_file = os.path.join(outdir, input_basename+'_CHOIR_MonomerSequence.fasta')
            with open(input_file, 'r') as infile, open(new_input_file, 'w') as outfile:
                outfile.write('>'+input_basename+'\n')
                n = 0
                for line in infile.readlines():
                    if not line.startswith('>'):
                        outfile.write(line)
                    else:
                        n += 1
                    if n == 2:
                        break
            args.sequence_mode = True
        else:
            raise pctools.FileFormatError(clrs['r']+'\n\n Input format must be either pdb or fasta\n Run ./ProtCHOIR -h for more information\n\n'+clrs['n'])
        if args.allow_monomers:
            assert args.sequence_mode is True, clrs['r']+'\n\n To allow building monomers you must use sequence mode. \n Run ProtCHOIR -h for more information\n\n'+clrs['n']

        # Start recording job progress
        with open('CHOIR_Progress.out', 'w') as f:
            f.write("Starting new ProtCHOIR run\n")

        # Pickle Runtime arguments
        pickle.dump(args, open('CHOIR_Args.pickle', 'wb'))

        # Show arguments used and create CHOIR.conf
        pctools.print_section(0, "Runtime Arguments")
        runtime_arguments = {}
        choir_args = os.path.join(outdir, "CHOIR.args")
        with open(choir_args, 'w') as f:
            for name, value in vars(args).items():
                runtime_arguments[name] = value
                print(name+"="+str(value))
                f.write(name+"="+str(value)+"\n")
        print('\nRuntime parameters written to: '+clrs['g']+os.path.basename(choir_args)+clrs['n']+'\n')

        # Initialize report
        report = {}
        report['runtime_arguments'] = runtime_arguments
        report['input_filename'] = os.path.basename(new_input_file)

        # Write errorprof placeholder summary
        placeholder_report = report.copy()
        report_data = ['input_filename', 'sequence_mode', 'templatedmodel', 'protomer_residues', 'tmspans', 'highest_scoring_state', 'homo_oligomeric_over_other_score', 'best_template', 'best_nchains', 'best_id', 'best_cov', 'best_qscore', 'model_oligomer_name', 'model_molprobity', 'gesamt_rmsd', 'protchoir_score', 'surface_score', 'interfaces_score', 'quality_score', 'total_runtime', 'exit']
        for data in report_data:
            if data not in placeholder_report:
                placeholder_report[data] = 'NA'
        with open(input_basename+'_CHOIR_Summary.tsv', 'w') as f:
            f.write('Input\tSeq.Mode\tTemplated\tLength\tTMSpans\tLikelyState\tH3OScore\tTemplate\tChains\tIdentity\tCoverage\tAv.QScore\tBestModel\tMolprobity\tRMSD\tProtCHOIR\tSurface\tInterfaces\tQuality\tRuntime\tExit\n')
            f.write('\t'.join([str(placeholder_report[data]) for data in report_data])+'\n')

        # Start analysis of protomer
        analyse_protomer_results, report, args = analyze_protomer(new_input_file, report, args)

        # If no suitable homo-oligomeric template wasfound, exit nicely.
        if analyse_protomer_results is None:
            finalize(report, input_basename, start_time, start_timestamp, args)
            pctools.print_sorry()
            sys.exit(0)

        # Else, proceed conditionally on runtime arguments
        elif analyse_protomer_results is not None and args.sequence_mode is True:
            residue_index_mapping = None
            minx = None
            maxx = None
            if args.skip_conservation:
                entropies = None
                z_entropies = None
                pdb_name, clean_input_file, largest_oligo_complexes, interfaces_dict, tmdata = analyse_protomer_results
            elif not args.skip_conservation:
                pdb_name, clean_input_file, largest_oligo_complexes, interfaces_dict, entropies, z_entropies, tmdata = analyse_protomer_results
                if entropies == z_entropies == minx == maxx == None:
                    args.skip_conservation = True


        elif analyse_protomer_results is not None and args.sequence_mode is False:
            if args.skip_conservation:
                minx = None
                maxx = None
                entropies = None
                z_entropies = None
                pdb_name, clean_input_file, largest_oligo_complexes, interfaces_dict, residue_index_mapping, tmdata = analyse_protomer_results
            elif not args.skip_conservation:
                pdb_name, clean_input_file, largest_oligo_complexes, interfaces_dict, entropies, z_entropies, residue_index_mapping, minx, maxx, tmdata = analyse_protomer_results
                if entropies == z_entropies == minx == maxx == None:
                    args.skip_conservation = True

        report['runtime_arguments']['skip_conservation'] = args.skip_conservation

        new_input_file = clean_input_file

        # Use information of complexes to build oligomers
        best_oligo_template, built_oligomers, report = make_oligomer(new_input_file, largest_oligo_complexes, report, args, residue_index_mapping=residue_index_mapping)

        # If no models were built, exit nicely.
        if built_oligomers is None:
            finalize(report, input_basename, start_time, start_timestamp, args)
            pctools.print_sorry()
            sys.exit(0)

        # Analyse built models
        reports = analyse_oligomers(new_input_file, best_oligo_template, built_oligomers, interfaces_dict, tmdata, report, args, entropies=entropies, z_entropies=z_entropies, minx=minx, maxx=maxx)
        finalize(reports, input_basename, start_time, start_timestamp, args)



# Execute
###############################################################################
if __name__ == "__main__":
    main()
