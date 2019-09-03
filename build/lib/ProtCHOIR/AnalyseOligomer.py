# Imports
###############################################################################
import os
import pickle
import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ProtCHOIR.Initialise import *
from matplotlib.lines import Line2D
import ProtCHOIR.Toolbox as pctools
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


# Classes
###############################################################################
class ComplexRadar():
    def __init__(self, fig, variables, ranges,
                 n_ordinate_levels=8):
        angles = np.arange(0, 360, 360./len(variables))
        axes = [fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True,
                label="axes{}".format(i))
                for i in range(len(variables))]
        l, text = axes[0].set_thetagrids(angles,
                                         labels=variables)
        labels = []

        for txt, angle in zip(text, angles):
            lab = axes[0].text(np.deg2rad(angle), -0.12, txt.get_text(),
                               transform=txt.get_transform(), ha=txt.get_ha(),
                               va=txt.get_va())

            lab.set_rotation(angle-90)
            lab.set_fontweight('bold')
            lab.set_fontsize(50)
            labels.append(lab)
            axes[0].set_xticklabels([])

        for ax in axes[1:]:
            ax.patch.set_visible(False)
            ax.grid(True)
            ax.xaxis.set_visible(False)

        for i, ax in enumerate(axes):
            grid = np.linspace(*ranges[i], num=n_ordinate_levels)
            gridlabel = ["{}".format(abs(round(x, 1)))
                         for x in grid]
            gridlabel[0] = ""  # clean up origin
            ax.set_rgrids(grid, labels=gridlabel, angle=angles[i])
            ax.set_ylim(*ranges[i])

        # variables for plotting
        self.angle = np.deg2rad(np.r_[angles, angles[0]])
        self.ranges = ranges
        self.ax = axes[0]

    def plot(self, data, *args, **kw):
        sdata = _scale_data(data, self.ranges)
        self.ax.plot(self.angle, np.r_[sdata, sdata[0]], *args, **kw)

    def fill(self, data, *args, **kw):
        sdata = _scale_data(data, self.ranges)
        self.ax.fill(self.angle, np.r_[sdata, sdata[0]], *args, **kw)


# Functions
###############################################################################
def _invert(x, limits):
    """inverts a value x on a scale from
    limits[0] to limits[1]"""
    return limits[1] - (x - limits[0])


def _scale_data(data, ranges):
    """scales data[1:] to ranges[0],
    inverts if the scale is reversed"""
    for d, (y1, y2) in zip(data[1:], ranges[1:]):
        assert (y1 <= d <= y2) or (y2 <= d <= y1)
    x1, x2 = ranges[0]
    d = data[0]
    if x1 > x2:
        d = _invert(d, (x1, x2))
        x1, x2 = x2, x1
    sdata = [d]
    for d, (y1, y2) in zip(data[1:], ranges[1:]):
        if y1 > y2:
            d = _invert(d, (y1, y2))
            y1, y2 = y2, y1
        sdata.append((d-y1) / (y2-y1) * (x2 - x1) + x1)
    return sdata


def plot_deltas(model_name, template_name, interfaces_comparison, args):
    # Reset Matplotlib parameters to default
    matplotlib.rcParams.update(matplotlib.rcParamsDefault)
    # Plot everything
    interfaces = []
    model_areas = []
    model_energies = []
    model_hbonds = []
    model_saltbridges = []
    model_dissulfides = []
    template_areas = []
    template_energies = []
    template_hbonds = []
    template_saltbridges = []
    template_dissulfides = []
    n = 1
    for interface, comparison_data in interfaces_comparison.items():
        interfaces.append(interface)
        model_areas.append(comparison_data['model area'])
        model_energies.append(comparison_data['model energy'])
        model_hbonds.append(comparison_data['model hb'])
        model_saltbridges.append(comparison_data['model sb'])
        model_dissulfides.append(comparison_data['model ss'])
        template_areas.append(comparison_data['template area'])
        template_energies.append(comparison_data['template energy'])
        template_hbonds.append(comparison_data['template hb'])
        template_saltbridges.append(comparison_data['template sb'])
        template_dissulfides.append(comparison_data['template ss'])
        n += 1
    p, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(18, 24), sharex=True, gridspec_kw={'height_ratios': [1]*4})

    if n == 1:
        plt.suptitle('Interface comparison between '+model_name+' and '+template_name, fontsize=32, fontweight='bold')
    else:
        plt.suptitle('Interfaces comparison between '+model_name+' and '+template_name, fontsize=32, fontweight='bold')

    width = 18/(4*n)
    # Interpolate axis 0-18
    indices = []
    for i in range(1, n):
        a = float(i) / n
        x = (1 - a) * 0 + a * 18
        indices.append(x-width/2)

    indices = np.array(indices)

    # Plot Areas
    ax1.set_ylabel("Interface Area (A"+r"$^2$)", fontsize=24)
    ax1.grid(True, linestyle=':', linewidth=0.9, zorder=0, color='k')
    ax1.bar(indices, template_areas, width=width, label='Template')
    ax1.bar(indices + width, model_areas, width, label='Model')
    plt.xticks(indices + width/2, interfaces)
    ax1.tick_params(labelsize=22)
    ax1.set_xlim(0, 18)
    ax1.legend(loc='best', fontsize=24)

    # Plot Energies
    ax2.set_ylabel("Energy (kcal/mol)", fontsize=24)
    ax2.grid(True, linestyle=':', linewidth=0.9, zorder=0, color='k')
    ax2.bar(indices, template_energies, width, label='Template')
    ax2.bar(indices + width, model_energies, width, label='Model')
    ax2.tick_params(labelsize=22)

    # Plot Hydrogen Bonds
    ax3.set_ylabel("Hydrogen Bonds", fontsize=24)
    ax3.grid(True, linestyle=':', linewidth=0.9, zorder=0, color='k')
    ax3.bar(indices, template_hbonds, width, label='Template')
    ax3.bar(indices + width, model_hbonds, width, label='Model')
    ax3.tick_params(labelsize=22)

    # Plot Salt Bridges
    ax4.set_ylabel("Salt Bridges", fontsize=24)
    ax4.grid(True, linestyle=':', linewidth=0.9, zorder=0, color='k')
    ax4.bar(indices, template_saltbridges, width, label='Template')
    ax4.bar(indices + width, model_saltbridges, width, label='Model')
    ax4.tick_params(labelsize=22)

    # # Plot Disulfides
    # ax5.set_ylabel("Disulfide Bonds", fontsize=20)
    # ax5.grid(True, linestyle=':', linewidth=0.7, zorder=0, color='k')
    # ax5.bar(indices, template_dissulfides, width, label='Template')
    # ax5.bar(indices + width, model_dissulfides, width, label='Model')
    # ax5.tick_params(labelsize=18)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    p.subplots_adjust(hspace=0.07)
    outfile = model_name+'_'+template_name+'_CHOIR_InterfacesPlots.png'
    plt.savefig(outfile, dpi=300)
    # Close figure
    plt.close()
    print('Analysis plots for interfaces comparison generated : '+clrs['g']+os.path.basename(outfile)+clrs['n']+'\n')

    return outfile


def plot_molprobity(model_name, model_molprobity, template_name, template_molprobity):
    # Reset Matplotlib parameters to default
    matplotlib.rcParams.update(matplotlib.rcParamsDefault)
    ramaout_max = max([template_molprobity['rama_out'], model_molprobity['rama_out']])
    ramafav_max = max([template_molprobity['rama_fav'], model_molprobity['rama_fav']])
    rotout_max = max([template_molprobity['rot_out'], model_molprobity['rot_out']])
    clashscore_max = max([template_molprobity['clashscore'], model_molprobity['clashscore']])
    cbdev_max = max([template_molprobity['cb_dev'], model_molprobity['cb_dev']])
    molprobity_max = max([template_molprobity['molprobity_score'], model_molprobity['molprobity_score']])

    metrics = ('Rama. Out', 'Rama. Fav.', 'Rot. out.', 'CB Dev.', 'Clashscore', 'Molprobity')

    template_data = (-template_molprobity['rama_out'],
                     template_molprobity['rama_fav'],
                     -template_molprobity['rot_out'],
                     -template_molprobity['cb_dev'],
                     -template_molprobity['clashscore'],
                     -template_molprobity['molprobity_score'])

    model_data = (-model_molprobity['rama_out'],
                  model_molprobity['rama_fav'],
                  -model_molprobity['rot_out'],
                  -model_molprobity['cb_dev'],
                  -model_molprobity['clashscore'],
                  -model_molprobity['molprobity_score'])

    ranges = [(-100, 0),
              (0, 100),
              (-100, 0),
              (-100, 0),
              (-2*clashscore_max, 0),
              (-2*molprobity_max, 0)]

    angles = np.arange(0, 360, 360./len(ranges))
    fig = plt.figure(figsize=(20, 20))
    plt.rcParams.update({'font.size': 35})
    radar = ComplexRadar(fig, metrics, ranges)
    radar.plot(template_data, color='olive')
    radar.fill(template_data, alpha=0.2, color='olive')
    radar.plot(model_data, color='orange')
    radar.fill(model_data, alpha=0.2, color='orange')
    plt.rcParams.update({'font.size': 45})

    # Draw legend
    legend_elements = [Line2D([0], [0], color='olive', lw=30, label='Template'),
                       Line2D([0], [0], color='orange', lw=30, label='Model')]
    fig.legend(handles=legend_elements, ncol=1, loc='upper left', frameon=False)
    outfile = model_name+'_'+template_name+'_CHOIR_RadarPlots.png'
    plt.savefig(outfile, dpi=300)
    # Close figure
    plt.close()
    print('Molprobity radar plot generated: '+clrs['g']+os.path.basename(outfile)+clrs['n']+'\n')

    return './'+outfile


def analyse_oligomers(input_file, template_hitchain, oligomers_list, interfaces_dict, report, args, entropies=None, z_entropies=None, minx=None, maxx=None):
    pctools.print_section(3, 'OLIGOMER ANALYSIS')
    # Define template for comparisons
    template = template_hitchain.split(':')[0]
    template_file = template+'_CHOIR_RelevantChains.pdb'
    reports = []
    if 'M' in args.assessment:
        template_molprobity = pctools.run_molprobity(template_file, args)
    n = 0
    for oligomer in oligomers_list:
        n += 1
        report['model_filename'] = oligomer
        model_oligomer_name = os.path.basename(oligomer).split("_CHOIR_")[0].replace('.', '_')
        pctools.print_subsection('3.'+str(n), model_oligomer_name)
        print('Analysing oligomer file: '+clrs['y']+oligomer+clrs['n']+'\n')
        report['model_oligomer_name'] = model_oligomer_name
        report['model_figures'] = pctools.pymol_screenshot(oligomer, args, putty=True)

        if 'I' in args.assessment:
            pctools.print_subsection('3.'+str(n)+'[I]', 'Interfaces Comparison: '+model_oligomer_name)
            pdb_name, structure, nchains = pctools.parse_any_structure(oligomer)
            nchains, seqs, chain_ids = pctools.extract_seqs(structure, 0)
            relevant_chains = []
            for seq in seqs:
                relevant_chains.append(seq[0])

            pisa_output, pisa_error, protomer_data = pctools.run_pisa(oligomer, '', args.verbosity, gen_monomer_data=True, gen_oligomer_data=True)
            protomer_surface_residues = pctools.get_areas(protomer_data)
            report['assemblied_protomer_plot'], report['assemblied_protomer_exposed_area'], report['assemblied_protomer_hydrophobic_area'], report['assemblied_protomer_conserved_area'], minx, maxx = pctools.plot_analysis(pdb_name, protomer_surface_residues, entropies, z_entropies, args, minx=minx, maxx=maxx)

            if args.sequence_mode is False:
                report['exposed_area_reduction'] = round(100 * (float(report['protomer_exposed_area']) - float(report['assemblied_protomer_exposed_area'])) / float(report['protomer_exposed_area']), 1)
                report['hydrophobic_area_reduction'] = round(100 * (float(report['protomer_hydrophobic_area']) - float(report['assemblied_protomer_hydrophobic_area'])) / float(report['protomer_hydrophobic_area']), 1)
                report['conserved_area_reduction'] = round(100 * (float(report['protomer_conserved_area']) - float(report['assemblied_protomer_conserved_area'])) / float(report['protomer_conserved_area']), 1)

            model_oligomer = oligomer.split('_CHOIR_CorrectedChains')[0]
            xml_out = model_oligomer+'_CHOIR_PisaInterfaces.xml'
            model_interfaces_list, interfaces_output = pctools.parse_interfaces(xml_out, relevant_chains, args.verbosity)
            template_interfaces_list = interfaces_dict[template_hitchain]

            if args.verbosity > 0:
                print(clrs['y']+'MODEL INTERFACES'+clrs['n'])
                for model_interface in model_interfaces_list:
                    print(clrs['y']+' <> '.join(model_interface['chains'])+clrs['n'])
                    print(clrs['y']+'Interface Area: '+clrs['n']+str(model_interface['interface area'])+' A^2')
                    print(clrs['y']+'Interface Solvation Energy: '+clrs['n']+str(model_interface['interface solvation energy'])+' kcal/mol')
                    print(clrs['y']+'Hydrogen Bonds: '+clrs['n']+str(model_interface['hydrogen bonds']))
                    print(clrs['y']+'Salt Bridges: '+clrs['n']+str(model_interface['salt bridges']))
                    print(clrs['y']+'Disulphide Bridges: '+clrs['n']+str(model_interface['disulphide bridges'])+"\n\n")

            interfaces_comparison = {}
            for template_interface in template_interfaces_list:
                for model_interface in model_interfaces_list:
                    if set(model_interface['chains']) == set(template_interface['chains']):
                        comparison_data = {}

                        delta_area = round(model_interface['interface area']-template_interface['interface area'], 2)
                        comparison_data['model area'] = model_interface['interface area']
                        comparison_data['template area'] = template_interface['interface area']
                        comparison_data['delta area'] = delta_area
                        delta_energy = round(model_interface['interface solvation energy']-template_interface['interface solvation energy'], 2)
                        comparison_data['model energy'] = model_interface['interface solvation energy']
                        comparison_data['template energy'] = template_interface['interface solvation energy']
                        comparison_data['delta energy'] = delta_energy
                        delta_hb = round(model_interface['hydrogen bonds']-template_interface['hydrogen bonds'], 2)
                        comparison_data['model hb'] = model_interface['hydrogen bonds']
                        comparison_data['template hb'] = template_interface['hydrogen bonds']
                        comparison_data['delta hb'] = delta_hb
                        delta_sb = round(model_interface['salt bridges']-template_interface['salt bridges'], 2)
                        comparison_data['model sb'] = model_interface['salt bridges']
                        comparison_data['template sb'] = template_interface['salt bridges']
                        comparison_data['delta sb'] = delta_sb
                        delta_ss = round(model_interface['disulphide bridges']-template_interface['disulphide bridges'], 2)
                        comparison_data['model ss'] = model_interface['disulphide bridges']
                        comparison_data['template ss'] = template_interface['disulphide bridges']
                        comparison_data['delta ss'] = delta_ss
                        interfaces_comparison[''.join(sorted(model_interface['chains']))] = comparison_data

                        print(clrs['y']+'INTERFACES COMPARISON'+clrs['n'])
                        print(' <> '.join(model_interface['chains']))
                        if delta_area >= 0:
                            emphasis_color = clrs['g']
                        else:
                            emphasis_color = clrs['r']
                        print('Delta Interface Area: '+emphasis_color+str(delta_area)+clrs['n']+' A^2')
                        if delta_energy <= 0:
                            emphasis_color = clrs['g']
                        else:
                            emphasis_color = clrs['r']
                        print('Delta Interface Solvation Energy: '+emphasis_color+str(delta_energy)+clrs['n']+' kcal/mol')
                        if delta_hb >= 0:
                            emphasis_color = clrs['g']
                        else:
                            emphasis_color = clrs['r']
                        print('Delta Hydrogen Bonds: '+emphasis_color+str(delta_hb)+clrs['n'])
                        if delta_sb >= 0:
                            emphasis_color = clrs['g']
                        else:
                            emphasis_color = clrs['r']
                        print('Delta Salt Bridges: '+emphasis_color+str(delta_sb)+clrs['n'])
                        if delta_ss >= 0:
                            emphasis_color = clrs['g']
                        else:
                            emphasis_color = clrs['r']
                        print('Delta Disulphide Bridges: '+emphasis_color+str(delta_ss)+clrs['n']+"\n")

            report['comparison_plots'] = os.path.basename(plot_deltas(model_oligomer_name, template, interfaces_comparison, args))

        if 'G' in args.assessment:
            pctools.print_subsection('3.'+str(n)+'[G]', 'GESAMT Comparison')
            qscore, rmsd, fasta_out = pctools.run_gesamt(template, template_file, model_oligomer_name, oligomer, None, args, delfasta=True)
            report['gesamt_qscore'] = str(qscore)
            report['gesamt_rmsd'] = str(rmsd)

        if 'M' in args.assessment:
            pctools.print_subsection('3.'+str(n)+'[M]', 'Molprobity Comparison')
            model_molprobity = pctools.run_molprobity(oligomer, args)
            report['model_clashscore'] = str(model_molprobity['clashscore'])
            report['model_molprobity'] = str(model_molprobity['molprobity_score'])
            print(clrs['y']+'MOLPROBITY COMPARISON'+clrs['n'])
            print('Criterion\tTempl.\tModel')
            print('Rama. Fav.\t'+str(template_molprobity['rama_fav'])+'\t'+str(model_molprobity['rama_fav']))
            print('Rama. Out.\t'+str(template_molprobity['rama_out'])+'\t'+str(model_molprobity['rama_out']))
            print('Rot. Out.\t'+str(template_molprobity['rot_out'])+'\t'+str(model_molprobity['rot_out']))
            print('CBeta Dev.\t'+str(template_molprobity['cb_dev'])+'\t'+str(model_molprobity['cb_dev']))
            print('Clashscore\t'+str(template_molprobity['clashscore'])+'\t'+str(model_molprobity['clashscore']))
            print('Molprob. Score\t'+str(template_molprobity['molprobity_score'])+'\t'+str(model_molprobity['molprobity_score']))
            report['molprobity_radar'] = plot_molprobity(model_oligomer_name, model_molprobity, template, template_molprobity)

        pickle.dump(report, open(model_oligomer_name+'_CHOIR_Report.pickle', 'wb'))

        reports.append(report.copy())

    return reports
