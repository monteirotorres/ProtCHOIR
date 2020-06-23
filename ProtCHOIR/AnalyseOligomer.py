# Imports
###############################################################################
import os
import math
import pickle
import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool
from ProtCHOIR.Initialise import *
from matplotlib.lines import Line2D
import ProtCHOIR.Toolbox as pctools
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
        if y1 == y2:
            sdata.append(y1)
        else:
            sdata.append((d-y1) / (y2-y1) * (x2 - x1) + x1)
    return sdata


def plot_deltas(model_name, template_name, interfaces_comparison, args):
    output = []
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
    p, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, figsize=(18, 24), sharex=True, gridspec_kw={'height_ratios': [1]*5})

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

    # Plot Disulfides
    ax5.set_ylabel("Disulfide Bonds", fontsize=24)
    ax5.grid(True, linestyle=':', linewidth=0.9, zorder=0, color='k')
    ax5.bar(indices, template_dissulfides, width, label='Template')
    ax5.bar(indices + width, model_dissulfides, width, label='Model')
    ax5.tick_params(labelsize=22)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    p.subplots_adjust(hspace=0.07)
    outfile = model_name+'_'+template_name+'_CHOIR_InterfacesPlots.png'
    plt.savefig(outfile, dpi=300)
    # Close figure
    plt.close()
    output = 'Analysis plots for interfaces comparison generated : '+clrs['g']+os.path.basename(outfile)+clrs['n']+'\n'

    return outfile, output


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
              (-2*cbdev_max, 0),
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
    output = 'Molprobity radar plot generated: '+clrs['g']+os.path.basename(outfile)+clrs['n']+'\n'

    return './'+outfile, output


def analyse_model(oligomer):
    output = []
    model_report = g_report.copy()
    model_report['model_filename'] = oligomer
    model_oligomer_name = os.path.basename(oligomer).split("_CHOIR_")[0].replace('.', '_')
    output.append(pctools.subsection('3', model_oligomer_name))
    output.append('Analysing oligomer file: '+clrs['y']+oligomer+clrs['n']+'\n')
    model_report['model_oligomer_name'] = model_oligomer_name
    if g_args.generate_report is True:
        model_report['model_figures'], pymol_output = pctools.pymol_screenshot(oligomer, g_args, putty=True)
        output.append(pymol_output)
    pdb_name, structure, nchains = pctools.parse_any_structure(oligomer)
    nchains, seqs, chain_ids = pctools.extract_seqs(structure, 0)
    relevant_chains = []
    for seq in seqs:
        relevant_chains.append(seq[0])

    pisa_output, pisa_error, protomer_data = pctools.run_pisa(oligomer, '', g_args.verbosity, gen_monomer_data=True, gen_oligomer_data=True)
    protomer_surface_residues = pctools.get_areas(protomer_data)
    model_report['assemblied_protomer_plot'], model_report['assemblied_protomer_exposed_area'], model_report['assemblied_protomer_hydrophobic_area'], model_report['assemblied_protomer_conserved_area'], minx, maxx, analysis_output = pctools.plot_analysis(pdb_name, protomer_surface_residues, g_entropies, g_z_entropies, g_tmdata, g_args, minx=g_minx, maxx=g_maxx)
    output.append(analysis_output)

    if 'I' in g_args.assessment and not g_args.allow_monomers:
        output.append(pctools.subsection('3'+'[I]', 'Interfaces Comparison: '+model_oligomer_name))
        if g_args.sequence_mode is False and g_args.skip_conservation is False:
            model_report['exposed_area_reduction'] = int(100 * (float(model_report['assemblied_protomer_exposed_area']) - float(model_report['protomer_exposed_area'])) / float(model_report['protomer_exposed_area']))
            model_report['hydrophobic_area_reduction'] = int(100 * (float(model_report['assemblied_protomer_hydrophobic_area']) - float(model_report['protomer_hydrophobic_area'])) / float(model_report['protomer_hydrophobic_area']))
            model_report['conserved_area_reduction'] = int(100 * (float(model_report['assemblied_protomer_conserved_area']) - float(model_report['protomer_conserved_area'])) / float(model_report['protomer_conserved_area']))

            if model_report['exposed_area_reduction'] < -5:
                if model_report['hydrophobic_area_reduction'] < 0:
                    hydophobic_surface_score = 10*(model_report['hydrophobic_area_reduction']/model_report['exposed_area_reduction'])/3
                else:
                    hydophobic_surface_score = 0
                if hydophobic_surface_score > 10:
                    hydophobic_surface_score = 10
                output.append('Hydrophobic surface score: '+str(hydophobic_surface_score))
                if model_report['conserved_area_reduction'] < 0:
                    conserved_surface_score = 10*(model_report['conserved_area_reduction']/model_report['exposed_area_reduction'])/3
                else:
                    conserved_surface_score = 0
                if conserved_surface_score > 10:
                    conserved_surface_score = 10
                output.append('Conserved surface score: '+str(conserved_surface_score))
                model_report['surface_score'] = round((hydophobic_surface_score+conserved_surface_score)/2, 2)
            else:
                output.append(clrs['r']+'Exposed area reduction too small.'+clrs['n'])
                model_report['surface_score'] = 0
            output.append('Final surface score: '+str(model_report['surface_score']))
        else:
            model_report['surface_score'] = 'NA'

        model_oligomer = oligomer.split('_CHOIR_CorrectedChains')[0]
        xml_out = model_oligomer+'_CHOIR_PisaInterfaces.xml'
        model_interfaces_list, interfaces_output = pctools.parse_interfaces(xml_out, relevant_chains, g_args.verbosity)
        template_interfaces_list = g_interfaces_dict[g_template_hitchain]

        if model_interfaces_list and template_interfaces_list:
            if g_args.verbosity > 0:
                output.append(clrs['y']+'MODEL INTERFACES'+clrs['n'])
                for model_interface in model_interfaces_list:
                    output.append(clrs['y']+' <> '.join(model_interface['chains'])+clrs['n'])
                    output.append(clrs['y']+'Interface Area: '+clrs['n']+str(model_interface['interface area'])+' A^2')
                    output.append(clrs['y']+'Interface Solvation Energy: '+clrs['n']+str(model_interface['interface solvation energy'])+' kcal/mol')
                    output.append(clrs['y']+'Hydrogen Bonds: '+clrs['n']+str(model_interface['hydrogen bonds']))
                    output.append(clrs['y']+'Salt Bridges: '+clrs['n']+str(model_interface['salt bridges']))
                    output.append(clrs['y']+'Disulphide Bridges: '+clrs['n']+str(model_interface['disulphide bridges'])+"\n\n")

            interfaces_comparison = {}
            for template_interface in template_interfaces_list:
                for model_interface in model_interfaces_list:
                    if set(model_interface['chains']) == set(template_interface['chains']):
                        comparison_data = {}
                        denominator = 12
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


                        output.append(clrs['y']+'INTERFACES COMPARISON'+clrs['n'])
                        output.append(' <> '.join(model_interface['chains']))
                        if delta_area >= 0:
                            emphasis_color = clrs['g']
                            relative_area = 100
                        else:
                            emphasis_color = clrs['r']
                            relative_area = round(model_interface['interface area'] * 100 / template_interface['interface area'], 2)
                        output.append('Delta Interface Area: '+emphasis_color+str(delta_area)+clrs['n']+' A^2 ('+str(relative_area)+'%)')

                        if delta_energy <= 0:
                            emphasis_color = clrs['g']
                            relative_energy = 100
                        else:
                            emphasis_color = clrs['r']
                            if model_interface['interface solvation energy'] < 0 and template_interface['interface solvation energy'] < 0:
                                relative_energy = round(model_interface['interface solvation energy'] * 100 / template_interface['interface solvation energy'], 2)
                            elif model_interface['interface solvation energy'] > 0 and template_interface['interface solvation energy'] < 0:
                                relative_energy = 0
                            elif model_interface['interface solvation energy'] < 0 and template_interface['interface solvation energy'] > 0:
                                relative_energy = 100
                            elif model_interface['interface solvation energy'] > 0 and template_interface['interface solvation energy'] > 0:
                                relative_energy = 0
                        output.append('Delta Interface Solvation Energy: '+emphasis_color+str(delta_energy)+clrs['n']+' kcal/mol ('+str(relative_energy)+'%)')

                        if model_interface['hydrogen bonds'] == template_interface['hydrogen bonds'] == 0:
                            relative_hb = 0
                            emphasis_color = clrs['r']
                            denominator -= 2
                        elif delta_hb >= 0:
                            relative_hb = 100
                            emphasis_color = clrs['g']
                        else:
                            emphasis_color = clrs['r']
                            relative_hb = round(model_interface['hydrogen bonds'] * 100 / template_interface['hydrogen bonds'], 2)
                        output.append('Delta Hydrogen Bonds: '+emphasis_color+str(delta_hb)+clrs['n']+' ('+str(relative_hb)+'%)')

                        if model_interface['salt bridges'] == template_interface['salt bridges'] == 0:
                            relative_sb = 0
                            emphasis_color = clrs['r']
                            denominator -= 3
                        elif delta_sb >= 0:
                            relative_sb = 100
                            emphasis_color = clrs['g']
                        else:
                            relative_sb = round(model_interface['salt bridges'] * 100 / template_interface['salt bridges'], 2)
                            emphasis_color = clrs['r']
                        output.append('Delta Salt Bridges: '+emphasis_color+str(delta_sb)+clrs['n']+' ('+str(relative_sb)+'%)')

                        if model_interface['disulphide bridges'] == template_interface['disulphide bridges'] == 0:
                            relative_ss = 0
                            emphasis_color = clrs['r']
                            denominator -= 4
                        elif delta_ss >= 0:
                            relative_ss = 100
                            emphasis_color = clrs['g']
                        else:
                            relative_ss = round(model_interface['disulphide bridges'] * 100 / template_interface['disulphide bridges'], 2)
                            emphasis_color = clrs['r']
                        output.append('Delta Disulphide Bridges: '+emphasis_color+str(delta_ss)+clrs['n']+' ('+str(relative_ss)+'%)\n')

                        if denominator == 0:
                                comparison_data['score'] = 0
                        else:
                            comparison_data['score'] = round((relative_area+2*relative_energy+2*relative_hb+3*relative_sb+4*relative_ss)/denominator, 2)
                        output.append('Interface score: '+str(comparison_data['score']))
                        interfaces_comparison[''.join(sorted(model_interface['chains']))] = comparison_data

            comparison_plots, interfaces_output = plot_deltas(model_oligomer_name, template, interfaces_comparison, g_args)
            model_report['comparison_plots'] = os.path.basename(comparison_plots)
            output.append(interfaces_output)
            summed_score = 0
            for interface, data in interfaces_comparison.items():
                summed_score += data['score']

            model_report['interfaces_score'] = round(summed_score/(10*len(interfaces_comparison)), 2)
            output.append('Final interfaces score: '+str(model_report['interfaces_score']))
        else:
            if 'surface_score' not in model_report:
                model_report['surface_score'] = 0
            model_report['interfaces_score'] = 0

    else:
        model_report['surface_score'] = 'NA'
        model_report['interfaces_score'] = 'NA'
        model_report['comparison_plots'] = 'NA'
        model_report['assemblied_protomer_exposed_area'] = 'NA'
        model_report['assemblied_protomer_hydrophobic_area'] = 'NA'
        model_report['assemblied_protomer_conserved_area'] = 'NA'

    if 'G' in g_args.assessment:
        output.append(pctools.subsection('3'+'[G]', 'GESAMT Comparison'))
        qscore, rmsd, fasta_out, gesamt_output = pctools.run_gesamt(template, template_file, model_oligomer_name, oligomer, None, g_args)
        output.append(gesamt_output)
        model_report['gesamt_qscore'] = str(qscore)
        model_report['gesamt_rmsd'] = str(rmsd)
    else:
        model_report['gesamt_qscore'] = 'NA'
        model_report['gesamt_rmsd'] = 'NA'


    if 'M' in g_args.assessment:
        output.append(pctools.subsection('3'+'[M]', 'Molprobity Comparison'))
        model_molprobity, molprobity_output = pctools.run_molprobity(oligomer, g_args)
        output.append(molprobity_output)
        model_report['model_clashscore'] = str(model_molprobity['clashscore'])
        model_report['model_molprobity'] = str(model_molprobity['molprobity_score'])
        output.append(clrs['y']+'MOLPROBITY COMPARISON'+clrs['n'])
        output.append('Criterion\tTempl.\tModel')
        output.append('Rama. Fav.\t'+str(template_molprobity['rama_fav'])+'\t'+str(model_molprobity['rama_fav']))
        output.append('Rama. Out.\t'+str(template_molprobity['rama_out'])+'\t'+str(model_molprobity['rama_out']))
        output.append('Rot. Out.\t'+str(template_molprobity['rot_out'])+'\t'+str(model_molprobity['rot_out']))
        output.append('CBeta Dev.\t'+str(template_molprobity['cb_dev'])+'\t'+str(model_molprobity['cb_dev']))
        output.append('Clashscore\t'+str(template_molprobity['clashscore'])+'\t'+str(model_molprobity['clashscore']))
        output.append('Molprob. Score\t'+str(template_molprobity['molprobity_score'])+'\t'+str(model_molprobity['molprobity_score']))
        molprobity_radar, radar_output = plot_molprobity(model_oligomer_name, model_molprobity, template, template_molprobity)
        output.append(radar_output)
        model_report['molprobity_radar'] = molprobity_radar
        delta_clashscore = (model_molprobity['clashscore'] - template_molprobity['clashscore'])/10
        output.append('Delta clashscore: '+str(delta_clashscore))
        if delta_clashscore >= 1:
            model_report['quality_score'] = round(10 - math.log(delta_clashscore**5, 10), 2)
        else:
            model_report['quality_score'] = 10
        output.append('Final quality score: '+str(model_report['quality_score']))
    else:
        model_report['model_clashscore'] = 'NA'
        model_report['model_molprobity'] = 'NA'
        model_report['quality_score'] = 'NA'

    if 'M' in g_args.assessment and 'I' in g_args.assessment and not g_args.allow_monomers:
        if g_args.sequence_mode is False and g_args.skip_conservation is False:
            model_report['protchoir_score'] = round(sum([model_report['interfaces_score'], model_report['surface_score'], model_report['quality_score']])/3, 2)
        else:
            model_report['protchoir_score'] = round(sum([model_report['interfaces_score'], model_report['quality_score']])/2, 2)
    elif 'M' in g_args.assessment:
        model_report['protchoir_score'] = model_report['quality_score']
    elif 'I' in g_args.assessment:
        if g_args.sequence_mode is False and g_args.skip_conservation is False:
            model_report['protchoir_score'] = round(sum([model_report['interfaces_score'], model_report['surface_score']])/2, 2)
        else:
            model_report['protchoir_score'] = model_report['interfaces_score']
    else:
        model_report['protchoir_score'] = 'NA'
    if str(model_report['protchoir_score']) == 'NA':
        model_report['score_color'] = 'grey'
    elif model_report['protchoir_score'] <= 5:
        model_report['score_color'] = 'red'
    elif 5 < model_report['protchoir_score'] <= 7:
        model_report['score_color'] = 'orange'
    elif model_report['protchoir_score'] > 7:
        model_report['score_color'] = 'green'

    pickle.dump(model_report, open(model_oligomer_name+'_CHOIR_model_report.pickle', 'wb'))

    return model_report, '\n'.join(output)




def analyse_oligomers(input_file, template_hitchain, oligomers_list, interfaces_dict, tmdata, report, args, entropies=None, z_entropies=None, minx=None, maxx=None):
    global g_template_hitchain
    global g_interfaces_dict
    global g_tmdata
    global g_report
    global g_args
    global g_entropies
    global g_z_entropies
    global g_minx
    global g_maxx
    global template
    global template_file
    global template_molprobity
    g_template_hitchain = template_hitchain
    g_interfaces_dict = interfaces_dict
    g_tmdata = tmdata
    g_report = report
    g_args = args
    g_entropies = entropies
    g_z_entropies = z_entropies
    g_minx = minx
    g_maxx = maxx
    pctools.print_section(3, 'OLIGOMER ANALYSIS')
    # Define template for comparisons
    template = template_hitchain.split(':')[0]
    template_file = template+'_CHOIR_RelevantChains.pdb'
    reports = []
    if 'M' in args.assessment:
        template_molprobity, molprobity_output = pctools.run_molprobity(template_file, args)
        print(molprobity_output)

    # Run the analysis for all models in parallel
    if args.multiprocess is True:
        p = Pool()
        for model_report, output in p.map_async(analyse_model, oligomers_list).get():
            print(output)
            reports.append(model_report)
        p.close()
        p.join()

    else:
        for oligomer in oligomers_list:
            model_report, output = analyse_model(oligomer)
            print(output)
            reports.append(model_report)

    return reports
