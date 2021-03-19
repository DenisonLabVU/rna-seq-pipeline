import argparse
import pandas as pd
import os
import fnmatch
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--Sample_list", help="A tab-delineated file with each sample base name on a new line.")
parser.add_argument("--Virus", help="Virus name. Options are MHV, MERS, SARS2.")
parser.add_argument("--Working_Dir", help="Absolute or relative path of directory with data.")
parser.add_argument("--Output_Dir", help="Absolute or relative path of directory for output folders and files. Default working directory.")
pasrser.add_arugment("--Experiment_Name", help="Experiment name for naming output reports.")
args = parser.parse_args()

#Make a report dataframe with sample column loaded
sample_file = str(args.Sample_List)
report = pd.DataFrame()
report['sample'] = (line.rstrip() for line in sample_file)
#Set other variables
virus = str(args.Virus)
wd = str(args.Working_Dir)
if args.Output_Dir:
    od = str(args.Output_Dir)
else:
    od = wd
exp = str(args.Experiment_Name)

#report['sample'] = ['DMSO-A', 'DMSO-B', 'DMSO-C', 'NHC-2A', 'NHC-2B', 'NHC-2C', 'NHC-4A', 'NHC-4B', 'NHC-4C', 'RDV-0125A', 'RDV-0125B', 'RDV-0125C', 'RDV-05A', 'RDV-05B', 'RDV-05C']
#
if not os.path.exists(od + '/Junction_Files'):
    os.makedirs(od + '/Junction_Files')
save_dir_file = od + 'Junction_Files/'
if not os.path.exists(od + '/Junction_Plots'):
    os.makedirs(od + '/Junction_Plots')
save_dir_plot = od + "Junction_Plots/"

for file in os.listdir(wd):
    if fnmatch.fnmatch(file, "*_virema_Virus_Recombination_Results.bed"):
        sample_name = str(file.split("_")[0])
        bed = pd.read_csv(wd + file, sep="\t", header=0, index_col=False, names=['genome', 'start', 'stop', 'type', 'depth', 'rgb1', 'rgb2', 'start_seq', 'stop_seq'])
        recombined_nts = bed['depth'].sum()
        bed = bed.sort_values(by=['depth'], ascending=True)
        total = bed['depth'].sum()
        bed['frequency'] = bed['depth'] / total
        bed['logfreq'] = np.log10(bed['frequency'])
        bed = bed.reset_index(drop=True)
        bed_forward = bed.loc[bed['start'] < bed['stop']]
        bed_forward = bed_forward.reset_index(drop=True)
        report.loc[report['sample'] == str(sample_name), ['recombined_nts']] = recombined_nts
        bed.to_csv(save_dir_file + sample_name + '_junctions.txt', sep='\t', index=False)
        bed_forward.to_csv(save_dir_file + sample_name + '_forward_junctions.txt', sep='\t', index=False)
        sns.set_style("ticks")
        fig = plt.figure(figsize=(4, 4))
        plt.ioff()
        plt.scatter(bed_forward.stop, bed_forward.start, c=bed_forward.logfreq, cmap='gist_rainbow', alpha=1, vmin=0, vmax=-5, s=15)
        plt.xlim([-1500, 33500])
        plt.ylim([-1500, 33500])
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.xlabel("3' Positon", fontsize=14)
        plt.ylabel("5' Position", fontsize=14)
        cax = fig.add_axes([0.15, 0.95, 0.70, 0.02])
        cbar = plt.colorbar(orientation="horizontal", cax=cax)
        cbar.ax.tick_params(labelsize=10)
        cbar.ax.set_title("log10(Frequency)", fontsize=12)
        plt.savefig(save_dir + sample_name + "_junctionplot.png", dpi=600, bbox_inches='tight')
        plt.close('all')
    if fnmatch.fnmatch(file, "*_virema_Virus_cuttingsites.f.bedgraph"):
        sample_name = str(file.split("_")[0])
        depth_f = pd.read_csv(wd + file, sep="\t", header = 0, index_col=False, names=['genome', 'position', 'position1', 'coverage'])
        f_nts = depth_f['coverage'].sum()
        report.loc[report['sample'] == sample_name, ['total_cutting_f_nts']] = f_nts
    if fnmatch.fnmatch(file, "*_virema_Virus_cuttingsites.r.bedgraph"):
        sample_name = str(file.split("_")[0])
        depth_r = pd.read_csv(wd + file, sep="\t", header=0, index_col=False, names=['genome', 'position', 'position1', 'coverage'])
        r_nts = depth_r['coverage'].sum()
        report.loc[report['sample'] == sample_name, ['total_cutting_r_nts']] = r_nts
    #if fnmatch.fnmatch(file, "_bbmap_coverage.txt"):
        #sample_name = file.split("_")[1]
        #depth = pd.read_csv(file, header =True)
        #total_depth = depth.coverage.sum()
        #report.loc[report['sample'] == sample_name, ['total_nts'] = total_depth
report['cutting_f_jfreq'] = (report['recombined_nts'] / report['total_cutting_f_nts']) * 1000000
report['cutting_r_jfreq'] = (report['recombined_nts'] / report['total_cutting_r_nts']) * 1000000
report.to_csv(wd + exp + "ViReMa_report.txt", sep="\t", index=False)