#!/bin/python3
##Last modified 3/9/21 by JGB
import argparse
import pandas as pd
import os
import fnmatch
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("Sample_list", help="A text file with each sample base name on a new line.")
parser.add_argument("Virus", help="Virus name. Options are MHV, MERS, SARS2.")
parser.add_argument("Working_Dir", help="Absolute or relative path of directory with data.")
parser.add_argument("Experiment_Name", help="Experiment name for naming output reports.")
parser.add_argument("--Output_Dir", help="Absolute or relative path of directory for output folders and files. Default working directory.")
parser.add_argument("--Shannon_Entropy", help="Path to folder with Virus_Recombination_Results.txt files for Shannon Entropy")
parser.add_argument("--Virus_Accession", help="NCBI virus accession number")
parser.add_argument("--Min_Coverage", help="Minimum counts to include in calculation of Shannon Entropy")
args = parser.parse_args()

#Make a report dataframe with sample column loaded
sample_list = [line.rstrip('\n') for line in open(str(args.Sample_list))]
report = pd.DataFrame()
report['sample'] = sample_list
#Set other variables
virus = str(args.Virus)
wd = str(args.Working_Dir)
if args.Output_Dir:
    od = str(args.Output_Dir)
else:
    od = wd
exp = str(args.Experiment_Name)

#report['sample'] = ['DMSO-A', 'DMSO-B', 'DMSO-C', 'NHC-2A', 'NHC-2B', 'NHC-2C', 'NHC-4A', 'NHC-4B', 'NHC-4C', 'RDV-0125A', 'RDV-0125B', 'RDV-0125C', 'RDV-05A', 'RDV-05B', 'RDV-05C']
#Make target folders for output files
if not os.path.exists(od + '/Junction_Files'):
    os.makedirs(od + '/Junction_Files')
save_dir_file = od + 'Junction_Files/'
if not os.path.exists(od + '/Junction_Plots'):
    os.makedirs(od + '/Junction_Plots')
save_dir_plot = od + "Junction_Plots/"

##Shannon Entropy script originally authored by Andrew Routh.
if args.Shannon_Entropy:
    if args.Min_Coverage:
        Min_Coverage = int(args.Min_Coverage)
    else:
        Min_Coverage = 0
    se_dir = str(args.Shannon_Entropy)
    Virus_Accession = str(args.Virus_Accession)
    se_output_normalized = pd.DataFrame(columns=['sample'])
    se_output_normalized['sample'] = sample_list
    for file in os.listdir(se_dir):
        if fnmatch.fnmatch(file, "*_virema_Virus_Recombination_Results.txt"):
            sample_name = str(file.split("_")[0])
            Dicts = {}
            with open(file, 'r') as file1:
                Data = file1.readline()
                while Data:
                    Name = Data[13:-1]
                    Dicts[Name] = file1.readline().split("\t")[:-1]
                    Data = file1.readline()
                    Data = file1.readline()
            DictKeys = {}
            n = 1
            for Gene in Dicts:
                Data = Dicts[Gene]
                if args.Virus_Accession in Gene:
                    coverage_file = pd.read(wd + sample_name + "_bbmap_coverage.txt", sep="\t", header=0)
                    Virus_Coverage = np.mean(coverage_file['Coverage'])
                    Total_Reads = Virus_Coverage
                else:
                    Total_Reads = 0
                    print("Running Shannon Entropy Calculation for " + sample_name + ". Unknown genome and not normalizing to coverage.")
                Sums = []
                Rec_Total = 0
                for i in Data:
                    data = i.split("_")
                    Freq = int(data[-1])
                    if Freq > Min_Coverage:
                        Rec_Total += Freq
                        Sums.append(Freq)
                    else:
                        pass
                Entropy=0
                # for i in Sums:
                #     Fraction = i / float(Rec_Total)
                #     Entropy -= math.log(Fraction, 2) * Fraction
                #     se_output[Gene] = Gene
                #     se_output[shannon_entropy] = Entropy
                # Entropy = 0
                for i in Sums:
                    Fraction = i / float(Rec_Total + Total_Reads)
                    Entropy -= math.log(Fraction, 2) * Fraction
                Fraction = Total_Reads / float(Rec_Total + Total_Reads)
                Entropy -= math.log(Fraction, 2) * Fraction
                se_output_normalized.loc[se_output_normalized["sample"] == sample_name, [str(Gene)]] = Entropy
se_output_normalized.to_csv(od + exp + "_shannon_entropy_normalized.txt", sep="\t", index=False)
#Isolate forward junctions and make junction plots.
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
        if (virus == 'MHV'):
            sns.set_style("ticks")
            fig = plt.figure(figsize=(4, 4))
            plt.ioff()
            plt.scatter(bed_forward.stop, bed_forward.start, c=bed_forward.logfreq, cmap='gist_rainbow', alpha=1, vmin=0, vmax=-6, s=15)
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
            plt.savefig(save_dir_plot + sample_name + "_junctionplot.png", dpi=600, bbox_inches='tight')
            plt.close('all')
        if (virus == 'MERS' or virus == 'SARS2'):
            sns.set_style("ticks")
            fig = plt.figure(figsize=(4, 4))
            plt.ioff()
            plt.scatter(bed_forward.stop, bed_forward.start, c=bed_forward.logfreq, cmap='gist_rainbow', alpha=1,
                        vmin=0, vmax=-6, s=15)
            plt.xlim([-500, 31500])
            plt.ylim([-500, 31500])
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)
            plt.xlabel("3' Positon", fontsize=14)
            plt.ylabel("5' Position", fontsize=14)
            cax = fig.add_axes([0.15, 0.95, 0.70, 0.02])
            cbar = plt.colorbar(orientation="horizontal", cax=cax)
            cbar.ax.tick_params(labelsize=10)
            cbar.ax.set_title("log10(Frequency)", fontsize=12)
            plt.savefig(save_dir_plot + sample_name + "_junctionplot.png", dpi=600, bbox_inches='tight')
            plt.close('all')
    if fnmatch.fnmatch(file, "*_virema_Virus_cuttingsites.f.bedgraph"):
        sample_name = str(file.split("_")[0])
        depth_f = pd.read_csv(wd + file, sep="\t", header = 0, index_col=False, names=['genome', 'position', 'position1', 'coverage'])
        f_nts = sum(depth_f['coverage'])
        report.loc[report['sample'] == sample_name, ['total_cutting_f_nts']] = f_nts
    if fnmatch.fnmatch(file, "*_virema_Virus_cuttingsites.r.bedgraph"):
        sample_name = str(file.split("_")[0])
        depth_r = pd.read_csv(wd + file, sep="\t", header=0, index_col=False, names=['genome', 'position', 'position1', 'coverage'])
        r_nts = sum(depth_r['coverage'])
        report.loc[report['sample'] == sample_name, ['total_cutting_r_nts']] = r_nts
    if fnmatch.fnmatch(file, "*_coverage.txt"):
        sample_name = file.split("_")[0]
        depth = pd.read_csv(wd + file, sep="\t", header =0)
        total_depth = sum(depth['Coverage'])
        report.loc[report['sample'] == sample_name, ['total_nts']] = total_depth
report['total_cutting_site_nts'] = report['total_cutting_f_nts'] + report['total_cutting_r_nts']
report['cutting_f_jfreq'] = (report['recombined_nts'] / report['total_cutting_f_nts']) * 1000000
report['cutting_r_jfreq'] = (report['recombined_nts'] / report['total_cutting_r_nts']) * 1000000
report['cutting_jfreq'] = (report['recombined_nts'] / report['total_cutting_site_nts']) * 1000000
report['jfreq'] = (report['recombined_nts'] / report['total_nts']) * 1000000
report.to_csv(od + exp + "_ViReMa_report.txt", sep="\t", index=False)