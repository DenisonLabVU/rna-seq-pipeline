#uhomologyplot
import numpy as np
import matplotlib.pyplot as plt

Files = '''MERS_LCP1_virema_MHcompileVirus_Recombination_Results.txt
MERS_LDP1_virema_MHcompileVirus_Recombination_Results.txt
MERS_LEP1_virema_MHcompileVirus_Recombination_Results.txt
SARSCoV2_A_virema_MHcompileVirus_Recombination_Results.txt
SARSCoV2_B_virema_MHcompileVirus_Recombination_Results.txt
SARSCoV2_C_virema_MHcompileVirus_Recombination_Results.txt
WT16hmA_virema_MHcompileVirus_Recombination_Results.txt
WT16hmB_virema_MHcompileVirus_Recombination_Results.txt
WT16hmC_virema_MHcompileVirus_Recombination_Results.txt
XN24hmA_virema_MHcompileVirus_Recombination_Results.txt
XN24hmB_virema_MHcompileVirus_Recombination_Results.txt
XN24hmC_virema_MHcompileVirus_Recombination_Results.txt
WT-1_virema_MHcompileVirus_Recombination_Results.txt
WT-2_virema_MHcompileVirus_Recombination_Results.txt
WT-3_virema_MHcompileVirus_Recombination_Results.txt
XN-1_virema_MHcompileVirus_Recombination_Results.txt
XN-2_virema_MHcompileVirus_Recombination_Results.txt
XN-3_virema_MHcompileVirus_Recombination_Results.txt
'''.split()


Dict = {}
N = 20

#MERS = 'JX869059.2_to_JX869059.2'
#SARS =  'MT020881.1_to_MT020881.1'
#MHV = 'AY910861.1_to_AY910861.1'

for i in Files:
    Dict[i] = np.array([0]*N)
    with open(i, 'r') as In:
        Data = In.readline()
        while Data:
            if 'RevStrand' in Data:
                Data = In.readline()
                Data = In.readline()
                Data = In.readline()
            else:
                Data = In.readline()
                Data = Data.split()
                for j in Data:
                    Fuzz = int(j.split('_')[1][1:])
                    #Count = int(j.split('_')[-1])
                    Count = 1
                    Dict[i][Fuzz] += Count
                Data = In.readline()
                Data = In.readline()
                print(i)
                print(Dict[i])
                Dict[i] = Dict[i]/np.sum(Dict[i])
                print(Dict[i])

def MakeTheoreticalDistribution(N):
    Dist = [0] * (N)
    Dist[0] = 1.0
    for i in range(1,N):
        Prob = 0.25**i
        Dist[i] = Prob
    for i in range(N)[::-1]:
        Dist[i] -= sum(Dist[i+1:])
    return Dist


#Dist = np.array(MakeTheoreticalDistribution(N))
#plt.bar(np.arange(N), Dist, width=0.2, color='b')
#x = 0
#for i in Files:
#    if 'MERS' in i:
#        x += 0.25
#        plt.bar(np.arange(N) + x, Dict[i], width=0.2, color='orange')
#plt.bar(np.arange(N), Dist, width=0.2, color='b')
#x = 0
#for i in Files:
#    if 'SARS' in i:
#        x += 0.25
#        plt.bar(np.arange(N) + x, Dict[i], width=0.2, color='red')
#plt.bar(np.arange(N), Dist, width=0.2, color='b')
#plt.bar(np.arange(N), Dist, width=0.2, color='b')
#x = 0
#for i in Files:
#    if 'WT' in i:
#        x += 0.25
#        plt.bar(np.arange(N) + x, Dict[i], width=0.2, color='green')
#plt.bar(np.arange(N), Dist, width=0.2, color='b')
#x = 0
#for i in Files:
#    if 'XN' in i:
#        x += 0.25
#        plt.bar(np.arange(N) + x, Dict[i], width=0.2, color='cyan')

