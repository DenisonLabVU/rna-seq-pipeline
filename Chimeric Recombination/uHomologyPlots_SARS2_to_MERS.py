#uhomologyplot
import numpy as np
import matplotlib.pyplot as plt

Files = '''/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/Shannon_Entropy/0-1A_SARS2_to_MERS_combined.txt
/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/Shannon_Entropy/0-1B_SARS2_to_MERS_combined.txt
/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/Shannon_Entropy/0-1C_SARS2_to_MERS_combined.txt
/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/Shannon_Entropy/0-01A_SARS2_to_MERS_combined.txt
/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/Shannon_Entropy/0-01B_SARS2_to_MERS_combined.txt
/Users/gribblj/Dropbox/Perlman_MERS_SARS2/RNAseq/Shannon_Entropy/0-01C_SARS2_to_MERS_combined.txt
'''.split()


Dict = {}
N = 20

#MERS_to_SARS2 = 'JX869059.2_to_MT020881.1'
SARS2_to_MERS =  'MT020881.1_to_JX869059.2'

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
    print(Dist)


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

