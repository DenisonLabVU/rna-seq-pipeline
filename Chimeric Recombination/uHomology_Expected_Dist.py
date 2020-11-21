import numpy as np

Files = '''
'''.split()

Dict = {}
N=20

MERS_to_SARS2 = 'JX869059.2_to_MT020881.1'
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