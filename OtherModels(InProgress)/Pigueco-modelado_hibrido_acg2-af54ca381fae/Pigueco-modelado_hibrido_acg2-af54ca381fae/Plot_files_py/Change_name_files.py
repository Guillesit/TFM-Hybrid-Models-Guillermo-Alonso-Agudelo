import os
foldername='modified/MetodoSimplificadoCeros6/OutputValuesPopulation/'
for filename in os.listdir(foldername):
    for j in range(10,100):
        if filename=="Out%s.txt"%j:
            new='Out0%s.txt'%j
            
            os.rename(foldername+filename,foldername+new )
    for j in range(0,10):
        if filename=="Out%s.txt"%j:
            new='Out00%s.txt'%j
            os.rename(foldername+filename,foldername+new )