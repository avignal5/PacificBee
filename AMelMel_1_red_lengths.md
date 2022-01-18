# Estimate read lengths





```
#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

import pandas as pd
import numpy as np
import os
import subprocess
import glob
#from matplotlib import pyplot as plt

fastqList = glob.glob('/home/gencel/vignal/work/PacificBeeAll/Run*/Analyse_H5toFastq*/*.gz')
for i in fastqList:
    path = i.split('/')
    run_smrt_lib = path[6].split('_')
    run = run_smrt_lib[1]
    smrt_lib = run_smrt_lib[2].split('.')
    smrt = smrt_lib[1]
    lib = path[8].split('.')[0]
    cmd = "zcat " + i + " | awk '{if(NR%4==2) print length($1)}' > /home/gencel/vignal/work/PacificBeeAll/ReadLength/ReadLength_" + run + "_" + smrt + "_" + lib + ".txt"
    #print(cmd)
    os.system(cmd)
```

