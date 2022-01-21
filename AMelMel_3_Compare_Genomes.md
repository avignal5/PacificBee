# Comparing AMelMel with HAv3.1



<!-- MDTOC maxdepth:6 firsth1:0 numbering:0 flatten:0 bullets:1 updateOnSave:1 -->

- [Alignment of AMelMel on HAv3.1](#alignment-of-amelmel-on-hav31)   
   - [Prepare the index files for the HAv3.1 genome assembly as target sequence](#prepare-the-index-files-for-the-hav31-genome-assembly-as-target-sequence)   
   - [Define an optimal scoring matrix](#define-an-optimal-scoring-matrix)   
   - [Align the query to the target_end](#align-the-query-to-the-target_end)   
   - [Convert to psl file](#convert-to-psl-file)   
   - [make the whole chromosome plots](#make-the-whole-chromosome-plots)   

<!-- /MDTOC -->

## Alignment of AMelMel on HAv3.1
* Done with the LAST software : https://gitlab.com/mcfrith/last
### Prepare the index files for the HAv3.1 genome assembly as target sequence
* Option -uNEAR, as the two sequences are very similar

```bash
lastdb  -P0 -uNEAR -R01 Apis_mellifera_HAV3_1_NEAR  \
        ~/CompareWithHav3/HAV3_1_lastdb/GCF_003254395.2_Amel_HAv3.1_genomic.fna
```
### Define an optimal scoring matrix

```bash
last-train -P4 --revsym --matsym --gapsym -E0.05 -C2 \
            ~/CompareWithHav3/HAV3_1_lastdb/Apis_mellifera_HAV3_1_NEAR \
            ~/CompareWithHav3/HAV3_1_lastdb/GCA_003314205.1_INRA_AMelMel_1.0_genomic.fna \
            > HAV3-AMelMel.mat
```
### Align the query to the target_end

```bash
lastal -m50 -E0.05 -C2 -p HAV3-AMelMel.mat Apis_mellifera_HAV3_1_NEAR \
                        GCA_003314205.1_INRA_AMelMel_1.0_genomic.fna \
                        > HAV3_1_AMelMel.maf
```
### Convert to psl file

```bash
maf-convert  psl HAV3_1_AMelMel.maf > HAV3_1_AMelMelb.psl
```

### make the whole chromosome plots

MAKE THE FOLLOWING SCRIPT MORE GENERIC

```python
#!/usr/bin/env python3
# _*_ coding: Utf-8 _*_
# coding: utf-8

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import collections  as mc
from matplotlib import colors
import re
import csv
import subprocess
from collections import defaultdict

with open("/Users/avignal/Documents/Stats/2016_PacificBee/AgpAndConversionFiles/Ours/GenBankHAV3_1_AMelMel_V2c.txt") as csvFile:
	data=csv.reader(csvFile, delimiter = "\t")
	for row in data:
		if re.search('^LG',row[0]):
			LG = row[0]
			HAV3 = row[1]
			AMelMel = row[2]
			xName = row[4]
			yName = row[3]
			fileName = "_".join(["HAV3_1_AMelMelV2",LG,".pdf"])
			path = "/Users/avignal/Documents/Stats/2016_PacificBee/2019AlignAMelMelHav3_1/WholeChromosomesV2c/"
			#print("\t".join([GenbankName, posBreak, str(start), str(end)]))


			subprocess.run(['/Users/avignal/Documents/Stats/2016_PacificBee/2019AlignAMelMelHav3/mafTotab13.py',
							'-i',
							'/Users/avignal/Documents/Stats/2016_PacificBee/2019AlignAMelMelHav3_1/WholeChromosomesV2c/HAV3_1_AMelMel.psl',
							'-c',
							HAV3,
							'-q',
							AMelMel,
							'-v',
							'/Users/avignal/Documents/Stats/2016_PacificBee/AgpAndConversionFiles/Ours/meanBkptsAMelMelV2c',
							'-z',
							'/Users/avignal/Documents/Stats/2016_PacificBee/AgpAndConversionFiles/FromNCBI_HAV3_1/meanBkptsHAv3_1',
							'-o',
							path + fileName,
							'-w',
							'True',
							'-x',
							xName,
							'-y',
							yName])

```
