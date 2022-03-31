# Comparing AMelMel with HAv3.1



<!-- MDTOC maxdepth:6 firsth1:1 numbering:0 flatten:0 bullets:1 updateOnSave:1 -->

- [Comparing AMelMel with HAv3.1](#comparing-amelmel-with-hav31)   
   - [Alignment of AMelMel on HAv3.1](#alignment-of-amelmel-on-hav31)   
      - [Prepare the index files for the HAv3.1 genome assembly as target sequence](#prepare-the-index-files-for-the-hav31-genome-assembly-as-target-sequence)   
      - [Define an optimal scoring matrix](#define-an-optimal-scoring-matrix)   
      - [Align the query to the target_end](#align-the-query-to-the-target_end)   
      - [Convert to psl file](#convert-to-psl-file)   
      - [make the whole chromosome plots](#make-the-whole-chromosome-plots)   
         - [Script for passing the parameters for all chromosomes to the plotting script PloAlignments.py](#script-for-passing-the-parameters-for-all-chromosomes-to-the-plotting-script-ploalignmentspy)   
      - [Make liftover gff and gtf annotation files](#make-liftover-gff-and-gtf-annotation-files)   
         - [Remove alignments due to repeats, to obtain a 1 to 1 alignment in chain format](#remove-alignments-due-to-repeats-to-obtain-a-1-to-1-alignment-in-chain-format)   

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
maf-convert  psl HAV3_1_AMelMel.maf > HAV3_1_AMelMel.psl
```

### make the whole chromosome plots

#### Script for passing the parameters for all chromosomes to the plotting script PloAlignments.py

* See plots in /Plot_Chrom_Alignments/Figures/HAv3_1_AMelMel_Chrs/


```python
#!/usr/bin/env python3
# _*_ coding: Utf-8 _*_
# coding: utf-8

#PlotAllChrs.py

"""
For plotting all chromosomes with correct parameters.
Edit paths and file names as needed.
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import collections  as mc
from matplotlib import colors
import re
import csv
import subprocess
from collections import defaultdict

with open("../Data/GenBankHAV3_1_AMelMel.txt") as csvFile:
	data=csv.reader(csvFile, delimiter = "\t")
	for row in data:
		if re.search('^LG',row[0]):
			LG = row[0]
			HAV3 = row[1]
			AMelMel = row[2]
			xName = row[4]
			yName = row[3]
			fileName = "_".join(["HAV3_1_AMelMelV2",LG,".pdf"])
			path = "../Figures/HAv3_1_AMelMel_Chrs/HAv3_1_AMelMel_Chrs/"
			subprocess.run(['./PlotAlignments.py',
							'-i',
							'../Data/HAV3_1_AMelMel.psl',
							'-c',
							HAV3,
							'-q',
							AMelMel,
							'-v',
							'../Data/ContigLimitsAMelMel.txt',
							'-z',
							'../Data/ContigLimitsHAv3_1.txt',
							'-o',
							path + fileName,
							'-w',
							'True',
							'-x',
							xName,
							'-y',
							yName])

```



### Make liftover gff and gtf annotation files

AMelMel is nor annotated by the NCBI, so we make here gtf and gff annotation files with coordinates for AMelMel1.1, using the HAv3.1 annotation files.

#### Remove alignments due to repeats, to obtain a 1 to 1 alignment in chain format

* LAST split 1:

```bash
last-split -m1 HAV3_1_AMelMel_v2c.maf > HAV3_1_AMelMel_v2c_1_.maf
```

* LAST split 2:

```bash
maf-swap HAV3_1_AMelMel_v2c_1_.maf | last-split -m1 > HAV3_1_AMelMel_v2c_2_.maf
```

* Swap HAv3.1 back as the reference:

```bash
maf-swap HAV3_1_AMelMel_v2c_2_.maf > HAV3_1_AMelMel_v2c_2Swapped.maf
```

* Convert to the UCSC .chain format:

```bash
maf-convert  chain HAV3_1_AMelMel_v2c_2Swapped.maf > HAV3_1_AMelMel_v2c_2Swapped.chain
```

* Downloaded the NCBI GCF_003254395.2_Amel_HAv3.1_genomic.gff file from the Apis mellifera reference genome webpage: https://www.ncbi.nlm.nih.gov/genome/?term=apis+mellifera
* Downloaded the NCBI genomic.gtf annotation file from the ncbi_dataset folder: https://www.ncbi.nlm.nih.gov/datasets/genomes/?taxon=7460



* Convert to AMelMel1.1 coordinates with the CrossMap software: http://crossmap.sourceforge.net/

```bash
CrossMap.py gff HAV3_1_AMelMel_v2c_2Swapped.chain GCF_003254395.2_Amel_HAv3.1_genomic.gff Amel_AMelMel1.1_genomic.gff
CrossMap.py gff HAV3_1_AMelMel_v2c_2Swapped.chain genomic.gtf Amel_AMelMel1.1_genomic.gtf
```

* For each annotation file, there is one output with the successfully mapped annotations and one output with the unmapped annotations:

```
Amel_AMelMel1.1_genomic.gff
Amel_AMelMel1.1_genomic.gff.ZKLHTWCY.unmap
Amel_AMelMel1.1_genomic.gtf
Amel_AMelMel1.1_genomic.gtf.JUF2VOH4.unmap

```
