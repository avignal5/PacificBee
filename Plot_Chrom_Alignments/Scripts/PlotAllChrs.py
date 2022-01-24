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
			path = "../Figures/HAv3_1_AMelMel_Chrs/"
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
